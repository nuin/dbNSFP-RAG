"""Shared chromosome streaming and durable checkpoints for database builders."""

import json
import logging
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from .chunker import chunk_dataframe
from .ingest import (
    filter_by_genes,
    list_chromosome_files,
    list_chromosome_files_in_zip,
    parse_chromosome_file,
    parse_chromosome_from_zip,
)

logger = logging.getLogger(__name__)


@dataclass
class BuildOptions:
    source: Path
    chromosomes: list[str] | None = None
    limit: int | None = None
    resume: bool = True
    use_grch37: bool = False
    genes: set[str] | None = None
    gene_data: dict | None = None


def load_checkpoint(checkpoint_path: Path) -> dict:
    if checkpoint_path.exists():
        with checkpoint_path.open() as stream:
            return json.load(stream)
    return {"completed_chromosomes": [], "total_variants": 0}


def save_checkpoint(checkpoint_path: Path, data: dict):
    temporary = checkpoint_path.with_suffix(".json.tmp")
    with temporary.open("w") as stream:
        json.dump(data, stream)
    temporary.replace(checkpoint_path)


def _is_zip(source):
    return source.suffix == ".zip" or (source.is_file() and zipfile.is_zipfile(source))


def source_files(options: BuildOptions) -> list[str]:
    files = (
        list_chromosome_files_in_zip(options.source)
        if _is_zip(options.source)
        else [str(path) for path in list_chromosome_files(options.source)]
    )
    if options.chromosomes:
        chromosomes = {str(c).removeprefix("chr") for c in options.chromosomes}
        files = [f for f in files if any(f"chr{c}.gz" in f or f"chr{c}." in f for c in chromosomes)]
    return files


def chromosome_records(options: BuildOptions, filename: str):
    chunks = (
        parse_chromosome_from_zip(options.source, filename)
        if _is_zip(options.source)
        else parse_chromosome_file(Path(filename))
    )
    for chunk in tqdm(chunks, desc=Path(filename).name, unit="batch"):
        if options.genes is not None:
            chunk = filter_by_genes(chunk, options.genes)
        yield from chunk_dataframe(
            chunk, use_grch37=options.use_grch37, gene_data=options.gene_data
        )


class BuildCheckpoint:
    def __init__(self, store, resume: bool):
        self.store = store
        self.sqlite = hasattr(store, "get_completed_chromosomes")
        self.path = store.db_path / "checkpoint.json" if not self.sqlite else None
        if not resume:
            self.completed = []
        elif self.sqlite:
            self.completed = store.get_completed_chromosomes()
        else:
            self.completed = load_checkpoint(self.path)["completed_chromosomes"]

    def save(self):
        total = self.store.count()
        if self.sqlite:
            self.store.set_completed_chromosomes(self.completed)
            self.store.set_total_variants_processed(total)
        else:
            save_checkpoint(
                self.path, {"completed_chromosomes": self.completed, "total_variants": total}
            )


def _write_batch(store, batch):
    if batch:
        ids, texts, metadata = zip(*batch)
        store.add_variants(list(ids), list(texts), list(metadata), show_progress=False)
        batch.clear()


def _process_chromosome(store, options, filename, seen_ids):
    batch = []
    total = store.count()
    for record in chromosome_records(options, filename):
        if record[0] in seen_ids:
            continue
        if options.limit and total + len(batch) >= options.limit:
            _write_batch(store, batch)
            return False
        seen_ids.add(record[0])
        batch.append(record)
        if len(batch) >= 10000:
            _write_batch(store, batch)
            total = store.count()
    _write_batch(store, batch)
    return True


def run_build(store, options: BuildOptions) -> int:
    if not options.resume:
        store.delete_all()
    checkpoint = BuildCheckpoint(store, options.resume)
    if not options.resume:
        checkpoint.save()
    # FAISS appends records, whereas SQLite replaces the same primary key.
    existing_ids = set(store.metadata["ids"]) if hasattr(store, "metadata") else set()
    try:
        for filename in source_files(options):
            match = re.search(r"chr([0-9XYMNP]+)", Path(filename).name)
            chromosome = match.group(1) if match else Path(filename).name
            if chromosome in checkpoint.completed:
                logger.info("Skipping completed chromosome %s", chromosome)
                continue
            if options.limit and store.count() >= options.limit:
                break
            seen_ids = (
                store.get_variant_ids(chromosome)
                if hasattr(store, "get_variant_ids")
                else existing_ids
            )
            complete = _process_chromosome(store, options, filename, seen_ids)
            if complete:
                checkpoint.completed.append(chromosome)
            checkpoint.save()
            logger.info("Chromosome %s processed; %s stored variants", chromosome, store.count())
            if not complete:
                break
    except KeyboardInterrupt:
        checkpoint.save()
        logger.warning("Build interrupted; completed chromosomes saved. Rerun to resume.")
        return 130
    return 0

"""Main entry point for dbNSFP RAG pipeline."""

import argparse
import json
import sys
from pathlib import Path

from tqdm import tqdm

from .chunker import chunk_dataframe
from .config import RAW_DATA_DIR, VECTORDB_DIR, VECTORDB_GRCH37_NGSGENES, DBNSFP_ZIP, DBNSFP_DIR
from .ingest import (
    list_chromosome_files_in_zip, parse_chromosome_from_zip,
    list_chromosome_files, parse_chromosome_file,
    sort_chr_files, filter_by_genes
)
from .panels import get_panel, list_panels
from .vectorstore import VariantVectorStore


def load_checkpoint(checkpoint_path: Path) -> dict:
    """Load checkpoint data."""
    if checkpoint_path.exists():
        with open(checkpoint_path) as f:
            return json.load(f)
    return {"completed_chromosomes": [], "total_variants": 0}


def save_checkpoint(checkpoint_path: Path, data: dict):
    """Save checkpoint data."""
    with open(checkpoint_path, "w") as f:
        json.dump(data, f)


def build_index(
    source: Path,
    chromosomes: list[str] | None = None,
    limit: int | None = None,
    resume: bool = True,
):
    """
    Build vector index from dbNSFP files with checkpoint/resume support.

    Args:
        source: Path to zip file or directory
        chromosomes: Optional list of chromosomes to process
        limit: Optional limit on total variants (for testing)
        resume: Whether to resume from checkpoint (default True)
    """
    print("=" * 60)
    print("dbNSFP RAG Pipeline - Index Builder")
    print("=" * 60)
    print(f"Source: {source}")

    # Initialize vector store
    print("\nInitializing vector store...")
    store = VariantVectorStore()
    existing_count = store.count()
    print(f"Existing variants in index: {existing_count:,}")

    # Checkpoint file
    checkpoint_path = store.db_path / "checkpoint.json"
    checkpoint = load_checkpoint(checkpoint_path) if resume else {"completed_chromosomes": [], "total_variants": 0}

    if checkpoint["completed_chromosomes"]:
        print(f"Resuming from checkpoint: {len(checkpoint['completed_chromosomes'])} chromosomes done")
        print(f"  Completed: {', '.join(checkpoint['completed_chromosomes'])}")

    # Get chromosome files (sorted by size, smallest first)
    chr_files = list_chromosome_files_in_zip(source)

    # Filter chromosomes if specified
    if chromosomes:
        chr_set = set(str(c).replace("chr", "") for c in chromosomes)
        chr_files = [f for f in chr_files if any(f"chr{c}.gz" in f or f"chr{c}." in f for c in chr_set)]

    print(f"\nWill process {len(chr_files)} chromosome files (smallest first)")

    total_processed = checkpoint["total_variants"]
    batch_size = 10000

    try:
        for chr_file in chr_files:
            # Extract chromosome name
            fname = Path(chr_file).name
            import re
            match = re.search(r"chr([0-9XYMNP]+)", fname)
            chr_name = match.group(1) if match else fname

            # Skip if already completed
            if chr_name in checkpoint["completed_chromosomes"]:
                print(f"\nSkipping {fname} (already completed)")
                continue

            print(f"\n{'='*60}")
            print(f"Processing {fname}")
            print(f"{'='*60}")

            chr_variants = 0
            batch_ids = []
            batch_texts = []
            batch_metas = []
            seen_ids = set()

            chunks = parse_chromosome_from_zip(source, chr_file)
            chunks = tqdm(chunks, desc=f"Processing {fname}", unit="batch")

            for chunk_df in chunks:
                chunk_data = chunk_dataframe(chunk_df)

                for var_id, text, meta in chunk_data:
                    if var_id in seen_ids:
                        continue
                    seen_ids.add(var_id)

                    batch_ids.append(var_id)
                    batch_texts.append(text)
                    batch_metas.append(meta)

                    if len(batch_ids) >= batch_size:
                        store.add_variants(batch_ids, batch_texts, batch_metas, show_progress=False)
                        chr_variants += len(batch_ids)
                        total_processed += len(batch_ids)
                        batch_ids, batch_texts, batch_metas = [], [], []

                        if limit and total_processed >= limit:
                            break

                if limit and total_processed >= limit:
                    break

            # Process remaining batch for this chromosome
            if batch_ids:
                store.add_variants(batch_ids, batch_texts, batch_metas, show_progress=False)
                chr_variants += len(batch_ids)
                total_processed += len(batch_ids)

            # Save checkpoint after each chromosome
            checkpoint["completed_chromosomes"].append(chr_name)
            checkpoint["total_variants"] = total_processed
            save_checkpoint(checkpoint_path, checkpoint)

            print(f"  Chromosome {chr_name}: {chr_variants:,} variants")
            print(f"  Total so far: {total_processed:,} variants")
            print(f"  Checkpoint saved!")

            if limit and total_processed >= limit:
                print(f"\nReached limit of {limit:,} variants")
                break

    except KeyboardInterrupt:
        print("\n\nInterrupted! Progress saved to checkpoint.")
        print(f"Resume with: uv run python -m src.main build")

    print("\n" + "=" * 60)
    print("Indexing Complete!")
    print(f"Total variants processed: {total_processed:,}")
    print(f"Total variants in index: {store.count():,}")
    print(f"Completed chromosomes: {len(checkpoint['completed_chromosomes'])}")
    print("=" * 60)

    return 0


def build_panel_index(
    source: Path,
    panel_name: str,
    use_grch37: bool = True,
    chromosomes: list[str] | None = None,
    limit: int | None = None,
    resume: bool = True,
):
    """
    Build vector index for a specific gene panel with optional GRCh37 coordinates.

    Args:
        source: Path to dbNSFP zip file
        panel_name: Name of gene panel (from panels.py)
        use_grch37: Use GRCh37/hg19 coordinates (default True)
        chromosomes: Optional list of chromosomes to process
        limit: Optional limit on total variants (for testing)
        resume: Whether to resume from checkpoint
    """
    print("=" * 60)
    print(f"dbNSFP RAG Pipeline - Panel Index Builder")
    print("=" * 60)
    print(f"Panel: {panel_name}")
    print(f"Genome build: {'GRCh37 (hg19)' if use_grch37 else 'GRCh38'}")
    print(f"Source: {source}")

    # Get panel genes
    try:
        panel_genes = get_panel(panel_name)
    except ValueError as e:
        print(f"Error: {e}")
        return 1

    print(f"Panel contains {len(panel_genes)} genes")

    # Set database path based on build and panel
    if use_grch37 and panel_name.lower() == "ngsgenes":
        db_path = VECTORDB_GRCH37_NGSGENES
    else:
        db_path = VECTORDB_DIR / f"{panel_name.lower()}-{'grch37' if use_grch37 else 'grch38'}"

    print(f"Database path: {db_path}")

    # Initialize vector store
    print("\nInitializing vector store...")
    store = VariantVectorStore(db_path=db_path)
    existing_count = store.count()
    print(f"Existing variants in index: {existing_count:,}")

    # Checkpoint file
    checkpoint_path = db_path / "checkpoint.json"
    checkpoint = load_checkpoint(checkpoint_path) if resume else {"completed_chromosomes": [], "total_variants": 0}

    if checkpoint["completed_chromosomes"]:
        print(f"Resuming from checkpoint: {len(checkpoint['completed_chromosomes'])} chromosomes done")
        print(f"  Completed: {', '.join(checkpoint['completed_chromosomes'])}")

    # Detect source type and get chromosome files
    import zipfile
    is_zip = source.suffix == ".zip" or (source.is_file() and zipfile.is_zipfile(source))

    if is_zip:
        chr_files = list_chromosome_files_in_zip(source)
    else:
        chr_files = [str(f) for f in list_chromosome_files(source)]

    # Filter chromosomes if specified
    if chromosomes:
        chr_set = set(str(c).replace("chr", "") for c in chromosomes)
        chr_files = [f for f in chr_files if any(f"chr{c}.gz" in f or f"chr{c}." in f for c in chr_set)]

    print(f"\nWill process {len(chr_files)} chromosome files")
    print(f"Source type: {'ZIP' if is_zip else 'Directory'}")

    total_processed = checkpoint["total_variants"]
    total_filtered = 0
    batch_size = 10000

    try:
        for chr_file in chr_files:
            # Extract chromosome name
            fname = Path(chr_file).name
            import re
            match = re.search(r"chr([0-9XYMNP]+)", fname)
            chr_name = match.group(1) if match else fname

            # Skip if already completed
            if chr_name in checkpoint["completed_chromosomes"]:
                print(f"\nSkipping {fname} (already completed)")
                continue

            print(f"\n{'='*60}")
            print(f"Processing {fname}")
            print(f"{'='*60}")

            chr_variants = 0
            chr_filtered = 0
            batch_ids = []
            batch_texts = []
            batch_metas = []
            seen_ids = set()

            if is_zip:
                chunks = parse_chromosome_from_zip(source, chr_file)
            else:
                chunks = parse_chromosome_file(Path(chr_file))
            chunks = tqdm(chunks, desc=f"Processing {fname}", unit="batch")

            for chunk_df in chunks:
                # Filter by panel genes
                original_len = len(chunk_df)
                chunk_df = filter_by_genes(chunk_df, panel_genes)
                chr_filtered += original_len - len(chunk_df)

                if chunk_df.empty:
                    continue

                # Convert to embeddings format with GRCh37 option
                chunk_data = chunk_dataframe(chunk_df, use_grch37=use_grch37)

                for var_id, text, meta in chunk_data:
                    if var_id in seen_ids:
                        continue
                    seen_ids.add(var_id)

                    batch_ids.append(var_id)
                    batch_texts.append(text)
                    batch_metas.append(meta)

                    if len(batch_ids) >= batch_size:
                        store.add_variants(batch_ids, batch_texts, batch_metas, show_progress=False)
                        chr_variants += len(batch_ids)
                        total_processed += len(batch_ids)
                        batch_ids, batch_texts, batch_metas = [], [], []

                        if limit and total_processed >= limit:
                            break

                if limit and total_processed >= limit:
                    break

            # Process remaining batch
            if batch_ids:
                store.add_variants(batch_ids, batch_texts, batch_metas, show_progress=False)
                chr_variants += len(batch_ids)
                total_processed += len(batch_ids)

            total_filtered += chr_filtered

            # Save checkpoint
            checkpoint["completed_chromosomes"].append(chr_name)
            checkpoint["total_variants"] = total_processed
            save_checkpoint(checkpoint_path, checkpoint)

            print(f"  Chromosome {chr_name}: {chr_variants:,} panel variants")
            print(f"  Filtered out: {chr_filtered:,} non-panel variants")
            print(f"  Total so far: {total_processed:,} variants")

            if limit and total_processed >= limit:
                print(f"\nReached limit of {limit:,} variants")
                break

    except KeyboardInterrupt:
        print("\n\nInterrupted! Progress saved to checkpoint.")
        print(f"Resume with: uv run python -m src.main build-panel --panel {panel_name}")

    print("\n" + "=" * 60)
    print("Panel Indexing Complete!")
    print(f"Panel: {panel_name} ({len(panel_genes)} genes)")
    print(f"Total variants indexed: {total_processed:,}")
    print(f"Total filtered out: {total_filtered:,}")
    print(f"Database: {db_path}")
    print("=" * 60)

    return 0


def show_stats():
    """Show index statistics."""
    store = VariantVectorStore()
    count = store.count()

    print("dbNSFP RAG Index Statistics (ClickHouse)")
    print("-" * 40)
    print(f"Total variants indexed: {count:,}")

    if count > 0:
        print("\nSample variants by gene:")
        for gene in ["BRCA1", "BRCA2", "TP53", "EGFR", "KRAS"]:
            results = store.search_by_gene(gene, k=100)
            if results:
                print(f"  {gene}: {len(results)} variants")


def clear_index():
    """Clear the vector index and checkpoint."""
    store = VariantVectorStore()
    count = store.count()
    checkpoint_path = store.db_path / "checkpoint.json"

    if count == 0 and not checkpoint_path.exists():
        print("Index is already empty")
        return 0

    confirm = input(f"Delete {count:,} variants and checkpoint? [y/N]: ")
    if confirm.lower() == "y":
        store.delete_all()
        if checkpoint_path.exists():
            checkpoint_path.unlink()
        print("Index and checkpoint cleared")
    else:
        print("Cancelled")

    return 0


def list_chromosomes(source: Path):
    """List available chromosomes in the source."""
    if source.suffix == ".zip":
        files = list_chromosome_files_in_zip(source)
        print(f"Chromosomes in {source.name}:")
        for f in files:
            print(f"  - {Path(f).name}")
    else:
        print("Listing only supported for zip files")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="dbNSFP RAG Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build full index (auto-resumes if stopped)
  python -m src.main build

  # Start fresh (ignore checkpoint)
  python -m src.main build --fresh

  # Index specific chromosomes only
  python -m src.main build --chr 17 --chr 13

  # Index sample for testing
  python -m src.main build --limit 50000

  # Show index stats
  python -m src.main stats

  # List chromosomes in zip (shows size order)
  python -m src.main list

  # Clear index and checkpoint
  python -m src.main clear
        """,
    )

    # Use directory if available, otherwise zip
    default_source = DBNSFP_DIR if DBNSFP_DIR and DBNSFP_DIR.exists() else DBNSFP_ZIP

    parser.add_argument(
        "--source",
        type=Path,
        default=default_source,
        help="Path to dbNSFP zip or directory",
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Build command
    build_parser = subparsers.add_parser("build", help="Build index")
    build_parser.add_argument(
        "--chr",
        action="append",
        dest="chromosomes",
        help="Chromosome(s) to index (can repeat, e.g., --chr 17 --chr 13)",
    )
    build_parser.add_argument(
        "--limit",
        type=int,
        help="Limit number of variants to process",
    )
    build_parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start fresh, ignoring any checkpoint",
    )

    # Stats command
    subparsers.add_parser("stats", help="Show index statistics")

    # List command
    subparsers.add_parser("list", help="List chromosomes in source")

    # Clear command
    subparsers.add_parser("clear", help="Clear the index")

    # Build-panel command
    panel_parser = subparsers.add_parser("build-panel", help="Build panel-specific index")
    panel_parser.add_argument(
        "--panel", "-p",
        required=True,
        help="Panel name (e.g., NGSgenes, hereditary_cancer)",
    )
    panel_parser.add_argument(
        "--build", "-b",
        choices=["grch37", "grch38"],
        default="grch37",
        help="Genome build (default: grch37)",
    )
    panel_parser.add_argument(
        "--chr",
        action="append",
        dest="chromosomes",
        help="Specific chromosome(s) to process",
    )
    panel_parser.add_argument(
        "--limit",
        type=int,
        help="Limit variants for testing",
    )
    panel_parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start fresh, ignoring checkpoint",
    )

    # List-panels command
    subparsers.add_parser("list-panels", help="List available gene panels")

    args = parser.parse_args()

    if args.command == "build":
        return build_index(
            args.source,
            chromosomes=args.chromosomes,
            limit=args.limit,
            resume=not args.fresh,
        )

    elif args.command == "stats":
        return show_stats()

    elif args.command == "list":
        return list_chromosomes(args.source)

    elif args.command == "clear":
        return clear_index()

    elif args.command == "build-panel":
        return build_panel_index(
            args.source,
            panel_name=args.panel,
            use_grch37=(args.build == "grch37"),
            chromosomes=args.chromosomes,
            limit=args.limit,
            resume=not args.fresh,
        )

    elif args.command == "list-panels":
        print("Available gene panels:")
        for name, count in list_panels().items():
            print(f"  {name}: {count} genes")
        return 0

    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())

"""Build checkpoints must represent durable, complete work."""

from pathlib import Path

import pytest

from src.build_pipeline import BuildOptions, run_build
from src.variantdb import VariantDatabase


def _records(count, chromosome="1"):
    for pos in range(1, count + 1):
        yield (
            f"{chromosome}_{pos}_A_G",
            f"variant {pos}",
            {"chr": chromosome, "hg38_chr": "1", "pos": pos, "ref": "A", "alt": "G"},
        )


@pytest.mark.parametrize("chromosome", ["1", "2"])
def test_limit_does_not_mark_partial_chromosome_complete(tmp_path, monkeypatch, chromosome):
    from src import build_pipeline

    monkeypatch.setattr(build_pipeline, "source_files", lambda options: ["chr1.gz"])
    monkeypatch.setattr(build_pipeline, "chromosome_records", lambda *args: _records(3, chromosome))
    db = VariantDatabase(tmp_path / "test.db")
    try:
        assert run_build(db, BuildOptions(Path("unused"), limit=2)) == 0
        assert db.count() == 2
        assert db.get_completed_chromosomes() == []
        assert run_build(db, BuildOptions(Path("unused"), limit=3)) == 0
        assert db.count() == 3
        assert db.get_completed_chromosomes() == ["1"]
    finally:
        db.close()


def test_interrupt_returns_nonzero_and_leaves_chromosome_resumable(tmp_path, monkeypatch):
    from src import build_pipeline

    def interrupted(*args):
        yield from _records(10001)
        raise KeyboardInterrupt

    monkeypatch.setattr(build_pipeline, "source_files", lambda options: ["chr1.gz"])
    monkeypatch.setattr(build_pipeline, "chromosome_records", interrupted)
    db = VariantDatabase(tmp_path / "test.db")
    try:
        assert run_build(db, BuildOptions(Path("unused"))) == 130
        assert db.get_completed_chromosomes() == []
        assert db.count() == 10000
        assert db.get_total_variants_processed() == db.count()
        monkeypatch.setattr(build_pipeline, "chromosome_records", lambda *args: _records(10002))
        assert run_build(db, BuildOptions(Path("unused"))) == 0
        assert db.count() == 10002
    finally:
        db.close()


def test_sqlite_cli_help_does_not_import_optional_faiss(monkeypatch):
    from src import main

    monkeypatch.setattr("sys.argv", ["dbnsfp", "build-sqlite", "--help"])
    with pytest.raises(SystemExit) as exc:
        main.main()
    assert exc.value.code == 0


def test_fresh_build_replaces_stale_annotations(tmp_path, monkeypatch):
    from src import build_pipeline

    monkeypatch.setattr(build_pipeline, "source_files", lambda options: ["chr1.gz"])
    monkeypatch.setattr(build_pipeline, "chromosome_records", lambda *args: _records(2))
    db = VariantDatabase(tmp_path / "test.db")
    try:
        db.add_variants(["1_1_A_G"], ["stale"], [{"chr": "1", "pos": 1, "ref": "A", "alt": "G"}])
        db.set_completed_chromosomes(["1"])
        assert run_build(db, BuildOptions(Path("unused"), resume=False)) == 0
        assert db.count() == 2
        assert db.get_by_id("1_1_A_G")["document"] == "variant 1"
    finally:
        db.close()


class AppendStore:
    """Exercise the append-only backend contract without loading an embedding model."""

    def __init__(self, directory):
        self.db_path = directory
        self.metadata = {"ids": [], "documents": [], "metas": []}

    def count(self):
        return len(self.metadata["ids"])

    def delete_all(self):
        self.metadata = {"ids": [], "documents": [], "metas": []}

    def add_variants(self, ids, texts, metadata, show_progress=False):
        self.metadata["ids"].extend(ids)
        self.metadata["documents"].extend(texts)
        self.metadata["metas"].extend(metadata)


def test_append_backend_resume_does_not_duplicate_saved_ids(tmp_path, monkeypatch):
    from src import build_pipeline

    monkeypatch.setattr(build_pipeline, "source_files", lambda options: ["chr1.gz"])
    monkeypatch.setattr(build_pipeline, "chromosome_records", lambda *args: _records(3))
    store = AppendStore(tmp_path)
    assert run_build(store, BuildOptions(Path("unused"), limit=2)) == 0
    assert store.count() == 2
    assert run_build(store, BuildOptions(Path("unused"), limit=3)) == 0
    assert store.metadata["ids"] == ["1_1_A_G", "1_2_A_G", "1_3_A_G"]


def test_append_backend_fresh_rebuilds_saved_annotations(tmp_path, monkeypatch):
    from src import build_pipeline

    monkeypatch.setattr(build_pipeline, "source_files", lambda options: ["chr1.gz"])
    monkeypatch.setattr(build_pipeline, "chromosome_records", lambda *args: _records(2))
    store = AppendStore(tmp_path)
    store.add_variants(["1_1_A_G"], ["stale"], [{}])
    assert run_build(store, BuildOptions(Path("unused"), resume=False)) == 0
    assert store.metadata["documents"] == ["variant 1", "variant 2"]


@pytest.mark.parametrize("archive", [False, True])
def test_sqlite_build_streams_annotations_and_populates_panels(tmp_path, monkeypatch, archive):
    import gzip
    import zipfile

    from src import builders

    columns = ["#chr", "pos(1-based)", "hg19_chr", "hg19_pos(1-based)", "ref", "alt", "genename"]
    rows = ["1\t100\t1\t90\tA\tG\tBRCA1", "1\t101\t1\t91\tA\tT\tOTHER"]
    data = gzip.compress(("\t".join(columns) + "\n" + "\n".join(rows) + "\n").encode())
    name = "dbNSFP5.4_variant.chr1.gz"
    source = tmp_path / "input"
    source.mkdir()
    if archive:
        source = source / "input.zip"
        with zipfile.ZipFile(source, "w") as output:
            output.writestr(name, data)
    else:
        (source / name).write_bytes(data)
    monkeypatch.setattr(builders, "load_gene_data", lambda: {})
    monkeypatch.setattr(builders, "get_all_panel_genes", lambda: {"BRCA1"})
    monkeypatch.setattr(builders, "PANELS", {"test": {"BRCA1"}})
    path = tmp_path / "test.db"
    assert builders.build_sqlite_index(source, db_path=path) == 0
    db = VariantDatabase(path)
    try:
        assert db.count() == 1
        assert db.get_by_id("1_90_A_G")["metadata"]["gene"] == "BRCA1"
        assert db.count_by_panel("test") == 1
        assert db.get_completed_chromosomes() == ["1"]
    finally:
        db.close()

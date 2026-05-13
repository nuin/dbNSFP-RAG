#!/usr/bin/env python3
"""Load cp_new per-gene TSVs into the SQLite variant database.

Reads every data/exports/cp_new/{GENE}.tsv produced by cp_new.py (+ optionally
annotated by annotate_spliceai.py), converts each row to (variant_id, document,
metadata) tuples using the same chunker the build-sqlite pipeline uses, and
inserts via VariantDatabase.add_variants. The two new columns
(gnomad_v41_faf95_grpmax, spliceai_ds_max_masked) ride along through
METADATA_COLUMNS.

After loading, registers panel membership via add_panel_membership(panel='cp_new').

Idempotent: re-running overwrites any existing rows with the same variant_id
(INSERT OR REPLACE) and skips duplicate panel rows (INSERT OR IGNORE in the
junction table).

Usage:
  uv run python scripts/load_cp_new_sqlite.py \
      --tsv-dir data/exports/cp_new

  # Use the default GRCh37-all-panels DB (production target):
  # data/sqlite/grch37-all-panels.db
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.chunker import chunk_dataframe
from src.config import SQLITE_DB_PATH
from src.gene_data import load_gene_data
from src.variantdb import VariantDatabase

PANEL_NAME = "cp_new"


def coerce_extras(meta: dict, row: pd.Series) -> None:
    """Move per-row gnomAD FAF and SpliceAI masked values into the metadata dict.

    chunker.variant_to_metadata doesn't know about these columns (they're added
    by cp_new.py + annotate_spliceai.py), so we copy them in by hand.
    """
    for src_col, dst_key in (
        ("gnomad_v41_faf95_grpmax", "gnomad_v41_faf95_grpmax"),
        ("spliceai_ds_max_masked", "spliceai_ds_max_masked"),
    ):
        if src_col not in row.index:
            continue
        val = row[src_col]
        if pd.isna(val) or val == "" or val == ".":
            continue
        try:
            meta[dst_key] = float(val)
        except (TypeError, ValueError):
            pass


def load_tsv(tsv_path: Path, gene_data: dict) -> list[tuple[str, str, dict]]:
    """Read one per-gene TSV, return [(var_id, text, meta), ...] for GRCh37 rows."""
    first = tsv_path.read_text().splitlines()[:1]
    if not first or not first[0].startswith("#chr"):
        return []
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, na_values=[".", ""], low_memory=False)
    if df.empty:
        return []
    # chunk_dataframe builds GRCh37 var_id from hg19_chr / hg19_pos(1-based)
    chunked = chunk_dataframe(df, use_grch37=True, gene_data=gene_data)

    # Need to align chunked rows back to source rows for the extras.
    # chunk_dataframe skips rows with missing hg19_pos -- so iterate the same way.
    extras_iter = iter(df.iterrows())
    out = []
    for var_id, text, meta in chunked:
        while True:
            try:
                _, src_row = next(extras_iter)
            except StopIteration:
                src_row = None
                break
            hg19_pos = src_row.get("hg19_pos(1-based)")
            if pd.isna(hg19_pos) or hg19_pos == ".":
                continue
            break
        if src_row is not None:
            coerce_extras(meta, src_row)
        out.append((var_id, text, meta))
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv-dir", type=Path, default=Path("data/exports/cp_new"),
                   help="Per-gene TSV directory (output of cp_new.py)")
    p.add_argument("--db", type=Path, default=SQLITE_DB_PATH,
                   help="SQLite database path")
    p.add_argument("--dry-run", action="store_true",
                   help="Parse + report counts, do not write to DB")
    args = p.parse_args()

    if not args.tsv_dir.is_dir():
        sys.exit(f"TSV directory not found: {args.tsv_dir}")

    tsvs = sorted(args.tsv_dir.glob("*.tsv"))
    if not tsvs:
        sys.exit(f"No *.tsv files in {args.tsv_dir}")

    print(f"DB: {args.db}")
    print(f"TSV dir: {args.tsv_dir} ({len(tsvs)} files)")
    print(f"Panel: {PANEL_NAME}")
    print(f"Mode: {'dry-run' if args.dry_run else 'WRITE'}")

    print("\nLoading gene constraint data (pLI/LOEUF/mis_z)...")
    gene_data = load_gene_data()
    print(f"  loaded for {len(gene_data):,} genes")

    db: VariantDatabase | None = None if args.dry_run else VariantDatabase(db_path=args.db)
    if db is not None:
        print(f"  existing variants in DB: {db.count():,}")

    total_loaded = 0
    cp_new_ids: list[str] = []

    for tsv in tqdm(tsvs, desc="load"):
        rows = load_tsv(tsv, gene_data)
        if not rows:
            continue
        ids = [r[0] for r in rows]
        texts = [r[1] for r in rows]
        metas = [r[2] for r in rows]
        cp_new_ids.extend(ids)
        if db is not None:
            db.add_variants(ids, texts, metas, show_progress=False)
        total_loaded += len(rows)

    print(f"\nVariant rows processed: {total_loaded:,}")
    print(f"Unique variant_ids accumulated: {len(set(cp_new_ids)):,}")

    if db is not None:
        print(f"Registering panel membership ({PANEL_NAME})...")
        db.add_panel_membership(list(set(cp_new_ids)), PANEL_NAME)
        print(f"DB total variants after load: {db.count():,}")
        # Verify panel membership
        panel_count = db._conn.execute(
            "SELECT COUNT(*) AS c FROM variant_panels WHERE panel_name = ?",
            (PANEL_NAME,),
        ).fetchone()["c"]
        print(f"variant_panels rows for {PANEL_NAME}: {panel_count:,}")
        db.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Build three deliverables from the cp_new pipeline output:

1. SeqNext-minimal TSV (4 columns: gene, transcript, hgvs_c, classification)
       -> data/exports/cp_new/cp_new_seqnext_minimal.tsv
   For direct upload into SeqNext. Optionally excludes flagged rows
   (--drop-flagged) so SeqNext never sees the 70 intergenic ones.

2. Complete annotation TSV (every column from every per-gene file, plus
   the classification + VV resolution)
       -> data/exports/cp_new/cp_new_all_annotations.tsv.gz
   One row per (variant, gene) pairing. ~54,000 rows. Gzipped because it's
   wide (60+ columns) and large enough to matter.

3. Visualization-friendly SQLite database
       -> data/exports/cp_new/cp_new_viz.db
   Two tables:
     - cp_new_variants: every annotated variant, wide layout, indexed on
       (chr_grch37, pos_grch37), gene, classification
     - cp_new_classifications: just the 4,567 classified rows from SeqNext,
       indexed on gene, transcript, vv_status
   Designed to drop into Datasette, DB Browser for SQLite, Metabase, or
   any tool that can read .db. No FAISS / no ML deps, just standard SQL.

Usage:
  uv run python scripts/build_cp_new_exports.py
  uv run python scripts/build_cp_new_exports.py --drop-flagged
"""

from __future__ import annotations

import argparse
import gzip
import sqlite3
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

TSV_DIR = Path("data/exports/cp_new")
SEQNEXT_DIR = TSV_DIR / "seqnext"
SEQNEXT_COMBINED = SEQNEXT_DIR / "cp_new_seqnext.tsv"
SEQNEXT_MINIMAL = TSV_DIR / "cp_new_seqnext_minimal.tsv"
ALL_ANNOTATIONS = TSV_DIR / "cp_new_all_annotations.tsv.gz"
VIZ_DB = TSV_DIR / "cp_new_viz.db"


def build_seqnext_minimal(drop_flagged: bool) -> int:
    """Write the 4-column SeqNext-ready TSV."""
    if not SEQNEXT_COMBINED.exists():
        sys.exit(f"Missing input: {SEQNEXT_COMBINED}")
    df = pd.read_csv(SEQNEXT_COMBINED, sep="\t", dtype=str).fillna("")
    n_before = len(df)
    if drop_flagged and "vv_status" in df.columns:
        df = df[df["vv_status"] == "ok"].copy()
    cols = ["gene", "transcript", "hgvs_c", "classification"]
    df[cols].to_csv(SEQNEXT_MINIMAL, sep="\t", index=False)
    print(f"  SeqNext minimal:  {len(df):>6,} / {n_before:,} rows  -> {SEQNEXT_MINIMAL}")
    return len(df)


def build_all_annotations() -> int:
    """Concatenate every per-gene TSV with classification info merged in."""
    # Load per-gene annotated TSVs
    tsvs = sorted(p for p in TSV_DIR.glob("*.tsv") if p.parent == TSV_DIR)
    if not tsvs:
        sys.exit(f"No per-gene TSVs in {TSV_DIR}")

    # Load classifications keyed on (chr,pos,ref,alt,transcript)
    cls_df = pd.read_csv(SEQNEXT_COMBINED, sep="\t", dtype=str).fillna("")
    cls_df["key"] = (
        cls_df["chr_grch37"] + "_" +
        cls_df["pos_grch37"] + "_" +
        cls_df["ref"] + "_" + cls_df["alt"]
    )
    cls_map = dict(zip(cls_df["key"], zip(cls_df["transcript"], cls_df["hgvs_c"],
                                           cls_df["classification"], cls_df.get("vv_status", ""))))

    all_chunks = []
    for tsv in tqdm(tsvs, desc="all-annotations"):
        first = tsv.read_text().splitlines()[:1]
        if not first or not first[0].startswith("#chr"):
            continue
        df = pd.read_csv(tsv, sep="\t", dtype=str, low_memory=False).fillna("")
        if df.empty:
            continue
        # Compose join key for classification lookup
        def lookup(row, idx):
            key = f"{row['hg19_chr']}_{row['hg19_pos(1-based)']}_{row['ref']}_{row['alt']}"
            v = cls_map.get(key.replace(".0_", "_"))
            if v is None:
                return ""
            return v[idx] if v[idx] is not None else ""
        # Strip trailing .0 if pandas read int as float
        df["hg19_pos(1-based)"] = df["hg19_pos(1-based)"].str.replace(r"\.0$", "", regex=True)
        df["seqnext_transcript"] = df.apply(lambda r: lookup(r, 0), axis=1)
        df["seqnext_hgvs_c"]     = df.apply(lambda r: lookup(r, 1), axis=1)
        df["classification"]     = df.apply(lambda r: lookup(r, 2), axis=1)
        df["vv_status"]          = df.apply(lambda r: lookup(r, 3), axis=1)
        all_chunks.append(df)

    if not all_chunks:
        sys.exit("No annotated chunks to write")

    combined = pd.concat(all_chunks, ignore_index=True)
    with gzip.open(ALL_ANNOTATIONS, "wt") as out:
        combined.to_csv(out, sep="\t", index=False)
    print(f"  All annotations:  {len(combined):>6,} rows, {len(combined.columns)} cols  -> {ALL_ANNOTATIONS}")
    return len(combined)


def build_viz_db() -> int:
    """Standalone SQLite DB with two tables, indexed for query."""
    if not ALL_ANNOTATIONS.exists():
        sys.exit(f"Run --all-annotations first: missing {ALL_ANNOTATIONS}")
    if VIZ_DB.exists():
        VIZ_DB.unlink()

    conn = sqlite3.connect(str(VIZ_DB))
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    # cp_new_variants: every annotated variant row, wide layout
    print("  loading all annotations...")
    with gzip.open(ALL_ANNOTATIONS, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False).fillna("")
    df.columns = [c.replace("#", "").replace("(", "_").replace(")", "").replace("-", "_")
                  for c in df.columns]
    df.to_sql("cp_new_variants", conn, if_exists="replace", index=False)
    n_variants = len(df)

    # cp_new_classifications: just the SeqNext rows (4 columns + coords)
    cls_df = pd.read_csv(SEQNEXT_COMBINED, sep="\t", dtype=str).fillna("")
    cls_df.to_sql("cp_new_classifications", conn, if_exists="replace", index=False)
    n_classified = len(cls_df)

    # Indexes for visualization tools
    conn.executescript("""
        CREATE INDEX idx_variants_gene ON cp_new_variants(genename);
        CREATE INDEX idx_variants_classification ON cp_new_variants(classification);
        CREATE INDEX idx_variants_coords ON cp_new_variants(hg19_chr, hg19_pos_1_based);
        CREATE INDEX idx_classifications_gene ON cp_new_classifications(gene);
        CREATE INDEX idx_classifications_transcript ON cp_new_classifications(transcript);
        CREATE INDEX idx_classifications_vv_status ON cp_new_classifications(vv_status);
    """)
    conn.execute("ANALYZE")
    conn.commit()
    conn.close()

    size_mb = VIZ_DB.stat().st_size / (1024 * 1024)
    print(f"  Visualization DB: {n_variants:>6,} variants, {n_classified:,} classified")
    print(f"                    {size_mb:.1f} MB at {VIZ_DB}")
    return n_variants


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--drop-flagged", action="store_true",
                   help="Exclude vv_status != 'ok' rows from the SeqNext-minimal export")
    p.add_argument("--only", choices=["seqnext", "annotations", "viz"],
                   help="Build just one of the three artifacts")
    args = p.parse_args()

    print(f"=== Building cp_new exports ===\n")

    if args.only in (None, "seqnext"):
        print("1. SeqNext minimal TSV (4 columns)")
        build_seqnext_minimal(args.drop_flagged)
        print()

    if args.only in (None, "annotations"):
        print("2. Complete annotations TSV (gzipped)")
        build_all_annotations()
        print()

    if args.only in (None, "viz"):
        print("3. Visualization SQLite DB")
        build_viz_db()
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())

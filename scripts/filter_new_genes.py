#!/usr/bin/env python3
"""Filter cp_new pipeline outputs to the subset of genes that are NEW in
CP_new (APR2026) vs the old CP (01JUN2021) panel.

Diffs the two BED files for gene symbols, then slices every SeqNext
output (combined, minimal, strict, per-gene files, review files) to only
the new-gene rows. Writes results under a `new_genes_only/` subfolder of
the SeqNext directory and into the bundle.

Usage:
  uv run python scripts/filter_new_genes.py
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

CP_OLD_BED = Path("/Users/nuin/Projects/ahs/new_bed/CP/C+_ALL_IDPE_01JUN2021.bed")
CP_NEW_BED = Path("/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed")

DATA = Path("data/exports/cp_new")
SEQNEXT_DIR = DATA / "seqnext"
OUT_DIR = SEQNEXT_DIR / "new_genes_only"

BUNDLE_OUT = Path("cp_new_bundle/outputs/new_genes_only")


def gene_set(bed: Path) -> set[str]:
    s: set[str] = set()
    with open(bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                s.add(p[4])
    return s


def filter_tsv(src: Path, dst: Path, new_genes: set[str], gene_col: str = "gene") -> int:
    if not src.exists() or src.stat().st_size < 100:
        return 0
    df = pd.read_csv(src, sep="\t", dtype=str, low_memory=False).fillna("")
    if gene_col not in df.columns:
        return 0
    sub = df[df[gene_col].isin(new_genes)]
    dst.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(dst, sep="\t", index=False)
    return len(sub)


def main() -> int:
    old = gene_set(CP_OLD_BED)
    new = gene_set(CP_NEW_BED)
    new_only = sorted(new - old)
    print(f"CP (old):     {len(old)} genes")
    print(f"CP_new:       {len(new)} genes")
    print(f"new in CP_new: {len(new_only)} genes\n")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    BUNDLE_OUT.mkdir(parents=True, exist_ok=True)

    # Write the gene list
    (OUT_DIR / "_new_genes.txt").write_text("\n".join(new_only) + "\n")
    shutil.copy(OUT_DIR / "_new_genes.txt", BUNDLE_OUT / "_new_genes.txt")

    # Filter the SeqNext files
    targets = [
        ("cp_new_seqnext_combined.tsv", SEQNEXT_DIR / "cp_new_seqnext.tsv"),
        ("cp_new_seqnext_minimal.tsv",  DATA / "cp_new_seqnext_minimal.tsv"),
        ("cp_new_seqnext_strict.tsv",   DATA / "cp_new_seqnext_strict.tsv"),
        ("_review_canonical_splice.tsv", SEQNEXT_DIR / "_review_canonical_splice.tsv"),
        ("_review_intergenic.tsv",       SEQNEXT_DIR / "_review_intergenic.tsv"),
    ]

    new_set = set(new_only)
    for out_name, src in targets:
        n = filter_tsv(src, OUT_DIR / out_name, new_set)
        shutil.copy(OUT_DIR / out_name, BUNDLE_OUT / out_name)
        print(f"  {out_name:<38} {n:>6,} rows")

    # Per-gene files: just copy the ones whose stem is a new gene
    per_gene_in = SEQNEXT_DIR
    per_gene_out = OUT_DIR / "per_gene"
    per_gene_out.mkdir(exist_ok=True)
    bundle_per_gene = BUNDLE_OUT / "per_gene"
    bundle_per_gene.mkdir(exist_ok=True)
    copied = 0
    for g in new_only:
        src = per_gene_in / f"{g}_seqnext.tsv"
        if src.exists():
            shutil.copy(src, per_gene_out / src.name)
            shutil.copy(src, bundle_per_gene / src.name)
            copied += 1
    print(f"  per_gene/*.tsv                          {copied:>6,} files copied")

    print(f"\nFiltered outputs at:")
    print(f"  {OUT_DIR}")
    print(f"  {BUNDLE_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

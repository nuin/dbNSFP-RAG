#!/usr/bin/env python3
"""Build the clinical upload file for the 68 net-new cp_new genes.

Filters from cp_new_seqnext_FINAL.tsv:
  1. Only the 68 genes that are NEW in CP_new (APR2026) vs old CP (01JUN2021).
     Established CP genes are excluded -- the lab is not auto-classifying
     them yet to avoid disrupting the existing SeqNext mutation database.
  2. Drop deep-intronic variants. Keep coding + splice-region (-15 to +6
     from exon boundary). Deep intronic positions need separate splice/
     regulatory evaluation; auto-call is unsafe.
  3. vv_status in {ok, clinvar_assertion, local_enumeration} -- drop
     flagged:intergenic rows where VV reports the variant isn't coding on
     the BED's NM_.
  4. Drop canonical splice positions (-1, -2, +1, +2). SpliceAI masked is
     near-zero by design at canonical sites; auto-call would be a filter
     artifact at those positions.

Output:
  cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

NEW_GENES_TXT = Path("cp_new_bundle/outputs/new_genes_only/_new_genes.txt")
FINAL = Path("data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv")
OUT = Path("cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv")

# ROI in splice region (-15 to +6 from nearest exon boundary)
ROI_NEG = -15
ROI_POS = 6

OFFSET_RE = re.compile(r"c\.\d+([+\-])(\d+)")
CANONICAL_SPLICE = re.compile(r"c\.\d+[+\-][12][ACGT]>")


def intronic_offset(hgvsc: str) -> int | None:
    """Return signed intronic offset, or None for coding c. (no offset)."""
    if not hgvsc: return None
    m = OFFSET_RE.search(hgvsc)
    if not m: return None
    return -int(m.group(2)) if m.group(1) == "-" else int(m.group(2))


def main() -> int:
    if not NEW_GENES_TXT.exists():
        sys.exit(f"Missing: {NEW_GENES_TXT}")
    if not FINAL.exists():
        sys.exit(f"Missing: {FINAL}")

    new_genes = set(NEW_GENES_TXT.read_text().split())
    print(f"68-gene whitelist loaded: {len(new_genes)} genes")

    df = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    print(f"FINAL combined rows: {len(df):,}")

    # 1. Restrict to 68 new genes
    df = df[df["gene"].isin(new_genes)].copy()
    print(f"  after new-genes filter: {len(df):,}")

    # 2. Drop vv_status != ok-equivalent
    valid_status = {"ok", "clinvar_assertion", "local_enumeration"}
    n_before = len(df)
    df = df[df["vv_status"].isin(valid_status)]
    print(f"  dropped vv_status not ok: {n_before - len(df):,}  remaining: {len(df):,}")

    # 3. Drop canonical splice (-1/-2/+1/+2 + nucleotide)
    n_before = len(df)
    df = df[~df["hgvs_c"].str.contains(CANONICAL_SPLICE, na=False)]
    print(f"  dropped canonical splice (-1/-2/+1/+2): {n_before - len(df):,}  remaining: {len(df):,}")

    # 4. Apply intronic ROI: keep coding OR offset in [-15, +6]
    n_before = len(df)
    df["_offset"] = df["hgvs_c"].apply(intronic_offset)
    in_roi = df["_offset"].isna() | df["_offset"].between(ROI_NEG, ROI_POS)
    df = df[in_roi].drop(columns=["_offset"])
    print(f"  dropped deep intronic (offset outside [-15,+6]): {n_before - len(df):,}  remaining: {len(df):,}")

    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {len(df):,} upload rows -> {OUT}")
    print(f"\nby source:")
    for s, n in df["source"].value_counts().items():
        print(f"  {s}: {n:,}")
    print(f"\nby gene (top 15):")
    for g, n in df["gene"].value_counts().head(15).items():
        print(f"  {g}: {n}")
    print(f"\nintronic-vs-coding:")
    df["_offset"] = df["hgvs_c"].apply(intronic_offset)
    print(f"  coding (no offset):       {df['_offset'].isna().sum():,}")
    print(f"  splice ROI [-15,+6]:      {df['_offset'].notna().sum():,}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

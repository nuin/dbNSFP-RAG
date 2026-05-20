#!/usr/bin/env python3
"""Build the clinical upload file for the 68 net-new cp_new genes.

Filters from cp_new_seqnext_FINAL.tsv:
  1. Only the 68 genes that are NEW in CP_new (APR2026) vs old CP (01JUN2021).
  2. Drop deep-intronic. Keep coding + splice-region (-15 to +6).
  3. vv_status in {ok, clinvar_assertion, local_enumeration}.
  4. Drop canonical splice (-1, -2, +1, +2) -- SpliceAI masked artifact.
  5. **Drop pure synonymous coding variants** -- they don't change the
     protein. Clinically uninformative; takes up SeqNext database space
     for no decision support. Keep synonymous-flagged variants ONLY if
     they sit in the splice region (-15..+6, where they could affect
     splicing). Detection sources, in order:
       - hgvs_c has intronic offset -> not pure coding, keep
       - Option A catalog row (source=synonymous_catalog) -> drop
       - pipeline classification mentions 'synonymous' but no
         'intronic_ROI' tag -> drop
       - ClinVar clinvar_name contains 'p.XXX=' (synonymous marker) -> drop

Output: cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

NEW_GENES_TXT = Path("cp_new_bundle/outputs/new_genes_only/_new_genes.txt")
FINAL = Path("data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv")
CV = Path("data/exports/cp_new/seqnext/clinvar_benign_seqnext.tsv")
OUT = Path("cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv")

# ROI in splice region (-15 to +6 from nearest exon boundary)
ROI_NEG = -15
ROI_POS = 6

OFFSET_RE = re.compile(r"c\.\d+([+\-])(\d+)")
CANONICAL_SPLICE = re.compile(r"c\.\d+[+\-][12][ACGT]>")
# 'p.Xxx123=' (3-letter) or 'p.X123=' (1-letter) -- ClinVar synonymous marker
SYN_HGVSP = re.compile(r"p\.[A-Za-z]{1,3}\d+=")


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

    # Load ClinVar names for synonymous detection on clinvar-source rows
    cv = pd.read_csv(CV, sep="\t", dtype=str).fillna("") if CV.exists() else pd.DataFrame()
    cv_name = {}
    if not cv.empty:
        cv["key"] = cv["gene"] + "|" + cv["transcript"] + "|" + cv["hgvs_c"]
        cv_name = dict(zip(cv["key"], cv["clinvar_name"]))
    df["_cv_name"] = (df["gene"] + "|" + df["transcript"] + "|" + df["hgvs_c"]).map(cv_name).fillna("")

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
    df = df[in_roi].copy()
    print(f"  dropped deep intronic (offset outside [-15,+6]): {n_before - len(df):,}  remaining: {len(df):,}")

    # 5. Drop pure synonymous coding (no clinical value)
    # Keep ONLY if the variant has an intronic offset (in ROI), regardless of
    # synonymous status -- intronic variants in splice region need review.
    n_before = len(df)
    is_intronic_roi = df["_offset"].notna()  # already filtered to ROI in step 4
    # Synonymous-coding detection across the 3 sources:
    is_catalog = df["source"] == "synonymous_catalog"
    pipeline_syn = (df["source"] == "pipeline") & \
                    df["classification"].str.contains("synonymous", case=False, na=False) & \
                    ~df["classification"].str.contains("intronic_ROI", case=False, na=False)
    clinvar_syn = (df["source"] == "clinvar") & df["_cv_name"].str.contains(SYN_HGVSP, na=False)
    pure_syn = (is_catalog | pipeline_syn | clinvar_syn) & ~is_intronic_roi
    df = df[~pure_syn].drop(columns=["_offset", "_cv_name"])
    print(f"  dropped pure synonymous coding (no clinical value): {n_before - len(df):,}  remaining: {len(df):,}")

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

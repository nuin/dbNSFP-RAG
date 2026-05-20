#!/usr/bin/env python3
"""Merge the ClinVar Benign/LB pull on top of the dbNSFP-pipeline SeqNext output.

Strategy:
  - Start with the existing classifier output (cp_new_seqnext.tsv).
  - For each ClinVar B/LB row that is NOT already in the pipeline output
    (key: gene + transcript + hgvs_c), add it with source='clinvar'.
  - For rows present in both, keep the pipeline call but stamp the classification
    string to note ClinVar agreement.

Output: data/exports/cp_new/seqnext/cp_new_seqnext_with_clinvar.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

DATA = Path("data/exports/cp_new")
PIPE = DATA / "seqnext" / "cp_new_seqnext.tsv"
CV = DATA / "seqnext" / "clinvar_benign_seqnext.tsv"
OUT = DATA / "seqnext" / "cp_new_seqnext_with_clinvar.tsv"


def norm_c(c: str) -> str:
    """Normalize HGVSc for matching (strip whitespace)."""
    return (c or "").strip()


def main() -> int:
    pipe = pd.read_csv(PIPE, sep="\t", dtype=str).fillna("")
    cv = pd.read_csv(CV, sep="\t", dtype=str).fillna("")
    print(f"pipeline rows: {len(pipe):,}")
    print(f"ClinVar B/LB rows: {len(cv):,}")

    pipe["source"] = "pipeline"
    cv["source"] = "clinvar"

    # Match key: (gene, transcript, c.)
    pipe["key"] = pipe["gene"] + "|" + pipe["transcript"] + "|" + pipe["hgvs_c"].apply(norm_c)
    cv["key"] = cv["gene"] + "|" + cv["transcript"] + "|" + cv["hgvs_c"].apply(norm_c)

    pipe_keys = set(pipe["key"])
    cv_only = cv[~cv["key"].isin(pipe_keys)].copy()
    cv_overlap = cv[cv["key"].isin(pipe_keys)].copy()

    print(f"\noverlap (pipeline + ClinVar agree): {len(cv_overlap):,}")
    print(f"ClinVar-only (recovers from pipeline misses): {len(cv_only):,}")

    # Stamp pipeline rows that overlap with ClinVar
    pipe = pipe.set_index("key")
    overlap_keys = set(cv_overlap["key"])
    pipe.loc[pipe.index.isin(overlap_keys), "classification"] = (
        pipe.loc[pipe.index.isin(overlap_keys), "classification"] + " [ClinVar agrees]"
    )
    pipe = pipe.reset_index(drop=True)

    # Append ClinVar-only rows (align columns)
    for c in pipe.columns:
        if c not in cv_only.columns:
            cv_only[c] = ""
    keep = [c for c in pipe.columns]
    cv_only = cv_only[keep]

    merged = pd.concat([pipe, cv_only], ignore_index=True)
    merged.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {len(merged):,} rows -> {OUT}")

    print("\nClinVar-only additions by gene (top 15):")
    for g, n in cv_only["gene"].value_counts().head(15).items():
        print(f"  {g}: {n}")

    # Save also as cp_new_seqnext_minimal/strict-equivalents for the new combined output
    cols_min = ["gene", "transcript", "hgvs_c", "classification",
                "PhastCons100way", "PhyloP100way", "REVEL", "SpliceAI_masked",
                "FAF95_grpmax", "CADD_phred", "AlphaMissense_pred", "ClinVar_sig"]
    cols_min = [c for c in cols_min if c in merged.columns]
    minimal_path = DATA / "cp_new_seqnext_minimal_with_clinvar.tsv"
    merged[cols_min].to_csv(minimal_path, sep="\t", index=False)
    print(f"Minimal:   {minimal_path}")

    # Strict: drop intergenic + canonical-splice from this combined output
    import re
    splice = re.compile(r"c\.\d+[+\-][12][ACGT]>")
    mask_ok = (merged["vv_status"] == "ok") | (merged["vv_status"] == "clinvar_assertion")
    mask_splice = merged["hgvs_c"].str.contains(splice, na=False)
    strict = merged[mask_ok & ~mask_splice].copy()
    strict_path = DATA / "cp_new_seqnext_strict_with_clinvar.tsv"
    strict[cols_min].to_csv(strict_path, sep="\t", index=False)
    print(f"Strict:    {strict_path}  ({len(strict):,} rows)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

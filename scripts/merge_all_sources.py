#!/usr/bin/env python3
"""Merge all classification sources into one definitive SeqNext output:

  1. Pipeline       (dbNSFP nsSNVs through W1/W2/W3 + VV resolver)
  2. ClinVar B/LB   (NCBI ClinVar bulk, >=1 star)
  3. Synonymous     (locally enumerated catalog from hg19 RefSeq + FASTA)

Dedup key: (gene, transcript, hgvs_c). Priority when multiple sources
report the same variant: pipeline > ClinVar > catalog. Each row tagged with
its origin in the 'source' column.

Outputs:
  data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv     -- everything
  data/exports/cp_new/seqnext/cp_new_seqnext_FINAL_strict.tsv -- minus splice/intergenic
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

DATA = Path("data/exports/cp_new/seqnext")
PIPE = DATA / "cp_new_seqnext.tsv"
CV   = DATA / "clinvar_benign_seqnext.tsv"
SYN  = DATA / "cp_new_seqnext_synonymous_catalog.tsv"
OUT  = DATA / "cp_new_seqnext_FINAL.tsv"
OUT_STRICT = DATA / "cp_new_seqnext_FINAL_strict.tsv"

CANONICAL_SPLICE = re.compile(r"c\.\d+[+\-][12][ACGT]>")


def load_or_empty(p: Path, source: str) -> pd.DataFrame:
    if not p.exists():
        print(f"  warn: missing {p}"); return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    if "source" not in df.columns:
        df["source"] = source
    df["_priority"] = {"pipeline": 0, "clinvar": 1, "synonymous_catalog": 2}.get(source, 99)
    return df


def main() -> int:
    pipe = load_or_empty(PIPE, "pipeline")
    cv = load_or_empty(CV, "clinvar")
    syn = load_or_empty(SYN, "synonymous_catalog")

    print(f"pipeline:   {len(pipe):,}")
    print(f"clinvar:    {len(cv):,}")
    print(f"synonymous: {len(syn):,}")

    all_df = pd.concat([pipe, cv, syn], ignore_index=True)
    all_df["key"] = all_df["gene"] + "|" + all_df["transcript"] + "|" + all_df["hgvs_c"]

    # Dedup: keep lowest _priority per key
    all_df = all_df.sort_values("_priority").drop_duplicates(subset=["key"], keep="first")
    all_df = all_df.drop(columns=["_priority","key"])
    print(f"\nFINAL combined (dedup'd): {len(all_df):,}")
    print(f"  by source:")
    for s, n in all_df["source"].value_counts().items():
        print(f"    {s}: {n:,}")

    # Column order: core SeqNext + scores + audit
    cols_core = ["gene","transcript","hgvs_c","classification"]
    cols_scores = ["PhastCons100way","PhyloP100way","REVEL","SpliceAI_masked",
                   "FAF95_grpmax","CADD_phred","AlphaMissense_pred","ClinVar_sig"]
    cols_audit = ["chr_grch37","pos_grch37","ref","alt","vv_status","source"]
    ordered = [c for c in cols_core + cols_scores + cols_audit if c in all_df.columns]
    for c in ordered:
        if c not in all_df.columns: all_df[c] = ""
    all_df[ordered].to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {OUT}")

    # Strict: drop intergenic + canonical splice
    ok = (all_df["vv_status"].isin(["ok", "clinvar_assertion", "local_enumeration"]))
    nosplice = ~all_df["hgvs_c"].str.contains(CANONICAL_SPLICE, na=False)
    strict = all_df[ok & nosplice].copy()
    strict[ordered].to_csv(OUT_STRICT, sep="\t", index=False)
    print(f"Strict (no intergenic, no canonical splice): {len(strict):,} -> {OUT_STRICT}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

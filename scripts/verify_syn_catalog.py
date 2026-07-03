#!/usr/bin/env python3
"""Verify every entry in the synonymous catalog actually is synonymous by
re-translating codon_ref and codon_alt with BioPython and checking they
produce the same amino acid.

Catches any catalog-side bug analogous to the codon_degeneracy issue we
just fixed in classify_cp_new.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from Bio.Data import CodonTable

CAT = Path("data/exports/cp_new/cp_new_synonymous_catalog.tsv")
TABLE = CodonTable.unambiguous_dna_by_id[1]  # standard table


def translate(codon: str) -> str | None:
    codon = (codon or "").upper().strip()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return None
    if codon in TABLE.stop_codons:
        return "*"
    return TABLE.forward_table.get(codon)


def main() -> int:
    if not CAT.exists(): sys.exit(f"missing {CAT}")
    df = pd.read_csv(CAT, sep="\t", dtype=str).fillna("")
    print(f"Catalog rows: {len(df):,}")

    df["_aaref"] = df["codon_ref"].apply(translate)
    df["_aaalt"] = df["codon_alt"].apply(translate)

    bad_translate = df[df["_aaref"].isna() | df["_aaalt"].isna()]
    print(f"  with un-translatable codon: {len(bad_translate):,}")

    classified = df[df["_aaref"].notna() & df["_aaalt"].notna()].copy()
    truly_syn = classified[classified["_aaref"] == classified["_aaalt"]]
    not_syn = classified[classified["_aaref"] != classified["_aaalt"]]
    print(f"  truly synonymous (aaref==aaalt): {len(truly_syn):,}")
    print(f"  NOT synonymous: {len(not_syn):,}")

    if len(not_syn):
        print(f"\n=== non-synonymous mis-tagged in catalog (top 10 by gene) ===")
        print(not_syn.groupby("gene").size().sort_values(ascending=False).head(10).to_string())
        print(f"\nsamples:")
        print(not_syn.head(10)[["gene","transcript","hgvs_c","codon_ref","codon_alt",
                               "_aaref","_aaalt"]].to_string(index=False))

        # Stop-codon hits
        stops = not_syn[(not_syn["_aaalt"] == "*") | (not_syn["_aaref"] == "*")]
        print(f"\n  of those, generating/abolishing stop codon: {len(stops):,}")

        bad_keys = set(zip(not_syn["gene"], not_syn["hgvs_c"]))
        out = Path("data/exports/cp_new/cp_new_syn_catalog_misclassified.tsv")
        not_syn.drop(columns=["_aaref","_aaalt"], errors="ignore").to_csv(out, sep="\t", index=False)
        print(f"\n  list written -> {out}")
    return 0 if len(not_syn) == 0 else 2


if __name__ == "__main__":
    sys.exit(main())

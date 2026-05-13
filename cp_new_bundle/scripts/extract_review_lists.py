#!/usr/bin/env python3
"""Extract the two manual-review lists from the cp_new SeqNext output.

Produces two TSVs in data/exports/cp_new/seqnext/ alongside the main output:

  _review_canonical_splice.tsv
      W2 'Benign synonymous' calls whose c. notation puts them at the
      canonical splice acceptor/donor positions (-1, -2, +1, +2). These
      pass the SpliceAI<=0.1 filter only because SpliceAI '-M 1' is
      designed to zero out scores AT canonical splice sites -- so the
      Benign call is an artifact of using masked mode at these positions.

  _review_intergenic.tsv
      Rows where VariantValidator returned vv_status=flagged:intergenic --
      i.e. the variant is non-coding on the BED's specified RefSeq NM_.
      The classification probably still holds biologically, but the c.
      notation in these rows is from a different transcript than the BED
      target and shouldn't be uploaded to SeqNext as-is.

Both files use the same column layout as cp_new_seqnext.tsv so they can be
diff'd / reviewed / fed back into downstream tools.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

DEFAULT_INPUT = Path("data/exports/cp_new/seqnext/cp_new_seqnext.tsv")
OUT_DIR = Path("data/exports/cp_new/seqnext")

# Canonical splice positions appear in HGVSc as e.g. c.123-1G>A, c.123-2A>C,
# c.456+1G>T, c.789+2T>G. Match c.<digits>{+,-}{1,2}<nucleotide>>
CANONICAL_SPLICE_RE = re.compile(r"c\.\d+[+\-][12][ACGT]>")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = p.parse_args()

    if not args.input.exists():
        sys.exit(f"SeqNext file not found: {args.input}")

    df = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df):,} rows from {args.input}")

    # --- canonical-splice W2 review ---
    is_w2_syn = df["classification"].str.contains("synonymous", case=False, na=False)
    is_canonical_splice = df["hgvs_c"].str.contains(CANONICAL_SPLICE_RE, na=False)
    splice_df = df[is_w2_syn & is_canonical_splice].copy()
    splice_df["review_reason"] = (
        "SpliceAI masked is zero at canonical splice sites by design; "
        "Benign call is filter artifact, not evidence of safety"
    )
    splice_path = args.out_dir / "_review_canonical_splice.tsv"
    splice_df.to_csv(splice_path, sep="\t", index=False)
    print(f"  canonical-splice rows: {len(splice_df):>5} -> {splice_path}")

    # --- intergenic flag review ---
    if "vv_status" in df.columns:
        intergenic_df = df[df["vv_status"] == "flagged:intergenic"].copy()
    else:
        intergenic_df = df[df["hgvs_c"].str.contains("VV_FLAGGED:intergenic", na=False)].copy()
    intergenic_df["review_reason"] = (
        "VariantValidator: variant is non-coding on the BED's RefSeq NM_; "
        "c. notation belongs to a different transcript -- drop or re-resolve "
        "against an alternative isoform"
    )
    interg_path = args.out_dir / "_review_intergenic.tsv"
    intergenic_df.to_csv(interg_path, sep="\t", index=False)
    print(f"  intergenic rows:        {len(intergenic_df):>5} -> {interg_path}")

    # --- breakdown by gene ---
    print("\nBreakdown by gene:")
    print("\n  canonical-splice W2 calls per gene:")
    for gene, n in splice_df["gene"].value_counts().items():
        print(f"    {gene:<15} {n}")
    print("\n  intergenic rows per gene:")
    for gene, n in intergenic_df["gene"].value_counts().items():
        print(f"    {gene:<15} {n}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

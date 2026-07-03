#!/usr/bin/env python3
"""Remove 5'UTR variants (HGVS c.-N, coordinate before the ATG) from the v3.3/v3.4
uploads (full + ROI). Lab request 2026-07: 5'UTR variants should not be in the
list. 3'UTR (c.*N) are left in unless separately requested.
"""
import re
from pathlib import Path
import pandas as pd

TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]
UTR5 = re.compile(r"^c\.-\d")   # c.-45A>T, c.-12-13C>T, etc.


def main():
    for path, pgname, tag in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        mask = df["hgvs_c"].str.match(UTR5, na=False)
        dropped = df[mask]
        df = df[~mask]
        # append to any existing dropped audit
        audit = p.with_name(p.stem + "_5utr_dropped.tsv")
        dropped.to_csv(audit, sep="\t", index=False)
        df.to_csv(p, sep="\t", index=False)
        print(f"  {p.name}: removed {len(dropped):,} 5'UTR -> {len(df):,} rows")
        # refresh per-gene
        pg = p.parent / pgname; pg.mkdir(exist_ok=True)
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{tag}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))
    print("done")


if __name__ == "__main__":
    main()

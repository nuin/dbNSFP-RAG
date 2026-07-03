#!/usr/bin/env python3
"""Ensure EVERY row has PhastCons100way + PhyloP100way, pulled uniformly from
the UCSC hg19 100-way tracks. Fills any remaining blanks (gnomad_common W1,
clinvar, etc.) using each row's hg19 position (chr_grch37/pos_grch37, or parsed
from the g_hg19 column for gnomad_common rows that lack native hg19 coords).

No dropping here -- for W1/W3/clinvar rows conservation is informational, not a
gate. (The synonymous W2 gate was already applied in fill_conservation.py.)
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd
import pyBigWig

PHASTCONS = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw"
PHYLOP    = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw"
TABLES = [
    "cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv",
]
G19 = re.compile(r"chr([\dXYM]+):g\.(\d+)")


def hg19_pos(row):
    """Return (chrom, pos) in hg19, from native coords or the g_hg19 column."""
    if row["chr_grch37"] and row["pos_grch37"]:
        return row["chr_grch37"], int(float(row["pos_grch37"]))
    m = G19.search(row.get("g_hg19", "") or "")
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def main():
    # Gather positions needing fill, grouped by GENE (tight ranges; two genes
    # far apart on the same chromosome must NOT be read as one giant span).
    needed = {}  # gene -> (chrom, set(pos))
    for path in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        miss = (df["PhastCons100way"]=="") | (df["PhyloP100way"]=="")
        for _, r in df[miss].iterrows():
            c, pos = hg19_pos(r)
            if c and pos:
                g = r["gene"]
                if g not in needed: needed[g] = (c, set())
                needed[g][1].add(pos)
    total = sum(len(s) for _, s in needed.values())
    print(f"positions needing conservation: {total:,} across {len(needed)} genes")

    pc = pyBigWig.open(PHASTCONS); pp = pyBigWig.open(PHYLOP)
    cache = {}
    for gene, (chrom, positions) in needed.items():
        lo, hi = min(positions), max(positions)
        c = f"chr{chrom}"
        try:
            pcv = pc.values(c, lo-1, hi); ppv = pp.values(c, lo-1, hi)
        except Exception as e:
            print(f"  WARN {gene} {c}:{lo}-{hi}: {e}"); continue
        for pos in positions:
            i = pos - lo
            a, b = pcv[i], ppv[i]
            cache[(chrom,pos)] = (None if a is None or a!=a else float(a),
                                  None if b is None or b!=b else float(b))
        print(f"  {gene} {c}: {len(positions):,} positions ({hi-lo+1:,} bp span)", flush=True)
    pc.close(); pp.close()

    for path in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        n = 0
        for i in df.index:
            r = df.loc[i]
            if r["PhastCons100way"] and r["PhyloP100way"]:
                continue
            c, pos = hg19_pos(r)
            v = cache.get((c,pos)) if c else None
            if not v: continue
            if not r["PhastCons100way"] and v[0] is not None:
                df.at[i,"PhastCons100way"] = f"{v[0]:.4f}"
            if not r["PhyloP100way"] and v[1] is not None:
                df.at[i,"PhyloP100way"] = f"{v[1]:.4f}"
            n += 1
        df.to_csv(p, sep="\t", index=False)
        still = ((df["PhastCons100way"]=="")).sum()
        print(f"  {p.name}: filled {n:,} | still-missing PhastCons {still:,}")

    # refresh per-gene
    for folder, stem, tag in [
        ("cp_new_bundle/outputs/new_genes_68_v3_3","UPLOAD_v3_3","v3_3"),
        ("cp_new_bundle/outputs/new_genes_68_v3_4","UPLOAD_v3_4","v3_4"),
    ]:
        for src, pgdir, suf in [(f"{folder}/{stem}.tsv", f"{folder}/per_gene", tag),
                                (f"{folder}/{stem}_roi.tsv", f"{folder}/per_gene_roi", f"{tag}_roi")]:
            if not Path(src).exists(): continue
            d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
            pg = Path(pgdir); pg.mkdir(exist_ok=True)
            for old in pg.glob("*.tsv"): old.unlink()
            for g in sorted(d["gene"].unique()):
                (pg/f"{g}__{suf}.tsv").write_text(d[d["gene"]==g].to_csv(sep="\t", index=False))
    print("per-gene refreshed")


if __name__ == "__main__":
    main()

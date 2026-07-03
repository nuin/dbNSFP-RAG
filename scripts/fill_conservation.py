#!/usr/bin/env python3
"""Fill PhastCons100way / PhyloP100way for synonymous rows marked
"PhastCons unmeasured" using the UCSC hg19 100-way conservation bigWigs
(read remotely, per-gene range). Lab needs nucleotide conservation to apply
BP4/BP7. Rows where the track has no value are DROPPED (per lab: "if they are
not available then I would rather not have them").

Also re-applies the W2 conservation gate now that real values exist:
  coding synonymous: keep if PhastCons < 1.0
  intronic (ROI):    keep if PhastCons < 1.0 AND PhyloP < 0.1

Updates all 4 tables (v3.3/v3.4 x full/roi) in place + per-gene refresh.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pyBigWig

PHASTCONS = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw"
PHYLOP    = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw"

TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]

import re
OFF = re.compile(r"c\.[\*\-]?\d+([+\-])(\d+)")
def is_intronic(h):
    return bool(OFF.search(h or ""))


def build_cons_cache(genes_positions):
    """genes_positions: dict gene -> (chrom, minpos, maxpos). Returns
    {(chrom,pos): (phastcons, phylop)} across all gene ranges."""
    pc = pyBigWig.open(PHASTCONS)
    pp = pyBigWig.open(PHYLOP)
    cache = {}
    for gene, (chrom, lo, hi) in genes_positions.items():
        c = f"chr{chrom}"
        try:
            pcv = pc.values(c, lo-1, hi)     # 0-based half-open; covers pos lo..hi (1-based)
            ppv = pp.values(c, lo-1, hi)
        except Exception as e:
            print(f"  WARN {gene} {c}:{lo}-{hi}: {e}", file=sys.stderr)
            continue
        for i, pos in enumerate(range(lo, hi+1)):
            a = pcv[i]; b = ppv[i]
            cache[(chrom, pos)] = (
                None if (a is None or a != a) else float(a),
                None if (b is None or b != b) else float(b),
            )
        print(f"  cached {gene} {c}:{lo}-{hi} ({hi-lo+1:,} bases)", flush=True)
    pc.close(); pp.close()
    return cache


def main():
    # Per-gene ranges from the UNION of all tables' "unmeasured" rows (tight
    # ranges per gene, and don't depend on any single table having the rows).
    gene_ranges = {}  # gene -> (chrom, minpos, maxpos)
    for path, _, _ in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        need = df[df["classification"].str.contains("PhastCons unmeasured", na=False)]
        for gene, g in need.groupby("gene"):
            pos = pd.to_numeric(g["pos_grch37"], errors="coerce").dropna().astype(int)
            if not len(pos): continue
            chrom = g["chr_grch37"].iloc[0]
            lo, hi = int(pos.min()), int(pos.max())
            if gene not in gene_ranges:
                gene_ranges[gene] = (chrom, lo, hi)
            else:
                _, plo, phi = gene_ranges[gene]
                gene_ranges[gene] = (chrom, min(plo, lo), max(phi, hi))
    print(f"Reading UCSC conservation for {len(gene_ranges)} gene ranges (remote bigWig)...")
    cache = build_cons_cache(gene_ranges)
    print(f"Cached {len(cache):,} positions\n")

    for path, pgname, tag in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        unm = df["classification"].str.contains("PhastCons unmeasured", na=False)
        # KEEP all rows (lab will filter). Fill conservation; relabel by whether
        # the position passes the W2/BP7 gate so column D is honest and the lab
        # can filter on it. Only truly-unavailable (no track value) are dropped.
        drop_idx, n_benign, n_conserved = [], 0, 0
        for i in df[unm].index:
            r = df.loc[i]
            key = (r["chr_grch37"], int(float(r["pos_grch37"]))) if r["pos_grch37"] else None
            vals = cache.get(key) if key else None
            if not vals or vals[0] is None:
                drop_idx.append(i)   # no conservation at all -> drop
                continue
            pc_val, pp_val = vals
            df.at[i, "PhastCons100way"] = f"{pc_val:.4f}"
            df.at[i, "PhyloP100way"] = "" if pp_val is None else f"{pp_val:.4f}"
            conserved = (pc_val >= 1.0) or (is_intronic(r["hgvs_c"]) and (pp_val is None or pp_val >= 0.1))
            if conserved:
                df.at[i, "classification"] = ("Synonymous, SpliceAI<=0.1 (assumed 0), "
                                              "PhastCons>=1.0 (conserved - review, not auto-benign)")
                n_conserved += 1
            else:
                df.at[i, "classification"] = ("Benign synonymous, SpliceAI<=0.1 (assumed 0), "
                                              "PhastCons<1.0")
                n_benign += 1
        df = df.drop(index=drop_idx)
        df.to_csv(p, sep="\t", index=False)
        print(f"  {p.name}: benign {n_benign:,} | conserved-kept {n_conserved:,} | "
              f"dropped(no-track) {len(drop_idx):,} -> {len(df):,} rows")

        # per-gene refresh
        folder = p.parent
        pg = folder / pgname; pg.mkdir(exist_ok=True)
        suffix = tag
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{suffix}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))
    print("done")


if __name__ == "__main__":
    main()

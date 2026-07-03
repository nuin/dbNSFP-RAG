#!/usr/bin/env python3
"""Build a VCF of variants from the v2 upload that still lack a measured
SpliceAI score, ready to be fed to `spliceai -M 1 -A grch38`.

Translates each row's hg19 coords to hg38 via the dbNSFP per-gene position
maps, sorts + dedupes by (chr, pos, ref, alt), writes a minimal sites-only VCF.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

UP = Path("cp_new_bundle/outputs/new_genes_68_v2/UPLOAD_v2_scored.tsv")
PER_GENE_DBNSFP = Path("data/exports/cp_new")
OUT_VCF = Path("data/exports/cp_new/cp_new_v2_needed.vcf")
CHR_ORDER = [*[str(i) for i in range(1, 23)], "X", "Y", "M"]


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError): pass
    return str(v)


def build_hg19_to_hg38(genes: list[str]) -> dict[tuple[str, int], int]:
    out = {}
    for g in genes:
        p = PER_GENE_DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        if "hg19_chr" not in df.columns: continue
        for _, r in df.iterrows():
            try:
                k = (_s(r["hg19_chr"]), int(float(r["hg19_pos(1-based)"])))
                v = int(float(r["pos(1-based)"]))
                out.setdefault(k, v)
            except (TypeError, ValueError, KeyError): continue
    return out


def main() -> int:
    if not UP.exists(): sys.exit(f"missing {UP}")
    df = pd.read_csv(UP, sep="\t", dtype=str).fillna("")

    need = df[
        df["SpliceAI_masked"].isin(["", "0 (assumed)"]) &
        (df["ref"] != "") & (df["alt"] != "") &
        (df["pos_grch37"] != "") & (df["chr_grch37"] != "")
    ].copy()
    print(f"Upload rows needing real SpliceAI: {len(need):,}")

    genes = sorted(need["gene"].unique())
    print(f"Building hg19->hg38 map from {len(genes)} dbNSFP TSVs ...")
    hg = build_hg19_to_hg38(genes)
    print(f"  {len(hg):,} positions mapped\n")

    sites = set()
    n_no_map = 0
    for _, r in need.iterrows():
        try:
            hg19_pos = int(float(r["pos_grch37"]))
        except (TypeError, ValueError): continue
        hg38_pos = hg.get((r["chr_grch37"], hg19_pos))
        if hg38_pos is None:
            n_no_map += 1
            continue
        # SpliceAI handles only simple SNVs cleanly; pass everything anyway
        sites.add((r["chr_grch37"], hg38_pos, r["ref"], r["alt"]))

    print(f"Unique (chr, hg38_pos, ref, alt) sites: {len(sites):,}")
    print(f"  rows without hg19->hg38 map (skipped): {n_no_map:,}\n")

    # Sort by chromosome then position
    def chr_key(c):
        return CHR_ORDER.index(c) if c in CHR_ORDER else 99
    sites_sorted = sorted(sites, key=lambda x: (chr_key(x[0]), x[1], x[2], x[3]))

    OUT_VCF.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_VCF, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=build_spliceai_input_vcf.py (cp_new v2 backfill)\n")
        f.write("##INFO=<ID=.,Number=0,Type=Flag,Description=\"placeholder\">\n")
        # write contig lines for the chromosomes present
        for c in sorted({s[0] for s in sites_sorted}, key=chr_key):
            f.write(f"##contig=<ID=chr{c}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, pos, ref, alt in sites_sorted:
            f.write(f"chr{c}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"Wrote {len(sites_sorted):,} sites -> {OUT_VCF}")
    print(f"\nNext step:")
    print(f"  cd {OUT_VCF.parent} && \\")
    print(f"  source ~/data/spliceai/venv-spliceai/bin/activate && \\")
    print(f"  spliceai -I {OUT_VCF.name} -O cp_new_v2_needed.spliceai.vcf \\")
    print(f"           -R ~/data/hg38/hg38.fa -A grch38 -M 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())

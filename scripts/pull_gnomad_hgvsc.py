#!/usr/bin/env python3
"""Pull gnomAD v4.1 VEP HGVSc for each gene region and build a
(gene, NM_, c.) membership set — the liftover-free, transcript-correct
gnomAD anchor the lab asked for.

gnomAD v4 VEP carries MANE_SELECT (the RefSeq NM_) next to each transcript's
HGVSc. For the MANE transcript the c. numbering == the NM_ numbering, so
(SYMBOL, MANE_SELECT, c.) is a clean key that needs no genome-build liftover
and can't be confused across transcripts (each canonical c. is unique to its
genomic variant).

Pulls BOTH exomes and genomes VCFs (a variant observed in either counts).

Output: data/exports/cp_new/gnomad_hgvsc/_observed_cdot.tsv
  columns: gene, nm, cdot  (one row per observed (gene,NM_,c.))
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

EXOMES = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
          "release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz")
GENOMES = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
           "release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz")
RANGES = Path("cp_new_bundle/outputs/_shared/gene_hg38_ranges.json")
OUT_DIR = Path("data/exports/cp_new/gnomad_hgvsc")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# VEP field indices (0-based) per the Format header
VEP_SYMBOL = 3
VEP_HGVSC = 10
VEP_MANE_SELECT = 25


def tabix(url, chrom, start, end):
    try:
        r = subprocess.run(["tabix", url, f"chr{chrom}:{start}-{end}"],
                           capture_output=True, text=True, timeout=1800, check=False)
    except subprocess.TimeoutExpired:
        return None
    if r.returncode != 0:
        return None
    return r.stdout.splitlines()


def parse_vep_cdots(lines, want_symbol):
    """Yield (gene, nm, cdot) for every transcript entry that has a MANE_SELECT
    NM_ and an HGVSc, restricted to the gene of interest."""
    out = set()
    for line in lines or []:
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        info = parts[7]
        vep = None
        for kv in info.split(";"):
            if kv.startswith("vep="):
                vep = kv[4:]
                break
        if not vep:
            continue
        for entry in vep.split(","):
            f = entry.split("|")
            if len(f) <= VEP_MANE_SELECT:
                continue
            symbol = f[VEP_SYMBOL]
            if symbol != want_symbol:
                continue
            mane = f[VEP_MANE_SELECT].strip()
            hgvsc = f[VEP_HGVSC].strip()
            if not mane or not hgvsc:
                continue
            # hgvsc like 'ENST00000265081.7:c.1992G>A' -> 'c.1992G>A'
            if ":c." in hgvsc:
                cdot = "c." + hgvsc.split(":c.", 1)[1]
            else:
                continue
            out.add((symbol, mane, cdot))
    return out


def pull_gene(gene, info, buffer=50):
    chrom = info["chr"]
    start = max(1, info["start"] - buffer)
    end = info["end"] + buffer
    found = set()
    for url_tmpl in (EXOMES, GENOMES):
        lines = tabix(url_tmpl.format(chrom=chrom), chrom, start, end)
        if lines is None:
            print(f"  WARN {gene} pull failed ({'exomes' if url_tmpl is EXOMES else 'genomes'})",
                  file=sys.stderr)
            continue
        found |= parse_vep_cdots(lines, gene)
    return gene, found


def main():
    ranges = json.loads(RANGES.read_text())
    genes = sys.argv[1:] or sorted(ranges)
    ranges = {g: ranges[g] for g in genes if g in ranges}
    print(f"Pulling gnomAD v4.1 VEP HGVSc (exomes+genomes) for {len(ranges)} genes...")

    all_rows = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(pull_gene, g, info): g for g, info in ranges.items()}
        for i, fut in enumerate(as_completed(futs), 1):
            gene, found = fut.result()
            for (sym, nm, cdot) in found:
                all_rows.append({"gene": sym, "nm": nm, "cdot": cdot})
            print(f"  [{i}/{len(ranges)}] {gene:<10} {len(found):>6} (gene,NM_,c.) "
                  f"({time.time()-t0:.0f}s)", flush=True)

    df = pd.DataFrame(all_rows).drop_duplicates()
    out = OUT_DIR / "_observed_cdot.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nTotal observed (gene,NM_,c.): {len(df):,}")
    print(f"-> {out}")
    print(f"\nby gene (top 10):")
    print(df["gene"].value_counts().head(10).to_string())


if __name__ == "__main__":
    main()

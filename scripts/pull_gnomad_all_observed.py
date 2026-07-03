#!/usr/bin/env python3
"""Pull ALL gnomAD v4.1 observed variants (any frequency) for each gene region,
to use as a gnomAD anchor: only variants observed in gnomAD stay in the upload.

Outputs data/exports/cp_new/gnomad_observed/{GENE}.observed.tsv with
(chr_hg38, pos_hg38, ref, alt, faf95, af) and a combined set file.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

GNOMAD_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
              "release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz")
RANGES = Path("cp_new_bundle/outputs/_shared/gene_hg38_ranges.json")
OUT_DIR = Path("data/exports/cp_new/gnomad_observed")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_info(info: str) -> dict:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def tabix(url, chrom, start, end):
    try:
        r = subprocess.run(["tabix", url, f"chr{chrom}:{start}-{end}"],
                           capture_output=True, text=True, timeout=900, check=False)
    except subprocess.TimeoutExpired:
        return []
    if r.returncode != 0: return []
    return [l for l in r.stdout.splitlines() if l and not l.startswith("#")]


def pull_gene(gene, info, buffer=50):
    chrom = info["chr"]
    start = max(1, info["start"] - buffer)
    end = info["end"] + buffer
    rows = []
    for line in tabix(GNOMAD_URL.format(chrom=chrom), chrom, start, end):
        parts = line.split("\t")
        if len(parts) < 8: continue
        c, pos, _id, ref, alt, _q, vfilter, info_s = parts[:8]
        inf = parse_info(info_s)
        faf = inf.get("fafmax_faf95_max_joint") or inf.get("fafmax_faf95_max") or ""
        af = inf.get("AF_joint") or inf.get("AF") or ""
        rows.append({
            "gene": gene, "chr_hg38": c.removeprefix("chr"), "pos_hg38": pos,
            "ref": ref, "alt": alt, "faf95": faf, "af": af, "filter": vfilter,
        })
    return gene, pd.DataFrame(rows)


def main():
    ranges = json.loads(RANGES.read_text())
    print(f"Pulling ALL gnomAD v4.1 observed variants for {len(ranges)} genes...")
    all_dfs = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(pull_gene, g, info): g for g, info in ranges.items()}
        for i, fut in enumerate(as_completed(futs), 1):
            gene, df = fut.result()
            if not df.empty:
                df.to_csv(OUT_DIR / f"{gene}.observed.tsv", sep="\t", index=False)
                all_dfs.append(df)
            if i % 10 == 0 or i == len(futs):
                tot = sum(len(d) for d in all_dfs)
                print(f"  [{i}/{len(futs)}] {gene:<10} cached {tot:,} observed ({time.time()-t0:.0f}s)")
    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(OUT_DIR / "_all_observed.tsv", sep="\t", index=False)
    print(f"\nTotal observed gnomAD variants: {len(combined):,}")
    print(f"-> {OUT_DIR}/_all_observed.tsv")


if __name__ == "__main__":
    main()

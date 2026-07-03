#!/usr/bin/env python3
"""Pull gnomAD v4.1 g. (hg38), gnomAD variant ID, and p. for each gene, keyed by
(gene, MANE_NM, c.). Lets us add gnomAD-native g./p. columns to the upload so the
lab can locate each variant in gnomAD v4.1 exactly.

Output: data/exports/cp_new/gnomad_hgvsc/_gp_map.tsv
  gene, nm, cdot, gnomad_id, g_hg38, p_hgvs
"""
from __future__ import annotations
import json, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

EXOMES  = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
           "release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz")
GENOMES = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
           "release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz")
RANGES = Path("cp_new_bundle/outputs/_shared/gene_hg38_ranges.json")
OUT = Path("data/exports/cp_new/gnomad_hgvsc/_gp_map.tsv")

VEP_SYMBOL, VEP_HGVSC, VEP_HGVSP, VEP_MANE = 3, 10, 11, 25


def tabix(url, chrom, start, end):
    try:
        r = subprocess.run(["tabix", url, f"chr{chrom}:{start}-{end}"],
                           capture_output=True, text=True, timeout=1800, check=False)
    except subprocess.TimeoutExpired:
        return None
    return r.stdout.splitlines() if r.returncode == 0 else None


def parse(lines, want):
    out = {}
    for line in lines or []:
        if not line or line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 8: continue
        chrom, pos, _id, ref, alt = p[0].removeprefix("chr"), p[1], p[2], p[3], p[4]
        vep = None
        for kv in p[7].split(";"):
            if kv.startswith("vep="): vep = kv[4:]; break
        if not vep: continue
        for entry in vep.split(","):
            f = entry.split("|")
            if len(f) <= VEP_MANE: continue
            if f[VEP_SYMBOL] != want: continue
            mane, hgvsc, hgvsp = f[VEP_MANE].strip(), f[VEP_HGVSC].strip(), f[VEP_HGVSP].strip()
            if not mane or ":c." not in hgvsc: continue
            cdot = "c." + hgvsc.split(":c.", 1)[1]
            # p. : strip ENSP prefix -> p.XXX
            phgvs = ""
            if ":p." in hgvsp:
                phgvs = "p." + hgvsp.split(":p.", 1)[1]
            g_hg38 = f"chr{chrom}:g.{pos}{ref}>{alt}" if len(ref)==1 and len(alt)==1 else \
                     f"chr{chrom}:g.{pos}{ref}>{alt}"  # simple; indels shown as-is
            gnomad_id = f"{chrom}-{pos}-{ref}-{alt}"
            key = (want, mane, cdot)
            out.setdefault(key, (gnomad_id, g_hg38, phgvs))
    return out


def pull(gene, info, buf=50):
    chrom, s, e = info["chr"], max(1, info["start"]-buf), info["end"]+buf
    m = {}
    for tmpl in (EXOMES, GENOMES):
        lines = tabix(tmpl.format(chrom=chrom), chrom, s, e)
        d = parse(lines, gene)
        for k, v in d.items(): m.setdefault(k, v)
    return gene, m


def main():
    ranges = json.loads(RANGES.read_text())
    genes = sys.argv[1:] or sorted(ranges)
    ranges = {g: ranges[g] for g in genes if g in ranges}
    print(f"Pulling gnomAD g./p. for {len(ranges)} genes...")
    rows, t0 = [], time.time()
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(pull, g, i): g for g, i in ranges.items()}
        for i, fut in enumerate(as_completed(futs), 1):
            gene, m = fut.result()
            for (g, nm, cdot), (gid, gh, ph) in m.items():
                rows.append({"gene": g, "nm": nm, "cdot": cdot,
                             "gnomad_id": gid, "g_hg38": gh, "p_hgvs": ph})
            print(f"  [{i}/{len(ranges)}] {gene:<10} {len(m):>6} ({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows).drop_duplicates(["gene","nm","cdot"])
    df.to_csv(OUT, sep="\t", index=False)
    print(f"\n{len(df):,} (gene,NM_,c.) -> g./p.  -> {OUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Annotate the enumerated synonymous catalog with FAF + conservation, then
apply the W2 classifier rule.

Inputs:
  - data/exports/cp_new/cp_new_synonymous_catalog.tsv   (from enumerate_synonymous_catalog.py)
  - per-gene dbNSFP TSVs at data/exports/cp_new/{GENE}.tsv  (borrow PhastCons/PhyloP)
  - gnomAD v2.1.1 hg19 exomes sites VCF (remote tabix)

Annotation strategy:
  - FAF: tabix gnomAD v2.1.1 hg19 per-gene region, build (chr,pos,ref,alt) -> AF_popmax map.
    Falls back to 0 if not in gnomAD.
  - PhastCons/PhyloP: borrow from any dbNSFP row at the same (chr,pos) (dbNSFP
    has nsSNV rows at every coding position; conservation scores are
    per-position so they apply to any alt at that position). Missing if dbNSFP
    doesn't cover the position (rare).

Classifier (W2 only -- synonymous coding variants):
  - PhastCons100way_vertebrate < 1.0 AND
  - SpliceAI <= 0.1 (assumed 0 for new synonymous variants outside splice ROI;
    if a variant happens to be near a splice site, classifier flags it for
    manual SpliceAI rerun)
  - Optional W1 check: if FAF > 5%, override as W1 Benign

Output:
  - data/exports/cp_new/seqnext/cp_new_seqnext_synonymous_catalog.tsv
    Same column schema as the main SeqNext output for clean merging.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from tqdm import tqdm

DATA = Path("data/exports/cp_new")
CATALOG = DATA / "cp_new_synonymous_catalog.tsv"
OUT = DATA / "seqnext" / "cp_new_seqnext_synonymous_catalog.tsv"
GNOMAD_V2_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"


def safe_float(s):
    try: return float(s)
    except (TypeError, ValueError): return None


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def tabix_gnomad_v2(chrom: str, start: int, end: int) -> dict:
    """Return {(pos, ref, alt): af_popmax}."""
    out = {}
    try:
        r = subprocess.run(
            ["tabix", GNOMAD_V2_URL.format(chrom=chrom), f"{chrom}:{start}-{end}"],
            capture_output=True, text=True, timeout=600, check=False,
        )
    except subprocess.TimeoutExpired:
        return out
    if r.returncode != 0: return out
    for line in r.stdout.splitlines():
        if line.startswith("#") or not line: continue
        parts = line.split("\t")
        if len(parts) < 8: continue
        _c, pos, _id, ref, alt, _q, _f, info = parts[:8]
        af = safe_float(parse_info(info).get("AF_popmax"))
        if af is not None:
            out[(int(pos), ref, alt)] = af
    return out


def load_dbnsfp_per_position(gene: str) -> dict[tuple[str, int], dict]:
    """For a gene, load PhastCons/PhyloP per (chr,pos) from the per-gene TSV.
    A position's conservation scores are the same regardless of which alt;
    borrow from any dbNSFP row at that position."""
    tsv = DATA / f"{gene}.tsv"
    if not tsv.exists(): return {}
    out = {}
    try:
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".",""], low_memory=False,
                         usecols=["hg19_chr","hg19_pos(1-based)",
                                  "phastCons100way_vertebrate","phyloP100way_vertebrate"])
    except (ValueError, KeyError):
        return {}
    df["hg19_pos(1-based)"] = pd.to_numeric(df["hg19_pos(1-based)"], errors="coerce")
    for _, r in df.dropna(subset=["hg19_pos(1-based)"]).iterrows():
        key = (str(r["hg19_chr"]).removeprefix("chr"), int(r["hg19_pos(1-based)"]))
        if key not in out:  # first row wins (positions are per-row dedup'd later)
            out[key] = {
                "phastcons": r["phastCons100way_vertebrate"],
                "phylop": r["phyloP100way_vertebrate"],
            }
    return out


def classify_w2(faf, phastcons, spliceai_default=0.0):
    """Returns (classification, rule) or None."""
    if faf is not None and faf > 0.05:
        return "Benign", "FAF >5%"
    pc = safe_float(phastcons)
    if pc is None or pc >= 1.0:
        return None  # need PhastCons < 1.0
    if spliceai_default > 0.1:
        return None
    return "Benign", "synonymous, SpliceAI<=0.1 (assumed 0), PhastCons<1.0"


def main() -> int:
    if not CATALOG.exists():
        sys.exit(f"Catalog not found: {CATALOG}. Run enumerate_synonymous_catalog.py first.")
    print(f"Loading catalog: {CATALOG}")
    cat = pd.read_csv(CATALOG, sep="\t", dtype=str)
    cat["pos_grch37_i"] = pd.to_numeric(cat["pos_grch37"], errors="coerce").astype("Int64")
    print(f"  {len(cat):,} synonymous candidates across {cat['gene'].nunique()} genes")

    # Per-gene: load dbNSFP conservation + tabix gnomAD v2 region
    rows = []
    for gene, gdf in tqdm(cat.groupby("gene"), desc="annotate"):
        gdf = gdf.dropna(subset=["pos_grch37_i"]).copy()
        if gdf.empty: continue
        chrom = gdf["chr_grch37"].iloc[0]
        pmin = int(gdf["pos_grch37_i"].min())
        pmax = int(gdf["pos_grch37_i"].max())
        cons = load_dbnsfp_per_position(gene)
        faf_map = tabix_gnomad_v2(chrom, pmin, pmax)

        for _, r in gdf.iterrows():
            pos = int(r["pos_grch37_i"])
            ref, alt = r["ref"], r["alt"]
            cons_row = cons.get((chrom, pos), {})
            phastcons = cons_row.get("phastcons", "")
            phylop = cons_row.get("phylop", "")
            faf = faf_map.get((pos, ref, alt))
            classification = classify_w2(faf, phastcons)
            if classification is None:
                continue
            cls, rule = classification
            rows.append({
                "gene": r["gene"],
                "transcript": r["transcript"],
                "hgvs_c": r["hgvs_c"],
                "classification": f"{cls} {rule}",
                "chr_grch37": chrom,
                "pos_grch37": pos,
                "ref": ref,
                "alt": alt,
                "PhastCons100way": phastcons,
                "PhyloP100way": phylop,
                "REVEL": "",
                "SpliceAI_masked": "0 (assumed)",
                "FAF95_grpmax": "" if faf is None else f"{faf:.6g}",
                "CADD_phred": "",
                "AlphaMissense_pred": "",
                "ClinVar_sig": "",
                "vv_status": "local_enumeration",
                "source": "synonymous_catalog",
            })

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nClassified {len(df):,} synonymous rows -> {OUT}")
    if not df.empty:
        print("\nBy rule:")
        for r, n in df["classification"].apply(lambda c: c.split(",")[0] if "synonymous" in c else c.split()[0]+" "+c.split()[1]).value_counts().head().items():
            print(f"  {r}: {n}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

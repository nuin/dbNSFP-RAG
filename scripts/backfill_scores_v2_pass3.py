#!/usr/bin/env python3
"""Third backfill pass: enrich the remaining ClinVar rows with what we can get.

For the ~15,566 ClinVar rows still without scores after passes 1+2 (dbNSFP +
syn catalog), do:
  1. Backfill ref/alt from ClinVar bulk (variant_summary.txt.gz)
  2. Coord-based dbNSFP retry (handles non-canonical-transcript missense)
  3. Per-gene no-filter gnomAD tabix -> populate FAF95 by exact (chr,pos,ref,alt)

Updates cp_new_seqnext_UPLOAD_v2_scored.tsv in place + refreshes per-gene files.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

UP = Path("cp_new_bundle/outputs/v2_gnomad_inclusive/cp_new_seqnext_UPLOAD_v2_scored.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/v2_gnomad_inclusive/per_gene")
PER_GENE_DBNSFP = Path("data/exports/cp_new")
CLINVAR_VS = Path.home() / "data/clinvar/variant_summary.txt.gz"
RANGES_JSON = Path("cp_new_bundle/outputs/new_genes_only/_gene_hg38_ranges.json")
GNOMAD_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
              "release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError): pass
    return str(v)


def parse_info(info: str) -> dict:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def tabix_region(url: str, region: str) -> list[str]:
    try:
        r = subprocess.run(["tabix", url, region], capture_output=True,
                           text=True, timeout=900, check=False)
    except subprocess.TimeoutExpired:
        return []
    if r.returncode != 0: return []
    return [l for l in r.stdout.splitlines() if l and not l.startswith("#")]


def load_clinvar_coords(genes: set[str]) -> dict[tuple, tuple]:
    """Map (gene, hgvs_c) -> (chr, pos, ref, alt) from ClinVar GRCh37 rows."""
    out = {}
    with gzip.open(CLINVAR_VS, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        i_asm = idx["Assembly"]
        i_chr = idx["Chromosome"]
        i_pos = idx["PositionVCF"]
        i_ref = idx["ReferenceAlleleVCF"]
        i_alt = idx["AlternateAlleleVCF"]
        i_name = idx["Name"]
        i_gene = idx["GeneSymbol"]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_alt: continue
            if parts[i_asm] != "GRCh37": continue
            gene = parts[i_gene]
            if gene not in genes: continue
            name = parts[i_name]
            # Extract NM_xxx:c.xxx from Name -- format: "NM_002080.5(GOT2):c.816C>T (p.Cys272=)"
            if "):c." not in name: continue
            try:
                tx_paren, rest = name.split("):c.", 1)
                hgvs_c = "c." + rest.split(" ", 1)[0]
            except Exception:
                continue
            ref = parts[i_ref] or ""
            alt = parts[i_alt] or ""
            chrom = parts[i_chr]
            pos = parts[i_pos]
            if not chrom or not pos or not ref or not alt: continue
            out.setdefault((gene, hgvs_c), (chrom, pos, ref, alt))
    return out


def load_dbnsfp_by_coord(gene: str) -> dict[tuple, dict]:
    p = PER_GENE_DBNSFP / f"{gene}.tsv"
    if not p.exists(): return {}
    df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
    out = {}
    SCORE_MAP = {
        "PhastCons100way":   "phastCons100way_vertebrate",
        "PhyloP100way":      "phyloP100way_vertebrate",
        "REVEL":             "REVEL_score",
        "CADD_phred":        "CADD_phred",
        "AlphaMissense_pred":"AlphaMissense_pred",
    }
    for _, r in df.iterrows():
        try:
            k = (_s(r["hg19_chr"]), int(float(r["hg19_pos(1-based)"])),
                 _s(r["ref"]), _s(r["alt"]))
        except (ValueError, TypeError, KeyError): continue
        if k in out: continue
        scores = {}
        for uc, dc in SCORE_MAP.items():
            v = _s(r.get(dc))
            if ";" in v:
                for p_ in v.split(";"):
                    if p_.strip() and p_.strip() != ".":
                        v = p_.strip(); break
                else: v = ""
            scores[uc] = v
        out[k] = scores
    return out


def fetch_gnomad_region_faf(chrom: str, start: int, end: int) -> dict[tuple, str]:
    """Return {(chr, hg38_pos, ref, alt): FAF95_str} for an entire region (no FAF filter)."""
    url = GNOMAD_URL.format(chrom=chrom)
    out = {}
    for line in tabix_region(url, f"chr{chrom}:{start}-{end}"):
        parts = line.split("\t")
        if len(parts) < 8: continue
        c, pos, _id, ref, alt, _q, _f, info = parts[:8]
        ifields = parse_info(info)
        faf = ifields.get("fafmax_faf95_max_joint") or ifields.get("fafmax_faf95_max")
        if not faf or faf == ".": continue
        out[(c.removeprefix("chr"), int(pos), ref, alt)] = faf
    return out


def build_hg38_to_hg19_map(genes: set[str]) -> dict[tuple, str]:
    m = {}
    for g in sorted(genes):
        p = PER_GENE_DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        if "#chr" not in df.columns: continue
        for _, r in df.iterrows():
            try:
                k = (_s(r["#chr"]), int(float(r["pos(1-based)"])))
                v = _s(r.get("hg19_pos(1-based)"))
                if v: m.setdefault(k, v)
            except (TypeError, ValueError, KeyError): continue
    return m


def main() -> int:
    if not UP.exists(): sys.exit(f"missing {UP}")
    if not CLINVAR_VS.exists(): sys.exit(f"missing {CLINVAR_VS}")
    if not RANGES_JSON.exists(): sys.exit(f"missing {RANGES_JSON}")

    df = pd.read_csv(UP, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df):,} rows from v2_scored upload")
    genes = set(df[df["source"]=="clinvar"]["gene"].unique())

    # ---------- step 1: backfill ref/alt from ClinVar bulk ----------
    print(f"\n[1/3] Parsing ClinVar bulk for (gene,hgvs_c) -> (chr,pos,ref,alt) ...")
    t0 = time.time()
    cv_coords = load_clinvar_coords(genes)
    print(f"  indexed {len(cv_coords):,} ClinVar GRCh37 entries  ({time.time()-t0:.1f}s)")

    n_filled = 0
    for i, r in df.iterrows():
        if r["source"] != "clinvar": continue
        if r["ref"] and r["alt"]: continue
        c = cv_coords.get((r["gene"], r["hgvs_c"]))
        if not c: continue
        chrom, pos, ref, alt = c
        df.at[i, "chr_grch37"] = chrom
        df.at[i, "pos_grch37"] = pos
        df.at[i, "ref"] = ref
        df.at[i, "alt"] = alt
        n_filled += 1
    print(f"  backfilled ref/alt on {n_filled:,} ClinVar rows")

    # ---------- step 2: coord-based dbNSFP retry ----------
    print(f"\n[2/3] Coord-based dbNSFP retry for non-canonical-transcript hits ...")
    score_cols = ["PhastCons100way","PhyloP100way","REVEL","CADD_phred","AlphaMissense_pred"]
    cache = {}
    n_coord = 0
    for i, r in df.iterrows():
        if any(r[c] for c in score_cols): continue
        if not r["ref"] or not r["alt"] or not r["pos_grch37"]: continue
        g = r["gene"]
        if g not in cache:
            cache[g] = load_dbnsfp_by_coord(g)
        try:
            k = (r["chr_grch37"], int(float(r["pos_grch37"])), r["ref"], r["alt"])
        except (TypeError, ValueError): continue
        hit = cache[g].get(k)
        if hit:
            for c in score_cols:
                if hit.get(c): df.at[i, c] = hit[c]
            n_coord += 1
    print(f"  scored {n_coord:,} additional rows via coord-based dbNSFP")

    # ---------- step 3: no-filter gnomAD region tabix for FAF backfill ----------
    # Build the list of (gene, chr_hg38, pos_hg38, ref, alt) we still need FAF for
    print(f"\n[3/3] Pulling gnomAD by gene region (no FAF filter) for FAF backfill ...")
    need = df[(df["source"]=="clinvar") & (df["FAF95_grpmax"]=="")
              & (df["ref"]!="") & (df["alt"]!="")].copy()
    print(f"  rows needing FAF: {len(need):,}")

    # Build hg19->hg38 map for the lookup (gnomAD is hg38)
    hg38_map = build_hg38_to_hg19_map(genes)
    # Build the inverse: (chr, hg19_pos) -> hg38_pos
    hg19_to_hg38 = {}
    for (c, hg38p), hg19_str in hg38_map.items():
        try: hg19_to_hg38[(c, int(float(hg19_str)))] = hg38p
        except (TypeError, ValueError): pass
    print(f"  hg19->hg38 map: {len(hg19_to_hg38):,} positions")

    ranges = json.loads(RANGES_JSON.read_text())

    def pull_gene(gene):
        info = ranges.get(gene)
        if not info: return gene, {}
        return gene, fetch_gnomad_region_faf(info["chr"], info["start"]-50, info["end"]+50)

    gnomad_by_gene = {}
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=6) as ex:
        futures = {ex.submit(pull_gene, g): g for g in sorted(genes) if g in ranges}
        for j, fut in enumerate(as_completed(futures), 1):
            gene, faf_map = fut.result()
            gnomad_by_gene[gene] = faf_map
            if j % 10 == 0 or j == len(futures):
                total = sum(len(m) for m in gnomad_by_gene.values())
                print(f"    [{j}/{len(futures)}] {gene:<10}  total positions cached: {total:,}")
    print(f"  region pulls done in {time.time()-t0:.0f}s")

    # Join: for each unmatched ClinVar row, lookup gnomAD by (chr_hg38, hg38_pos, ref, alt)
    n_faf = 0
    for i in need.index:
        r = df.loc[i]
        gene = r["gene"]
        gn = gnomad_by_gene.get(gene)
        if not gn: continue
        try:
            hg19_pos = int(float(r["pos_grch37"]))
        except (TypeError, ValueError): continue
        hg38p = hg19_to_hg38.get((r["chr_grch37"], hg19_pos))
        if not hg38p: continue
        faf = gn.get((r["chr_grch37"], hg38p, r["ref"], r["alt"]))
        if faf:
            df.at[i, "FAF95_grpmax"] = faf
            n_faf += 1
    print(f"  backfilled FAF on {n_faf:,} rows")

    # ---------- emit ----------
    df.to_csv(UP, sep="\t", index=False)
    print(f"\nUpdated -> {UP}")

    # Refresh per-gene
    for g in sorted(df["gene"].unique()):
        sub = df[df["gene"] == g]
        (PER_GENE_OUT / f"{g}__v2.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {df['gene'].nunique()} per-gene files")

    # Final stats
    print(f"\n=== Final score coverage (clinvar source only) ===")
    cv = df[df["source"]=="clinvar"]
    n = len(cv)
    for c in score_cols + ["SpliceAI_masked","FAF95_grpmax"]:
        good = (~cv[c].isin(["", "."])).sum()
        print(f"  {c:<22} {100*good/n:>5.1f}%  ({good:,}/{n:,})")
    # ref/alt coverage
    with_coords = ((cv["ref"]!="") & (cv["alt"]!="")).sum()
    print(f"  {'ref/alt populated':<22} {100*with_coords/n:>5.1f}%  ({with_coords:,}/{n:,})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

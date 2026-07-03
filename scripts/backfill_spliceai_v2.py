#!/usr/bin/env python3
"""Backfill measured SpliceAI masked DS_MAX onto v2 upload rows.

We generated SpliceAI Illumina-masked output for cp_new genes back in May
(`data/exports/cp_new/cp_new.spliceai.vcf` + chunked files). But the values
were never joined into the per-gene dbNSFP TSVs, so 0 rows in the v2 upload
carry real SpliceAI scores. Synonymous_catalog rows show "0 (assumed)" instead
of the measured value.

This script:
  1. Loads SpliceAI VCFs -> {(chr_hg38, pos_hg38, ref, alt): ds_max_masked}
  2. Uses our existing hg19->hg38 position map to translate upload coords
  3. Replaces empty / "0 (assumed)" SpliceAI cells with measured values
  4. Refreshes per-gene files
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd

UP = Path("cp_new_bundle/outputs/new_genes_68_v2/UPLOAD_v2_scored.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/new_genes_68_v2/per_gene")
PER_GENE_DBNSFP = Path("data/exports/cp_new")
VCF_DIR = Path("data/exports/cp_new")
VCF_PATTERNS = ("cp_new.spliceai.vcf", "cp_new.bed.chunk_*.spliceai.vcf",
                "cp_new_v2_needed.spliceai.vcf")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError): pass
    return str(v)


def parse_spliceai_vcf(path: Path) -> dict[tuple[str, int, str, str], float]:
    """Return {(chr_no_prefix, pos, ref, alt): max DS across all transcripts/scores}."""
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    out = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8: continue
            chrom, pos, _id, ref, alt, _q, _f, info = parts[:8]
            chrom = chrom.removeprefix("chr")
            try: pos_i = int(pos)
            except ValueError: continue
            # SpliceAI INFO: SpliceAI=ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            spliceai = None
            for kv in info.split(";"):
                if kv.startswith("SpliceAI="):
                    spliceai = kv[len("SpliceAI="):]
                    break
            if not spliceai: continue
            best = 0.0
            for entry in spliceai.split(","):
                fields = entry.split("|")
                if len(fields) < 6: continue
                try:
                    ds = max(float(x) for x in fields[2:6] if x not in ("", "."))
                except ValueError:
                    continue
                if ds > best:
                    best = ds
            key = (chrom, pos_i, ref, alt)
            # Keep the higher score if we see the variant twice
            if key not in out or best > out[key]:
                out[key] = best
    return out


def build_hg19_to_hg38(genes: list[str]) -> dict[tuple[str, int], int]:
    """Walk per-gene dbNSFP TSVs, build {(chr, hg19_pos): hg38_pos}."""
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
    print(f"Loaded {len(df):,} rows from v2 upload")

    # Load all SpliceAI VCFs
    vcf_files = []
    for pat in VCF_PATTERNS:
        vcf_files.extend(VCF_DIR.glob(pat))
    vcf_files = sorted(set(vcf_files))
    print(f"\nLoading SpliceAI VCFs:")
    spliceai_map = {}
    for vcf in vcf_files:
        before = len(spliceai_map)
        d = parse_spliceai_vcf(vcf)
        spliceai_map.update(d)
        print(f"  {vcf.name}: {len(d):,} variants  ({len(spliceai_map)-before:,} new)")
    print(f"  total SpliceAI entries: {len(spliceai_map):,}")

    # Build hg19 -> hg38 map
    genes = sorted(df["gene"].unique())
    print(f"\nBuilding hg19->hg38 position map from {len(genes)} dbNSFP TSVs ...")
    hg19_to_hg38 = build_hg19_to_hg38(genes)
    print(f"  {len(hg19_to_hg38):,} positions mapped")

    # Backfill
    n_set = 0
    n_replaced = 0
    n_no_coord = 0
    n_no_match = 0
    for i, r in df.iterrows():
        existing = r["SpliceAI_masked"]
        # Replace empty OR "0 (assumed)" placeholder
        if existing and existing not in ("", "0 (assumed)"):
            continue
        if not r["chr_grch37"] or not r["pos_grch37"] or not r["ref"] or not r["alt"]:
            n_no_coord += 1
            continue
        try:
            hg19_pos = int(float(r["pos_grch37"]))
        except (TypeError, ValueError):
            continue
        hg38_pos = hg19_to_hg38.get((r["chr_grch37"], hg19_pos))
        if not hg38_pos:
            n_no_match += 1
            continue
        ds = spliceai_map.get((r["chr_grch37"], hg38_pos, r["ref"], r["alt"]))
        if ds is None:
            n_no_match += 1
            continue
        df.at[i, "SpliceAI_masked"] = f"{ds:.4f}"
        if existing == "0 (assumed)": n_replaced += 1
        else: n_set += 1

    print(f"\n=== Backfill stats ===")
    print(f"  newly populated:           {n_set:,}")
    print(f"  replaced '0 (assumed)':    {n_replaced:,}")
    print(f"  no coords on row:          {n_no_coord:,}")
    print(f"  no SpliceAI match:         {n_no_match:,}")

    df.to_csv(UP, sep="\t", index=False)
    print(f"\nUpdated -> {UP}")

    # Refresh per-gene
    for g in genes:
        sub = df[df["gene"] == g]
        (PER_GENE_OUT / f"{g}__v2.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {len(genes)} per-gene files")

    # Final SpliceAI coverage
    print(f"\n=== SpliceAI coverage by source ===")
    for src in sorted(df["source"].unique()):
        sub = df[df["source"] == src]
        n = len(sub)
        non_empty = (~sub["SpliceAI_masked"].isin(["", "."])).sum()
        real = (~sub["SpliceAI_masked"].isin(["", ".", "0 (assumed)"])).sum()
        print(f"  {src:<22} {n:>6,}  populated: {non_empty:>6,} ({100*non_empty/n:.1f}%)  measured: {real:>6,} ({100*real/n:.1f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

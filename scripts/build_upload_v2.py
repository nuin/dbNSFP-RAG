#!/usr/bin/env python3
"""Build v2 of the cp_new SeqNext upload, addressing lab feedback:

  1. Independent gnomAD pull of FAF>5% variants (already done by
     pull_gnomad_common.py + reannotate_gnomad_common.py).
  2. Stop filtering synonymous coding rows -- they want them in as Benign.

Inputs:
  - data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv         (existing merge)
  - cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv (new pull)
  - cp_new_bundle/outputs/new_genes_only/_new_genes.txt          (68-gene whitelist)

Outputs (new folder cp_new_bundle/outputs/v2_gnomad_inclusive/):
  - cp_new_seqnext_UPLOAD_v2.tsv             -- merged upload (HGVSc-required)
  - cp_new_seqnext_UPLOAD_v2_dropped.tsv     -- audit trail
  - cp_new_gnomad_common_standalone.tsv      -- just the gnomAD pull, HGVSc-resolved
  - cp_new_gnomad_common_needs_hgvs.tsv      -- noncoding rows pending VV resolution
  - per_gene/{GENE}__v2.tsv                  -- per-gene split (68 files)

Filters preserved from v1:
  - 68 new genes only
  - vv_status in {ok, clinvar_assertion, local_enumeration, gnomad_common with HGVSc}
  - Drop canonical splice (-1/-2/+1/+2 nucleotide)
  - Apply intronic ROI [-15, +6]

Filters DROPPED from v1 (per lab feedback):
  - Pure synonymous coding -- now KEPT as Benign
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

NEW_GENES = Path("cp_new_bundle/outputs/new_genes_only/_new_genes.txt")
FINAL = Path("data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv")
GNOMAD = Path("cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv")
OUT_DIR = Path("cp_new_bundle/outputs/v2_gnomad_inclusive")
OUT_DIR.mkdir(exist_ok=True)
PER_GENE = OUT_DIR / "per_gene"
PER_GENE.mkdir(exist_ok=True)

OFFSET_RE = re.compile(r"c\.\d+([+\-])(\d+)")
CANON_RE = re.compile(r"c\.\d+[+\-][12][ACGT]>")


def intronic_offset(hgvsc: str):
    if not hgvsc: return None
    m = OFFSET_RE.search(hgvsc)
    if not m: return None
    return -int(m.group(2)) if m.group(1) == "-" else int(m.group(2))


def main() -> int:
    if not NEW_GENES.exists(): sys.exit(f"missing {NEW_GENES}")
    if not FINAL.exists(): sys.exit(f"missing {FINAL}")
    if not GNOMAD.exists(): sys.exit(f"missing {GNOMAD}")

    new_genes = set(NEW_GENES.read_text().split())
    final = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    gn = pd.read_csv(GNOMAD, sep="\t", dtype=str).fillna("")

    print(f"Final (all sources): {len(final):,}")
    print(f"gnomAD common pull:  {len(gn):,}")
    print()

    # Keep FINAL columns + add gnomAD common rows to match schema
    upload_cols = ["gene","transcript","hgvs_c","classification",
                   "PhastCons100way","PhyloP100way","REVEL","SpliceAI_masked",
                   "FAF95_grpmax","CADD_phred","AlphaMissense_pred","ClinVar_sig",
                   "chr_grch37","pos_grch37","ref","alt","vv_status","source"]

    # Ensure FINAL has the schema
    for c in upload_cols:
        if c not in final.columns: final[c] = ""

    # Adapt gnomAD common to schema
    gn_adapted = pd.DataFrame()
    for c in upload_cols:
        if c in gn.columns:
            gn_adapted[c] = gn[c]
        else:
            gn_adapted[c] = ""
    gn_adapted["source"] = "gnomad_common"

    # Combined pool: FINAL + gnomAD common, dedupe with two-pass keys
    combined = pd.concat([final[upload_cols], gn_adapted[upload_cols]], ignore_index=True)
    # Priority: pipeline > clinvar > synonymous_catalog > gnomad_common
    priority = {"pipeline": 0, "clinvar": 1, "synonymous_catalog": 2, "gnomad_common": 3}
    combined["_prio"] = combined["source"].map(priority).fillna(9)

    # Pass 1: dedup by (gene,chr,pos,ref,alt) for rows that have coords
    coord_key = (combined["gene"] + "|" + combined["chr_grch37"] + ":"
                 + combined["pos_grch37"] + ":" + combined["ref"] + ">" + combined["alt"])
    has_coords = (combined["chr_grch37"] != "") & (combined["pos_grch37"] != "")
    combined["_k1"] = coord_key.where(has_coords, "")

    # Pass 2: dedup by (gene,transcript,hgvs_c) -- catches ClinVar (no coords)
    # vs other sources for the same variant
    combined["_k2"] = (combined["gene"] + "|" + combined["transcript"] + "|"
                      + combined["hgvs_c"])

    combined = combined.sort_values(["_prio"])
    # Drop dups by coord key (skip empty key)
    mask1 = combined["_k1"].duplicated(keep="first") & (combined["_k1"] != "")
    combined = combined[~mask1]
    # Drop dups by hgvs key
    mask2 = combined["_k2"].duplicated(keep="first") & (combined["hgvs_c"] != "")
    combined = combined[~mask2]
    combined = combined.drop(columns=["_prio","_k1","_k2"]).reset_index(drop=True)
    print(f"Combined+deduped: {len(combined):,}")

    # 1. 68 new genes
    df = combined[combined["gene"].isin(new_genes)].copy()
    print(f"  after new-genes filter:  {len(df):,}")

    # 2. valid vv_status
    valid_status = {"ok","clinvar_assertion","local_enumeration","gnomad_common"}
    df = df[df["vv_status"].isin(valid_status) | (df["source"]=="gnomad_common")]
    # gnomad_common rows are auto-valid where HGVSc is resolved
    df = df[~((df["source"]=="gnomad_common") & (df["hgvs_c"]==""))]
    print(f"  after vv_status filter:  {len(df):,}")

    # 3. Drop canonical splice (-1,-2,+1,+2)
    df = df[~df["hgvs_c"].str.contains(CANON_RE, na=False)]
    print(f"  after canonical splice:  {len(df):,}")

    # 4. Intronic ROI [-15, +6]
    df["_off"] = df["hgvs_c"].apply(intronic_offset)
    df = df[df["_off"].isna() | df["_off"].between(-15, 6)].copy()
    print(f"  after intronic ROI:      {len(df):,}")

    # NOTE: removed the "drop pure synonymous coding" filter

    df = df.drop(columns=["_off"])

    # Write outputs
    UP_OUT = OUT_DIR / "cp_new_seqnext_UPLOAD_v2.tsv"
    df.to_csv(UP_OUT, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows -> {UP_OUT}")

    # Audit: what's new vs v1?
    v1 = pd.read_csv("cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv",
                     sep="\t", dtype=str).fillna("")
    v1_k = set(v1["gene"] + "|" + v1["chr_grch37"] + ":" + v1["pos_grch37"] + ":" + v1["ref"] + ">" + v1["alt"])
    df["_k"] = df["gene"] + "|" + df["chr_grch37"] + ":" + df["pos_grch37"] + ":" + df["ref"] + ">" + df["alt"]
    new_to_v2 = df[~df["_k"].isin(v1_k)]
    print(f"\nNet new vs v1: {len(new_to_v2):,} rows")
    print(f"  by source:")
    print("  " + new_to_v2["source"].value_counts().to_string().replace("\n", "\n  "))

    df = df.drop(columns=["_k"])

    # Standalone gnomad_common (HGVSc resolved)
    sa = gn[gn["hgvs_c"] != ""].copy()
    SA_OUT = OUT_DIR / "cp_new_gnomad_common_standalone.tsv"
    sa.to_csv(SA_OUT, sep="\t", index=False)
    print(f"\nStandalone gnomAD common (HGVSc resolved): {len(sa):,}  -> {SA_OUT}")

    # Needs HGVSc (noncoding gnomAD commons, for batch VV later)
    nv = gn[gn["hgvs_c"] == ""].copy()
    NV_OUT = OUT_DIR / "cp_new_gnomad_common_needs_hgvs.tsv"
    nv.to_csv(NV_OUT, sep="\t", index=False)
    print(f"Pending VV resolution:                       {len(nv):,}  -> {NV_OUT}")

    # Per-gene split
    for g in sorted(new_genes):
        sub = df[df["gene"] == g]
        (PER_GENE / f"{g}__v2.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"\nPer-gene files -> {PER_GENE}/  ({len(new_genes)} files)")

    # Summary stats
    print(f"\n=== v2 upload summary ===")
    def root(c):
        for k in ["Benign","Likely_benign","Pathogenic","Likely_pathogenic","Uncertain"]:
            if c.startswith(k): return k
        return "?"
    df["_root"] = df["classification"].apply(root)
    print("by source:")
    print(df["source"].value_counts().to_string())
    print("\nby class:")
    print(df["_root"].value_counts().to_string())
    print(f"\ngene coverage: {df['gene'].nunique()} / 68")
    missing = sorted(new_genes - set(df['gene'].unique()))
    if missing:
        print(f"zero-coverage genes: {missing}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

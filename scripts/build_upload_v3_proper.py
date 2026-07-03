#!/usr/bin/env python3
"""Build v3 properly: drop ClinVar BEFORE dedup so syn_catalog rows that were
previously deduped against ClinVar are restored.

Original v3 built from v2's already-deduped output. Bug: variants that existed
in both ClinVar (W1 from ClinVar evidence) and synonymous_catalog (W1 from
gnomAD FAF) were deduped to ClinVar. When v3 dropped ClinVar, those variants
disappeared even though valid syn_catalog rows existed.

Example: GOT2 c.816C>T -- syn_catalog row says "Benign FAF >5%" (FAF=0.854).
v2 deduped this in favor of ClinVar; v3 dropped ClinVar; net result the variant
is missing from v3 despite being a perfectly valid W1 call.

Fix: start from FINAL (pre-dedup), drop ClinVar source first, then dedup.

Pipeline:
  1. Load FINAL (all sources, pre-filter, 160k rows for full panel or 46k for new genes)
  2. Add gnomAD-common pull rows
  3. DROP all source=clinvar (the lab feedback)
  4. Dedup remaining: pipeline > synonymous_catalog > gnomad_common
  5. Restrict to 68 new genes
  6. Apply standard upload filters (canonical splice, intronic ROI)
  7. Re-validate W2/W3 using measured SpliceAI (masked column; unmasked re-applied later)
  8. Emit
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

NEW_GENES = Path("cp_new_bundle/outputs/_shared/new_genes.txt")
FINAL = Path("cp_new_bundle/outputs/new_genes_68_v1/FINAL.tsv")  # same data, just moved during reorg
GNOMAD = Path("cp_new_bundle/outputs/new_genes_68_v1/gnomad_common.tsv")
V3_DIR = Path("cp_new_bundle/outputs/new_genes_68_v3")
PER_GENE = V3_DIR / "per_gene"
PER_GENE.mkdir(parents=True, exist_ok=True)

W2_SPLICEAI_MAX = 0.1
W3_SPLICEAI_MAX = 0.1
OFFSET_RE = re.compile(r"c\.\d+([+\-])(\d+)")
CANON_RE = re.compile(r"c\.\d+[+\-][12][ACGT]>")


def _f(v):
    if v in ("", None, ".", "0 (assumed)"):
        return 0.0 if v == "0 (assumed)" else None
    try: return float(v)
    except (ValueError, TypeError): return None


def intronic_offset(hgvsc: str):
    if not hgvsc: return None
    m = OFFSET_RE.search(hgvsc)
    if not m: return None
    return -int(m.group(2)) if m.group(1) == "-" else int(m.group(2))


def main() -> int:
    if not FINAL.exists(): sys.exit(f"missing {FINAL}")
    if not GNOMAD.exists(): sys.exit(f"missing {GNOMAD}")
    if not NEW_GENES.exists(): sys.exit(f"missing {NEW_GENES}")

    new_genes = set(NEW_GENES.read_text().split())
    final = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    gn = pd.read_csv(GNOMAD, sep="\t", dtype=str).fillna("")
    print(f"FINAL pre-dedup: {len(final):,}")
    print(f"gnomAD common:   {len(gn):,}")
    print(f"  by source in FINAL:")
    print("  " + final["source"].value_counts().to_string().replace("\n","\n  "))

    upload_cols = ["gene","transcript","hgvs_c","classification",
                   "PhastCons100way","PhyloP100way","REVEL","SpliceAI_masked",
                   "FAF95_grpmax","CADD_phred","AlphaMissense_pred","ClinVar_sig",
                   "chr_grch37","pos_grch37","ref","alt","vv_status","source"]
    for c in upload_cols:
        if c not in final.columns: final[c] = ""
    gn_adapted = pd.DataFrame({c: gn.get(c, "") for c in upload_cols})
    gn_adapted["source"] = "gnomad_common"

    combined = pd.concat([final[upload_cols], gn_adapted[upload_cols]], ignore_index=True)
    print(f"\nCombined (pre-filter): {len(combined):,}")

    # ----- step 1: drop ClinVar BEFORE dedup -----
    n_before = len(combined)
    combined = combined[combined["source"] != "clinvar"].reset_index(drop=True)
    print(f"  dropped clinvar source: {n_before - len(combined):,}  remaining: {len(combined):,}")

    # ----- step 2: dedup (pipeline > syn_catalog > gnomad_common) -----
    priority = {"pipeline": 0, "synonymous_catalog": 1, "gnomad_common": 2}
    combined["_prio"] = combined["source"].map(priority).fillna(9)
    combined["_k1"] = (combined["gene"] + "|" + combined["chr_grch37"] + ":"
                     + combined["pos_grch37"] + ":" + combined["ref"] + ">" + combined["alt"])
    combined.loc[(combined["chr_grch37"] == "") | (combined["pos_grch37"] == ""), "_k1"] = ""
    combined["_k2"] = combined["gene"] + "|" + combined["transcript"] + "|" + combined["hgvs_c"]
    combined = combined.sort_values("_prio")
    mask1 = combined["_k1"].duplicated(keep="first") & (combined["_k1"] != "")
    combined = combined[~mask1]
    mask2 = combined["_k2"].duplicated(keep="first") & (combined["hgvs_c"] != "")
    combined = combined[~mask2]
    combined = combined.drop(columns=["_prio","_k1","_k2"]).reset_index(drop=True)
    print(f"  after dedup: {len(combined):,}")

    # ----- step 3: 68 new genes -----
    df = combined[combined["gene"].isin(new_genes)].copy()
    print(f"  after 68-gene filter: {len(df):,}")

    # ----- step 4: canonical splice drop -----
    df = df[~df["hgvs_c"].str.contains(CANON_RE, na=False)]
    print(f"  after canonical splice: {len(df):,}")

    # ----- step 5: intronic ROI [-15, +6] -----
    df["_off"] = df["hgvs_c"].apply(intronic_offset)
    df = df[df["_off"].isna() | df["_off"].between(-15, 6)].copy()
    df = df.drop(columns=["_off"])
    print(f"  after intronic ROI:    {len(df):,}")

    # ----- step 6: re-validate W2/W3 with measured masked SpliceAI -----
    def is_w2(c): return "synonymous, SpliceAI" in c
    def is_w3(c): return c.startswith("Likely_benign FAF >0.1%")

    n_before = len(df)
    keep = []
    for i, r in df.iterrows():
        sa = _f(r["SpliceAI_masked"])
        if sa is None: sa = 0.0  # treat unknown as zero (will be re-checked when unmasked applied)
        if is_w2(r["classification"]) and sa > W2_SPLICEAI_MAX: continue
        if is_w3(r["classification"]) and sa > W3_SPLICEAI_MAX: continue
        keep.append(i)
    df = df.loc[keep]
    print(f"  after W2/W3 measured-SpliceAI revalidation: {len(df):,}  "
          f"(dropped {n_before - len(df):,})")

    # ----- emit -----
    out = V3_DIR / "UPLOAD_v3.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows -> {out}")

    for g in sorted(new_genes):
        sub = df[df["gene"] == g]
        (PER_GENE / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))

    # Stats
    print(f"\n=== v3 (proper) summary ===")
    print(f"by source:")
    print(df["source"].value_counts().to_string())
    def root(c):
        for k in ["Benign","Likely_benign","Pathogenic","Likely_pathogenic","Uncertain"]:
            if c.startswith(k): return k
        return "?"
    print(f"\nby class:")
    print(df["classification"].apply(root).value_counts().to_string())
    print(f"\ngene coverage: {df['gene'].nunique()} / 68")
    missing = sorted(new_genes - set(df['gene'].unique()))
    if missing:
        print(f"zero-coverage genes ({len(missing)}): {missing}")

    # GOT2 sanity
    print(f"\n=== GOT2 lab-flagged variants in v3 (proper) ===")
    for h in ["c.816C>T","c.213T>C","c.228T>G"]:
        m = df[(df["gene"]=="GOT2") & (df["hgvs_c"]==h)]
        if len(m):
            r = m.iloc[0]
            print(f"  {h:<12} -> [{r['source']}] {r['classification']}")
        else:
            print(f"  {h:<12} -> NOT IN UPLOAD")
    return 0


if __name__ == "__main__":
    sys.exit(main())

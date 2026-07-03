#!/usr/bin/env python3
"""Build v3.2: address lab's MSH3 feedback (2026-06-18).

Two patches to v3.1:

  1. Add 11,662 high-FAF (>5%) gnomAD variants that we resolved via MyVariant
     (dels/dups, intronic, etc.) as `Benign FAF >5%`.
  2. Re-run W3 missense evaluation directly on per-gene dbNSFP TSVs using the
     fixed multi-transcript score parser. The original pipeline (May)
     mis-parsed ';'-joined REVEL values and silently dropped most missense
     variants that should have passed W3.

Input:
  - cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv (v3.1)
  - cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_resolved.tsv
  - data/exports/cp_new/{GENE}.tsv (per-gene dbNSFP for W3 catch-up)

Output:
  - overwrites UPLOAD_v3.tsv with v3.2 (more rows, no row removed)
  - per-gene files refreshed
  - audit of net-new rows: NEW_IN_V3_2.tsv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

V3 = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
RESOLVED = Path("cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_resolved.tsv")
SYN_RAW = Path("data/exports/cp_new/cp_new_synonymous_catalog.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/new_genes_68_v3/per_gene")
PER_GENE_DBNSFP = Path("data/exports/cp_new")
NEW_AUDIT = Path("cp_new_bundle/outputs/new_genes_68_v3/NEW_IN_V3_2.tsv")

W3_FAF_MIN = 0.001
W3_REVEL_MAX = 0.29
W3_SPLICEAI_MAX = 0.1


def _f(v):
    """Parse possibly-';-joined' value to float; return first real value."""
    if v is None: return None
    s = str(v).strip()
    if not s or s == "." or s == "0 (assumed)":
        return 0.0 if s == "0 (assumed)" else None
    if ";" in s:
        for p in s.split(";"):
            p = p.strip()
            if p and p != ".":
                try: return float(p)
                except: continue
        return None
    try: return float(s)
    except: return None


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except: pass
    return str(v)


def _pick_first(s: str) -> str:
    """Return first non-empty/non-dot value from ';'-joined string."""
    if not s: return ""
    if ";" not in s: return "" if s in (".",) else s
    for p in s.split(";"):
        p = p.strip()
        if p and p != ".":
            return p
    return ""


def w3_catchup(v3_df: pd.DataFrame, project_nm: dict) -> list[dict]:
    """Walk per-gene dbNSFP files, find rows that pass W3 missense and aren't
    already in v3 by (gene, hgvs_c)."""
    existing = set(zip(v3_df["gene"], v3_df["hgvs_c"]))
    new_rows = []
    for gene in sorted(v3_df["gene"].unique()):
        target_nm = project_nm.get(gene)
        if not target_nm: continue
        p = PER_GENE_DBNSFP / f"{gene}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        n_added = 0
        for _, r in df.iterrows():
            # Need a hgvs_c. Pick the FIRST hgvs from snpEff list -- the canonical
            # transcript's HGVSc. (snpEff orders by canonical first.)
            hgvs = _pick_first(_s(r.get("HGVSc_snpEff")))
            if not hgvs: continue
            if (gene, hgvs) in existing: continue
            # missense detection: must have aaref != aaalt (and not stop)
            aaref = _s(r.get("aaref")).upper()
            aaalt = _s(r.get("aaalt")).upper()
            if not aaref or not aaalt: continue
            if aaref == aaalt: continue  # synonymous
            if aaalt in ("X", "*"): continue  # nonsense
            # W3 thresholds
            faf = _f(r.get("gnomad_v41_faf95_grpmax"))
            revel = _f(r.get("REVEL_score"))
            spliceai = _f(r.get("spliceai_ds_max_masked"))
            if spliceai is None: spliceai = 0.0
            if faf is None or faf <= W3_FAF_MIN: continue
            if faf > 0.05: continue  # would be W1 not W3 (handled separately)
            if revel is None or revel >= W3_REVEL_MAX: continue
            if spliceai > W3_SPLICEAI_MAX: continue
            new_rows.append({
                "gene": gene,
                "transcript": target_nm,
                "hgvs_c": hgvs,
                "classification": f"Likely_benign FAF >0.1%, REVEL<{W3_REVEL_MAX}, SpliceAI<={W3_SPLICEAI_MAX}",
                "PhastCons100way": _pick_first(_s(r.get("phastCons100way_vertebrate"))),
                "PhyloP100way": _pick_first(_s(r.get("phyloP100way_vertebrate"))),
                "REVEL": str(revel),
                "SpliceAI_masked": _pick_first(_s(r.get("spliceai_ds_max_masked"))),
                "SpliceAI_unmasked": "",
                "FAF95_grpmax": str(faf),
                "CADD_phred": _pick_first(_s(r.get("CADD_phred"))),
                "AlphaMissense_pred": _pick_first(_s(r.get("AlphaMissense_pred"))),
                "ClinVar_sig": _s(r.get("clinvar_clnsig")),
                "chr_grch37": _s(r.get("hg19_chr")),
                "pos_grch37": _s(r.get("hg19_pos(1-based)")),
                "ref": _s(r.get("ref")),
                "alt": _s(r.get("alt")),
                "vv_status": "pipeline_w3_catchup",
                "source": "pipeline",
            })
            existing.add((gene, hgvs))
            n_added += 1
        if n_added:
            print(f"  W3 catch-up [{gene}]: +{n_added} missense LB")
    return new_rows


def main() -> int:
    if not V3.exists(): sys.exit(f"missing {V3}")
    if not RESOLVED.exists(): sys.exit(f"missing {RESOLVED}")
    v3 = pd.read_csv(V3, sep="\t", dtype=str).fillna("")
    print(f"v3.1 input: {len(v3):,} rows")

    # Build project NM_ map
    project_nm = {g: tx for g, tx in
                  v3[v3["transcript"]!=""][["gene","transcript"]].drop_duplicates().values}
    print(f"Project NM_ map: {len(project_nm)} genes")

    # --- Patch 1: add myvariant-resolved high-FAF rows ---
    res = pd.read_csv(RESOLVED, sep="\t", dtype=str).fillna("")
    print(f"\nPatch 1: myvariant-resolved high-FAF: {len(res):,} candidates")
    # Dedup against v3
    v3_keys = set(zip(v3["gene"], v3["hgvs_c"]))
    res["_existing"] = res.apply(lambda r: (r["gene"], r["hgvs_c"]) in v3_keys, axis=1)
    new_res = res[~res["_existing"]].copy()
    print(f"  net-new (not already in v3): {len(new_res):,}")
    # Normalize transcript to project NM_ where versions disagree
    for i, r in new_res.iterrows():
        target = project_nm.get(r["gene"])
        if target: new_res.at[i, "transcript"] = target
    # Drop columns not in v3 schema
    keep_cols = list(v3.columns)
    for c in keep_cols:
        if c not in new_res.columns:
            new_res[c] = ""
    new_res = new_res[keep_cols]
    new_res_dicts = new_res.to_dict("records")

    # --- Patch 2: W3 missense catch-up from per-gene dbNSFP ---
    print(f"\nPatch 2: W3 missense catch-up (re-evaluating dbNSFP with fixed score parser)")
    w3_rows = w3_catchup(v3, project_nm)
    print(f"  W3 catch-up total: +{len(w3_rows):,} missense LB")

    # --- Patch 3: W2 syn catalog catch-up ---
    # annotate_classify_synonymous.py dropped rows where PhastCons couldn't be
    # borrowed from dbNSFP. Walk the raw catalog and add any missing (gene, hgvs_c).
    # When PhastCons is unavailable, classify as Benign W2 with a note that
    # conservation wasn't measured. SpliceAI defaults to assumed-0.
    print(f"\nPatch 3: W2 syn catalog catch-up (raw catalog -> upload, no PhastCons gate)")
    existing_after_w3 = set(zip(v3["gene"], v3["hgvs_c"]))
    for r in w3_rows:
        existing_after_w3.add((r["gene"], r["hgvs_c"]))
    for r in new_res_dicts:
        existing_after_w3.add((r["gene"], r["hgvs_c"]))

    syn_cat = pd.read_csv(SYN_RAW, sep="\t", dtype=str).fillna("")
    syn_cat = syn_cat[syn_cat["gene"].isin(set(v3["gene"].unique()))]
    print(f"  syn catalog rows for new genes: {len(syn_cat):,}")

    w2_rows = []
    for _, r in syn_cat.iterrows():
        key = (r["gene"], r["hgvs_c"])
        if key in existing_after_w3: continue
        target_nm = project_nm.get(r["gene"])
        if not target_nm: continue
        if r["transcript"] != target_nm: continue  # only project canonical tx
        w2_rows.append({
            "gene": r["gene"],
            "transcript": target_nm,
            "hgvs_c": r["hgvs_c"],
            "classification": "Benign synonymous, SpliceAI<=0.1 (assumed 0), PhastCons unmeasured",
            "PhastCons100way": "",
            "PhyloP100way": "",
            "REVEL": "",
            "SpliceAI_masked": "0 (assumed)",
            "SpliceAI_unmasked": "",
            "FAF95_grpmax": "",
            "CADD_phred": "",
            "AlphaMissense_pred": "",
            "ClinVar_sig": "",
            "chr_grch37": r["chr_grch37"],
            "pos_grch37": r["pos_grch37"],
            "ref": r["ref"],
            "alt": r["alt"],
            "vv_status": "local_enumeration",
            "source": "synonymous_catalog",
        })
        existing_after_w3.add(key)
    print(f"  W2 catch-up: +{len(w2_rows):,} syn variants (will need SpliceAI backfill after)")

    # Combine
    new_combined = pd.DataFrame(new_res_dicts + w3_rows + w2_rows)
    print(f"\nAdding {len(new_combined):,} new rows to v3.1")

    v3_new = pd.concat([v3, new_combined], ignore_index=True)
    print(f"v3.2 total: {len(v3_new):,} rows")

    # Final dedup by (gene, hgvs_c) just in case
    v3_new = v3_new.drop_duplicates(subset=["gene","hgvs_c"], keep="first").reset_index(drop=True)
    print(f"after dedup: {len(v3_new):,}")

    v3_new.to_csv(V3, sep="\t", index=False)
    new_combined.to_csv(NEW_AUDIT, sep="\t", index=False)
    print(f"\nWrote -> {V3}")
    print(f"Net-new audit -> {NEW_AUDIT}")

    # Per-gene refresh
    for g in sorted(v3_new["gene"].unique()):
        sub = v3_new[v3_new["gene"] == g]
        (PER_GENE_OUT / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {v3_new['gene'].nunique()} per-gene files")

    # MSH3 sanity
    print(f"\n=== MSH3 after v3.2 ===")
    msh3 = v3_new[v3_new["gene"]=="MSH3"]
    print(f"  total MSH3 rows: {len(msh3):,}  (was {len(v3[v3['gene']=='MSH3']):,} in v3.1)")
    expected = ["c.162_179del","c.199_207del","c.359-7G>A","c.181_189del","c.178_186del",
                "c.1897-8A>G","c.178G>C","c.190C>G","c.173C>T","c.1258A>G","c.1571A>C",
                "c.1313C>T","c.205C>T","c.2740A>G","c.1522A>G","c.146C>G","c.1160T>A",
                "c.162T>C","c.1992G>A","c.204T>G","c.111C>T","c.2685C>T","c.1194C>T",
                "c.96A>C","c.3009C>T"]
    hit = sum(1 for h in expected if h in msh3["hgvs_c"].values)
    print(f"  lab's expected variants present: {hit} / {len(expected)}")
    for h in expected:
        m = msh3[msh3["hgvs_c"]==h]
        if len(m):
            r = m.iloc[0]
            print(f"   ✓ {h:<18} [{r['source']:<10}] {r['classification'][:55]}")
        else:
            print(f"   ✗ {h:<18} STILL MISSING")
    return 0


if __name__ == "__main__":
    sys.exit(main())

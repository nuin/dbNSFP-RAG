#!/usr/bin/env python3
"""Build v3.3 — gnomAD-anchored, lab feedback 2026-06-19.

Lab: "Still need more filtering. There seem to be nonsense variants in the list,
pipeline fails need to be removed, the W2 increase is the issue, anchor on gnomad."

Filters applied to v3.2:
  1. DROP all rule_fails:* rows (pipeline fails)
  2. DROP nonsense variants (aaalt = X/stop) verified against dbNSFP AA
  3. ANCHOR on gnomAD: only keep variants observed in gnomAD v4.1
     (this removes the ~84k enumerated-but-never-observed W2 rows = the W2 increase)

A variant is "observed in gnomAD" if (chr_hg38, pos_hg38, ref, alt) is in the
gnomad_observed pull. We translate the upload's hg19 coords to hg38 first.

Input:  cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv  (v3.2)
Output: cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv  (overwrites -> v3.3)
        + per-gene refresh + DROPPED audit
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# INPUT: the shared v3.2 working state (also read by build_upload_v3_4.py).
# Do NOT overwrite it -- write the v3.3 deliverable to its own folder so the
# parallel-named send folders (v3_3 / v3_4) stay in sync with the build.
V3_IN = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
V3 = Path("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv")
DROPPED = Path("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_dropped.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/new_genes_68_v3_3/per_gene")
V3.parent.mkdir(parents=True, exist_ok=True)
PER_GENE_OUT.mkdir(parents=True, exist_ok=True)
PER_GENE_DBNSFP = Path("data/exports/cp_new")
OBSERVED = Path("data/exports/cp_new/gnomad_observed/_all_observed.tsv")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except: pass
    return str(v)


def make_lifter():
    """hg19->hg38 via pyliftover. dbNSFP's coord map only covers non-synonymous
    positions, so synonymous catalog rows (which dbNSFP excludes by design)
    can't be translated from it -- they were FALSELY dropped by the gnomAD
    anchor in the first v3.3 attempt. pyliftover handles all positions."""
    from pyliftover import LiftOver
    lo = LiftOver("hg19", "hg38")
    cache = {}
    def lift(chrom, pos1):
        key = (chrom, pos1)
        if key in cache: return cache[key]
        res = lo.convert_coordinate(f"chr{chrom}", pos1 - 1)  # 1-based -> 0-based
        out = (res[0][1] + 1) if res else None  # 0-based -> 1-based
        cache[key] = out
        return out
    return lift


def build_nonsense_set(genes):
    """Set of (gene, hgvs_c) that are nonsense (aaalt = X) per dbNSFP."""
    ns = set()
    for g in genes:
        p = PER_GENE_DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        for _, r in df.iterrows():
            aaalt = _s(r.get("aaalt")).upper()
            if aaalt in ("X", "*"):
                for h in _s(r.get("HGVSc_snpEff")).split(";"):
                    h = h.strip()
                    if h and h != ".":
                        ns.add((g, h))
    return ns


def main() -> int:
    if not V3.exists(): sys.exit(f"missing {V3}")
    if not OBSERVED.exists(): sys.exit(f"missing {OBSERVED} -- run pull_gnomad_all_observed.py first")

    df = pd.read_csv(V3_IN, sep="\t", dtype=str).fillna("")
    n0 = len(df)
    print(f"v3.2 input: {n0:,} rows")
    genes = sorted(df["gene"].unique())

    dropped_frames = []

    # --- Filter 1: drop rule_fails ---
    rf = df["classification"].str.startswith("rule_fails")
    dropped_frames.append(df[rf].assign(_drop="rule_fails"))
    df = df[~rf]
    print(f"  after dropping rule_fails: {len(df):,}  (-{rf.sum():,})")

    # --- Filter 2: drop nonsense (verified via dbNSFP aaalt=X) ---
    # NOTE: key on (gene, hgvs_c) is loose -- the same c. string can appear on
    # different transcripts for two distinct genomic variants (one synonymous on
    # the canonical, one nonsense on an alternate). So we EXEMPT synonymous_catalog
    # rows: they are BioPython-verified synonymous on the project NM_ and cannot
    # be nonsense. The nonsense filter exists only to catch pipeline mis-tagging.
    print(f"  building nonsense set from dbNSFP...")
    nonsense = build_nonsense_set(genes)
    print(f"    nonsense (gene,hgvs_c) keys: {len(nonsense):,}")
    ns_mask = df.apply(lambda r: (r["source"] != "synonymous_catalog")
                       and ((r["gene"], r["hgvs_c"]) in nonsense), axis=1)
    dropped_frames.append(df[ns_mask].assign(_drop="nonsense"))
    df = df[~ns_mask]
    print(f"  after dropping nonsense: {len(df):,}  (-{ns_mask.sum():,})")

    # --- Filter 3: anchor on gnomAD observed ---
    print(f"  loading gnomAD observed set...")
    obs = pd.read_csv(OBSERVED, sep="\t", dtype=str).fillna("")
    obs_set = set(zip(obs["chr_hg38"], obs["pos_hg38"].astype(str), obs["ref"], obs["alt"]))
    print(f"    gnomAD observed variants: {len(obs_set):,}")

    print(f"  translating upload hg19 -> hg38 via pyliftover...")
    lift = make_lifter()

    def is_observed(r):
        # gnomad_common source is by definition observed.
        if r["source"] == "gnomad_common":
            return True
        chrom = r["chr_grch37"]
        if not chrom or not r["pos_grch37"] or not r["ref"] or not r["alt"]:
            return False
        try:
            hg38p = lift(chrom, int(float(r["pos_grch37"])))
        except: return False
        if hg38p is None: return False
        return (chrom, str(hg38p), r["ref"], r["alt"]) in obs_set

    obs_mask = df.apply(is_observed, axis=1)
    dropped_frames.append(df[~obs_mask].assign(_drop="not_in_gnomad"))
    df = df[obs_mask]
    print(f"  after gnomAD anchor: {len(df):,}  (-{(~obs_mask).sum():,})")

    # --- emit ---
    df.to_csv(V3, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows -> {V3}  (v3.2 {n0:,} -> v3.3 {len(df):,})")

    dropped = pd.concat(dropped_frames, ignore_index=True)
    dropped.to_csv(DROPPED, sep="\t", index=False)
    print(f"Dropped audit -> {DROPPED}  ({len(dropped):,} rows)")
    print(f"  drop reasons:")
    print("  " + dropped["_drop"].value_counts().to_string().replace("\n","\n  "))

    for g in sorted(df["gene"].unique()):
        sub = df[df["gene"] == g]
        (PER_GENE_OUT / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {df['gene'].nunique()} per-gene files")

    # Summary
    print(f"\n=== v3.3 summary ===")
    print(f"by source:")
    print(df["source"].value_counts().to_string())
    def root(c):
        if c.startswith("Benign FAF"): return "Benign W1 (FAF>5%)"
        if c.startswith("Benign synonymous"): return "Benign W2 (syn)"
        if c.startswith("Likely_benign"): return "Likely_benign W3"
        return c[:40]
    print(f"\nby classification:")
    print(df["classification"].apply(root).value_counts().to_string())
    print(f"\ngene coverage: {df['gene'].nunique()} / 68")

    # MSH3 lab variants check
    print(f"\n=== MSH3 lab variants in v3.3 ===")
    msh3 = df[df["gene"]=="MSH3"]
    print(f"  MSH3 rows: {len(msh3):,}")
    expected = ["c.162_179del","c.199_207del","c.359-7G>A","c.178_186del","c.1897-8A>G",
                "c.178G>C","c.190C>G","c.173C>T","c.1258A>G","c.1571A>C","c.1313C>T",
                "c.205C>T","c.2740A>G","c.1522A>G","c.146C>G","c.1160T>A",
                "c.162T>C","c.1992G>A","c.204T>G","c.111C>T","c.2685C>T","c.1194C>T",
                "c.96A>C","c.3009C>T"]
    hit = sum(1 for h in expected if h in msh3["hgvs_c"].values)
    print(f"  lab-expected present: {hit} / {len(expected)}")
    for h in expected:
        if h not in msh3["hgvs_c"].values:
            print(f"   ✗ lost: {h}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

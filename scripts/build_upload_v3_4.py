#!/usr/bin/env python3
"""Build v3.4 — liftover-free gnomAD anchor via (gene, NM_, c.).

The lab's point: gnomAD v4's HGVS c. is enough; no genome-build liftover needed.
Correct, with one caveat they glossed: gnomAD's displayed c. defaults to MANE
Select, which isn't always the project NM_. So we extract gnomAD's HGVSc for the
specific MANE_SELECT NM_ (pull_gnomad_hgvsc.py) and match on (gene, NM_, c.).

This is liftover-free (c. is transcript-relative, build-independent) AND
transcript-correct (matched to the project NM_, not gnomAD's default display),
and immune to the cross-transcript c.-collision that broke the position anchor
(each canonical c. is unique to its genomic variant).

Same other filters as v3.3:
  - drop rule_fails
  - drop nonsense (pipeline mis-tags only; syn_catalog exempt)

Input:  v3.2 build (run build_upload_v3.py..v3_2.py + backfill first)
Anchor: data/exports/cp_new/gnomad_hgvsc/_observed_cdot.tsv
Output: cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

V32 = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")  # state after v3.2 + spliceai
OBSERVED_CDOT = Path("data/exports/cp_new/gnomad_hgvsc/_observed_cdot.tsv")
OBSERVED_POS = Path("data/exports/cp_new/gnomad_observed/_all_observed.tsv")
OUT_DIR = Path("cp_new_bundle/outputs/new_genes_68_v3_4")
OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / "per_gene").mkdir(exist_ok=True)
(OUT_DIR / "docs").mkdir(exist_ok=True)
V34 = OUT_DIR / "UPLOAD_v3_4.tsv"
DROPPED = OUT_DIR / "UPLOAD_v3_4_dropped.tsv"
PER_GENE_DBNSFP = Path("data/exports/cp_new")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except: pass
    return str(v)


def nm_base(nm: str) -> str:
    return nm.split(".")[0] if nm else nm


def build_nonsense_set(genes):
    ns = set()
    for g in genes:
        p = PER_GENE_DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        for _, r in df.iterrows():
            if _s(r.get("aaalt")).upper() in ("X", "*"):
                for h in _s(r.get("HGVSc_snpEff")).split(";"):
                    h = h.strip()
                    if h and h != ".":
                        ns.add((g, h))
    return ns


def main() -> int:
    if not V32.exists(): sys.exit(f"missing {V32} -- run v3.2 chain first")
    if not OBSERVED_CDOT.exists():
        sys.exit(f"missing {OBSERVED_CDOT} -- run pull_gnomad_hgvsc.py first")

    df = pd.read_csv(V32, sep="\t", dtype=str).fillna("")
    n0 = len(df)
    print(f"v3.2 input: {n0:,} rows")
    genes = sorted(df["gene"].unique())
    dropped_frames = []

    # Filter 1: rule_fails
    rf = df["classification"].str.startswith("rule_fails")
    dropped_frames.append(df[rf].assign(_drop="rule_fails"))
    df = df[~rf]
    print(f"  after rule_fails: {len(df):,}  (-{rf.sum():,})")

    # Filter 2: nonsense (syn_catalog exempt -- BioPython-verified synonymous)
    nonsense = build_nonsense_set(genes)
    ns_mask = df.apply(lambda r: (r["source"] != "synonymous_catalog")
                       and ((r["gene"], r["hgvs_c"]) in nonsense), axis=1)
    dropped_frames.append(df[ns_mask].assign(_drop="nonsense"))
    df = df[~ns_mask]
    print(f"  after nonsense: {len(df):,}  (-{ns_mask.sum():,})")

    # Filter 3: gnomAD c.-based anchor (liftover-free) for genes where the
    # project NM_ == gnomAD MANE_SELECT. For the handful where they differ,
    # the c. numbering can't be reconciled, so fall back to position matching
    # (pyliftover) against the observed-coordinate set.
    obs = pd.read_csv(OBSERVED_CDOT, sep="\t", dtype=str).fillna("")
    obs_set = set(zip(obs["gene"], obs["nm"].map(nm_base), obs["cdot"]))
    print(f"  gnomAD observed (gene,NM_,c.): {len(obs_set):,}")

    # Determine per-gene whether project NM_ matches gnomAD MANE
    proj_nm = {g: nm_base(t) for g, t in
               df[df["transcript"]!=""][["gene","transcript"]].drop_duplicates().values}
    gnomad_mane = obs.groupby("gene")["nm"].agg(lambda s: set(s.map(nm_base))).to_dict()
    mismatch_genes = {g for g, pnm in proj_nm.items()
                      if pnm not in gnomad_mane.get(g, set())}
    print(f"  transcript-mismatch genes (use position fallback): {sorted(mismatch_genes)}")

    # Position fallback infra (only built if needed)
    pos_set = None
    lift = None
    if mismatch_genes and OBSERVED_POS.exists():
        posdf = pd.read_csv(OBSERVED_POS, sep="\t", dtype=str).fillna("")
        pos_set = set(zip(posdf["chr_hg38"], posdf["pos_hg38"].astype(str),
                          posdf["ref"], posdf["alt"]))
        from pyliftover import LiftOver
        _lo = LiftOver("hg19", "hg38")
        _cache = {}
        def lift(chrom, pos1):
            k = (chrom, pos1)
            if k in _cache: return _cache[k]
            res = _lo.convert_coordinate(f"chr{chrom}", pos1 - 1)
            out = (res[0][1] + 1) if res else None
            _cache[k] = out
            return out

    def observed(r):
        if r["source"] == "gnomad_common":
            return True  # came from gnomAD by definition
        g = r["gene"]
        if g in mismatch_genes:
            # position fallback via pyliftover
            if pos_set is None: return False
            if not r["chr_grch37"] or not r["pos_grch37"] or not r["ref"] or not r["alt"]:
                return False
            try:
                hg38p = lift(r["chr_grch37"], int(float(r["pos_grch37"])))
            except: return False
            if hg38p is None: return False
            return (r["chr_grch37"], str(hg38p), r["ref"], r["alt"]) in pos_set
        # c.-based anchor (liftover-free)
        return (g, nm_base(r["transcript"]), r["hgvs_c"]) in obs_set

    obs_mask = df.apply(observed, axis=1)
    dropped_frames.append(df[~obs_mask].assign(_drop="not_in_gnomad_cdot"))
    df = df[obs_mask]
    print(f"  after gnomAD c.-anchor: {len(df):,}  (-{(~obs_mask).sum():,})")

    # Emit
    df.to_csv(V34, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows -> {V34}  (v3.2 {n0:,} -> v3.4 {len(df):,})")

    dropped = pd.concat(dropped_frames, ignore_index=True)
    dropped.to_csv(DROPPED, sep="\t", index=False)
    print(f"Dropped audit -> {DROPPED}  ({len(dropped):,})")
    print("  reasons:")
    print("  " + dropped["_drop"].value_counts().to_string().replace("\n","\n  "))

    for g in sorted(df["gene"].unique()):
        (OUT_DIR / "per_gene" / f"{g}__v3_4.tsv").write_text(
            df[df["gene"]==g].to_csv(sep="\t", index=False))
    print(f"Per-gene: {df['gene'].nunique()} files")

    print(f"\n=== v3.4 summary ===")
    print(df["source"].value_counts().to_string())
    def root(c):
        if c.startswith("Benign FAF"): return "Benign W1 (FAF>5%)"
        if c.startswith("Benign synonymous"): return "Benign W2 (syn)"
        if c.startswith("Likely_benign"): return "Likely_benign W3"
        return c[:30]
    print()
    print(df["classification"].apply(root).value_counts().to_string())
    print(f"\ngene coverage: {df['gene'].nunique()} / 68")

    # MSH3 check
    msh3 = df[df["gene"]=="MSH3"]
    expected = ["c.162_179del","c.199_207del","c.359-7G>A","c.181_189del","c.178_186del",
                "c.1897-8A>G","c.178G>C","c.190C>G","c.173C>T","c.1258A>G","c.1571A>C",
                "c.1313C>T","c.205C>T","c.2740A>G","c.1522A>G","c.146C>G","c.1160T>A",
                "c.162T>C","c.1992G>A","c.204T>G","c.111C>T","c.2685C>T","c.1194C>T",
                "c.96A>C","c.3009C>T"]
    hit = sum(1 for h in expected if h in msh3["hgvs_c"].values)
    print(f"\nMSH3: {len(msh3):,} rows; lab-expected present: {hit}/{len(expected)}")
    for h in expected:
        if h not in msh3["hgvs_c"].values:
            print(f"   ✗ {h}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

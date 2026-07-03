#!/usr/bin/env python3
"""Build v3.1: fix the codon_degeneracy bug that mis-tagged missense/nonsense
variants as Benign synonymous, and drop rule_fails:no_rule_applies rows.

Lab found (2026-06-16) HOXB13 c.108C>A, c.124C>A, c.144T>A etc. wrongly tagged
"Benign synonymous". Trace: classify_cp_new.py treated codon_degeneracy in {2,4}
as 'synonymous', but 2-fold degenerate sites have only 1 of 3 alts silent --
the other 2 are missense.

Fix:
  1. For every pipeline-source row labeled "Benign synonymous", look up its real
     aaref/aaalt from the per-gene dbNSFP TSV.
  2. If aaref != aaalt (or aaalt is X = stop): NOT synonymous -> reclassify
     - if missense and passes W3: relabel as Likely_benign
     - if nonsense (aaalt=X): drop
     - otherwise: drop
  3. Also drop all rule_fails:no_rule_applies rows (lab does not want them).

Input:  cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv
Output: overwrites in place + dropped audit appended
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

V3 = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
DROPPED = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3_dropped.tsv")
PER_GENE = Path("cp_new_bundle/outputs/new_genes_68_v3/per_gene")
PER_GENE_DBNSFP = Path("data/exports/cp_new")

W3_FAF_MIN = 0.001
W3_REVEL_MAX = 0.29
SPLICEAI_MAX = 0.1


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError): pass
    return str(v)


def _f(v):
    if v in ("", None, ".", "0 (assumed)"):
        return 0.0 if v == "0 (assumed)" else None
    try: return float(v)
    except (ValueError, TypeError): return None


def load_dbnsfp_aa_by_hgvsc(gene: str) -> dict[str, tuple[str, str]]:
    """Map hgvs_c (snpEff form) -> (aaref, aaalt) for the gene."""
    p = PER_GENE_DBNSFP / f"{gene}.tsv"
    if not p.exists(): return {}
    df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
    out = {}
    for _, r in df.iterrows():
        aaref = _s(r.get("aaref")).upper()
        aaalt = _s(r.get("aaalt")).upper()
        if not aaref or not aaalt: continue
        for h in _s(r.get("HGVSc_snpEff")).split(";"):
            h = h.strip()
            if h and h != ".":
                out.setdefault(h, (aaref, aaalt))
    return out


def looks_like_misclassified_syn(classification: str, source: str) -> bool:
    """True if labeled as W2 Benign synonymous AND from pipeline source."""
    return source == "pipeline" and "synonymous" in classification.lower()


def main() -> int:
    if not V3.exists(): sys.exit(f"missing {V3}")
    df = pd.read_csv(V3, sep="\t", dtype=str).fillna("")
    n0 = len(df)
    print(f"Input v3: {n0:,} rows")

    # Per-gene AA lookup cache
    cache = {}
    def aa(gene, hgvsc):
        if gene not in cache:
            cache[gene] = load_dbnsfp_aa_by_hgvsc(gene)
        return cache[gene].get(hgvsc)

    # --- Pass 1: fix mis-tagged pipeline synonymous rows ---
    syn_pipe = df[df.apply(lambda r: looks_like_misclassified_syn(
        r["classification"], r["source"]), axis=1)]
    print(f"\nPass 1: re-checking {len(syn_pipe):,} pipeline 'Benign synonymous' rows")

    n_true_syn = 0
    n_reclassed_w3 = 0
    n_dropped_missense = 0
    n_dropped_nonsense = 0
    n_no_dbnsfp = 0

    drop_audit_rows = []
    for i, r in syn_pipe.iterrows():
        hit = aa(r["gene"], r["hgvs_c"])
        if hit is None:
            # Fallback: if REVEL or AlphaMissense are populated, this is
            # definitely missense (those predictors are missense-only).
            if r["REVEL"] or r["AlphaMissense_pred"]:
                df.at[i, "classification"] = ("rule_fails:missense_inferred_from_predictors_"
                                              "(unverified_aa)")
                n_dropped_missense += 1
                drop_audit_rows.append(dict(df.loc[i]))
                continue
            n_no_dbnsfp += 1
            continue
        aaref, aaalt = hit
        if aaref == aaalt:
            n_true_syn += 1
            continue
        # Not actually synonymous
        if aaalt == "X" or aaalt == "*":
            # Nonsense -- definitely not benign
            df.at[i, "classification"] = (f"rule_fails:nonsense_was_mis_tagged_"
                                         f"({aaref}>{aaalt})")
            n_dropped_nonsense += 1
            drop_audit_rows.append(dict(df.loc[i]))
            continue

        # Missense: try W3
        faf = _f(r["FAF95_grpmax"])
        revel = _f(r["REVEL"])
        sa = _f(r["SpliceAI_unmasked"]) if r.get("SpliceAI_unmasked") else _f(r["SpliceAI_masked"])
        if sa is None: sa = 0.0
        if (faf is not None and faf > W3_FAF_MIN
            and revel is not None and revel < W3_REVEL_MAX
            and sa <= SPLICEAI_MAX):
            df.at[i, "classification"] = (f"Likely_benign FAF >0.1%, "
                                         f"REVEL<{W3_REVEL_MAX}, "
                                         f"SpliceAI<={SPLICEAI_MAX}")
            n_reclassed_w3 += 1
        else:
            reason_bits = []
            if faf is None: reason_bits.append("faf_none")
            elif faf <= W3_FAF_MIN: reason_bits.append(f"faf_{faf:.5f}")
            if revel is None: reason_bits.append("revel_none")
            elif revel >= W3_REVEL_MAX: reason_bits.append(f"revel_{revel:.3f}")
            if sa > SPLICEAI_MAX: reason_bits.append(f"spliceai_{sa:.3f}")
            df.at[i, "classification"] = (f"rule_fails:missense_was_mis_tagged_"
                                         f"({aaref}>{aaalt})_{'|'.join(reason_bits)}")
            n_dropped_missense += 1
            drop_audit_rows.append(dict(df.loc[i]))

    print(f"  truly synonymous (kept as-is):       {n_true_syn:,}")
    print(f"  reclassified to W3 LB (missense):    {n_reclassed_w3:,}")
    print(f"  missense -> rule_fails:               {n_dropped_missense:,}")
    print(f"  nonsense -> rule_fails:               {n_dropped_nonsense:,}")
    print(f"  no dbNSFP AA lookup (kept as-is):    {n_no_dbnsfp:,}")

    # --- Pass 2: drop ALL rule_fails:no_rule_applies rows ---
    n_before = len(df)
    nra_mask = df["classification"] == "rule_fails:no_rule_applies"
    dropped_nra = df[nra_mask].copy()
    df = df[~nra_mask].reset_index(drop=True)
    print(f"\nPass 2: dropped {n_before - len(df):,} 'rule_fails:no_rule_applies' rows")

    # --- Pass 3: drop missense/nonsense rule_fails rows (no value in upload) ---
    n_before = len(df)
    misn_mask = df["classification"].str.startswith(
        ("rule_fails:missense_was_mis_tagged",
         "rule_fails:nonsense_was_mis_tagged",
         "rule_fails:missense_inferred_from_predictors"))
    dropped_mis = df[misn_mask].copy()
    df = df[~misn_mask].reset_index(drop=True)
    print(f"Pass 3: dropped {n_before - len(df):,} mis-tagged missense/nonsense rule_fails rows")

    # --- emit ---
    df.to_csv(V3, sep="\t", index=False)
    print(f"\nWrote {len(df):,} rows -> {V3}  (delta: {len(df)-n0:+,})")

    # Append to dropped audit
    if drop_audit_rows or len(dropped_nra) or len(dropped_mis):
        all_dropped = pd.concat(
            [pd.read_csv(DROPPED, sep="\t", dtype=str) if DROPPED.exists() else pd.DataFrame(),
             pd.DataFrame(drop_audit_rows) if drop_audit_rows else pd.DataFrame(),
             dropped_nra, dropped_mis],
            ignore_index=True
        )
        all_dropped.to_csv(DROPPED, sep="\t", index=False)
        print(f"Dropped audit updated -> {DROPPED}  (+{len(drop_audit_rows)+len(dropped_nra)+len(dropped_mis):,} rows)")

    # Per-gene refresh
    for g in sorted(df["gene"].unique()):
        sub = df[df["gene"] == g]
        (PER_GENE / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {df['gene'].nunique()} per-gene files")

    # Stats
    print(f"\n=== v3.1 summary ===")
    print(f"by source:")
    print(df["source"].value_counts().to_string())
    print(f"\nby classification root:")
    def root(c):
        if c.startswith("Benign"): return "Benign (rule W1/W2)"
        if c.startswith("Likely_benign"): return "Likely_benign (rule W3)"
        if c.startswith("rule_fails"): return "rule_fails (kept for review)"
        return c
    print(df["classification"].apply(root).value_counts().to_string())

    # HOXB13 sanity
    print(f"\n=== HOXB13 problem variants in v3.1 ===")
    for h in ["c.108C>A","c.108C>G","c.124C>A","c.124C>G","c.144T>A","c.144T>G"]:
        m = df[(df["gene"]=="HOXB13") & (df["hgvs_c"]==h)]
        if len(m):
            print(f"  {h:<12} -> {m.iloc[0]['classification']}")
        else:
            print(f"  {h:<12} -> DROPPED")
    return 0


if __name__ == "__main__":
    sys.exit(main())

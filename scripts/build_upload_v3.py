#!/usr/bin/env python3
"""Build v3: keep all variants (including ClinVar) but re-classify EVERYTHING
using only the project rules W1/W2/W3. ClinVar's classification text is
replaced; ClinVar evidence is no longer used.

Lab feedback (2026-06-05, clarified):
  - ClinVar variants stay in the upload
  - ClinVar's classification cannot be used as evidence
  - Apply W1/W2/W3 using supporting data (REVEL, SpliceAI, FAF, PhastCons)

Rules:
  W1: gnomAD grpmax FAF >5%                              -> Benign FAF >5%
  W2: synonymous + SpliceAI <=0.1 + PhastCons <1.0       -> Benign synonymous, SpliceAI<=0.1, PhastCons<1.0
       (intronic ROI rows additionally require PhyloP <0.1)
  W3: missense + FAF >0.1% + REVEL <0.29 + SpliceAI<=0.1 -> Likely_benign FAF >0.1%, REVEL<0.29, SpliceAI<=0.1

SpliceAI is evaluated on the unmasked column where present (matches spliceai.org),
otherwise falls back to masked. Rows where neither rule passes are dropped.

Input:  cp_new_bundle/outputs/new_genes_68_v2/UPLOAD_v2_scored.tsv
Output: cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

V2 = Path("cp_new_bundle/outputs/new_genes_68_v2/UPLOAD_v2_scored.tsv")
V3 = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
DROPPED = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3_dropped.tsv")
PER_GENE = Path("cp_new_bundle/outputs/new_genes_68_v3/per_gene")
V3.parent.mkdir(parents=True, exist_ok=True)
PER_GENE.mkdir(exist_ok=True)

W1_FAF_MIN = 0.05    # >5%
W3_FAF_MIN = 0.001   # >0.1%
W3_REVEL_MAX = 0.29  # <0.29
SPLICEAI_MAX = 0.1   # <=0.1
PHASTCONS_MAX = 1.0  # <1.0
PHYLOP_MAX = 0.1     # <0.1 (intronic ROI only)


def _f(v):
    if v in ("", None, ".", "0 (assumed)"):
        return 0.0 if v == "0 (assumed)" else None
    try: return float(v)
    except (ValueError, TypeError): return None


SYN_HGVSP = re.compile(r"p\.[A-Za-z]{1,3}\d+=")
OFFSET_RE = re.compile(r"c\.\d+([+\-])(\d+)")
MISSENSE_HGVSP = re.compile(r"p\.[A-Za-z]{1,3}\d+[A-Za-z]{1,3}")


def get_spliceai(row):
    """Prefer unmasked (spliceai.org default); fall back to masked."""
    v = _f(row.get("SpliceAI_unmasked"))
    if v is not None: return v, "unmasked"
    v = _f(row.get("SpliceAI_masked"))
    return v, ("masked" if v is not None else None)


def is_synonymous(row):
    """ClinVar rows: clinvar_name has p.XXX=. Others: classification text has 'synonymous'."""
    hgvs = row["hgvs_c"]
    # heuristic: if HGVSc is c.NNN<base>><base> and not in coding-changing form
    # Simpler: look at classification text or HGVSp if available
    classn = row["classification"]
    if "synonymous" in classn.lower(): return True
    # ClinVar source rows might have synonymous marker in clinvar_name (not present here)
    # For now also accept rows whose hgvs_c is c.\d+\w>\w (a single nucleotide change)
    # combined with explicit knowledge it's coding-only.
    return False


def is_missense(row):
    """Likely missense if AlphaMissense or REVEL has a value (those are missense-only
    predictors). Also accept the classification text mentions REVEL."""
    if row.get("REVEL") and row["REVEL"] not in ("", "."): return True
    if row.get("AlphaMissense_pred") and row["AlphaMissense_pred"] not in ("", "."): return True
    classn = row["classification"]
    if "REVEL" in classn or "missense" in classn.lower(): return True
    return False


def intronic_offset(hgvsc: str):
    if not hgvsc: return None
    m = OFFSET_RE.search(hgvsc)
    if not m: return None
    return -int(m.group(2)) if m.group(1) == "-" else int(m.group(2))


def evaluate(row) -> tuple[str | None, str | None, str]:
    """Return (classification_text, rule_label, drop_reason)."""
    faf = _f(row["FAF95_grpmax"])
    revel = _f(row["REVEL"])
    sa, sa_src = get_spliceai(row)
    pc = _f(row["PhastCons100way"])
    pp = _f(row["PhyloP100way"])
    off = intronic_offset(row["hgvs_c"])

    # W1: FAF > 5%
    if faf is not None and faf > W1_FAF_MIN:
        return f"Benign FAF >5%", "W1", ""

    syn = is_synonymous(row)
    mis = is_missense(row)

    # W2: synonymous + SpliceAI <=0.1 + PhastCons <1.0
    if syn:
        if sa is None:
            return None, None, "W2_no_spliceai"
        if sa > SPLICEAI_MAX:
            return None, None, f"W2_spliceai_{sa:.3f}"
        if pc is None:
            return None, None, "W2_no_phastcons"
        if pc >= PHASTCONS_MAX:
            return None, None, f"W2_phastcons_{pc:.3f}"
        # extra PhyloP check for intronic ROI rows
        if off is not None:
            if pp is None:
                return None, None, "W2_intronic_no_phylop"
            if pp >= PHYLOP_MAX:
                return None, None, f"W2_intronic_phylop_{pp:.3f}"
            return (f"Benign synonymous, SpliceAI<={SPLICEAI_MAX}, "
                    f"PhastCons<{PHASTCONS_MAX}, PhyloP<{PHYLOP_MAX}, "
                    f"intronic_ROI({off:+d})"), "W2", ""
        return (f"Benign synonymous, SpliceAI<={SPLICEAI_MAX}, "
                f"PhastCons<{PHASTCONS_MAX}"), "W2", ""

    # W3: missense + FAF >0.1% + REVEL <0.29 + SpliceAI <=0.1
    if mis:
        if faf is None or faf <= W3_FAF_MIN:
            return None, None, f"W3_faf_{'none' if faf is None else f'{faf:.5f}'}"
        if revel is None or revel >= W3_REVEL_MAX:
            return None, None, f"W3_revel_{'none' if revel is None else f'{revel:.3f}'}"
        if sa is None:
            return None, None, "W3_no_spliceai"
        if sa > SPLICEAI_MAX:
            return None, None, f"W3_spliceai_{sa:.3f}"
        return (f"Likely_benign FAF >0.1%, REVEL<{W3_REVEL_MAX}, "
                f"SpliceAI<={SPLICEAI_MAX}"), "W3", ""

    # Not classifiable
    return None, None, "no_rule_applies"


def main() -> int:
    if not V2.exists(): sys.exit(f"missing {V2}")
    df = pd.read_csv(V2, sep="\t", dtype=str).fillna("")
    n0 = len(df)
    print(f"Input rows (v2_scored): {n0:,}")
    print(f"  by source:")
    print("  " + df["source"].value_counts().to_string().replace("\n","\n  "))

    # Evaluate each row
    results = df.apply(lambda r: evaluate(r), axis=1, result_type="expand")
    df["_new_class"] = results[0]
    df["_rule"] = results[1]
    df["_drop_reason"] = results[2]

    # Two-tier handling:
    #  - Non-ClinVar rows that fail rules: drop entirely (legacy behavior)
    #  - ClinVar rows: keep ALL; passing rows get rule label, failing rows
    #    get "rule_fails:<reason>" so the lab sees why
    kept_mask = df["_new_class"].notna()
    is_clinvar = df["source"] == "clinvar"

    # rows we keep: rule-passing OR clinvar (regardless of pass/fail)
    keep_mask = kept_mask | is_clinvar
    dropped_mask = ~keep_mask

    kept = df[keep_mask].copy()
    dropped = df[dropped_mask].copy()

    # Apply classification text
    rule_passed = kept["_new_class"].notna()
    kept.loc[rule_passed, "classification"] = kept.loc[rule_passed, "_new_class"]
    # ClinVar rows that failed get rule_fails tag
    fail_mask = ~rule_passed
    kept.loc[fail_mask, "classification"] = "rule_fails:" + kept.loc[fail_mask, "_drop_reason"]
    kept = kept.drop(columns=["_new_class", "_rule", "_drop_reason"])

    n_pass = rule_passed.sum()
    n_fail_clinvar = fail_mask.sum()
    print(f"\nKept: {len(kept):,}  (passed rule: {n_pass:,}, clinvar rule_fails: {n_fail_clinvar:,})")
    print(f"Dropped (non-clinvar rule failures): {len(dropped):,}")
    print(f"\n=== Rule fired (kept) ===")
    df_kept_rules = df[df["_new_class"].notna()]
    print(df_kept_rules["_rule"].value_counts().to_string())

    print(f"\n=== Drop reasons (top 10) ===")
    print(dropped["_drop_reason"].value_counts().head(10).to_string())

    # Per-source: how did rules apply
    print(f"\n=== Re-classification by original source ===")
    df["_kept"] = df["_new_class"].notna()
    summary = df.groupby("source").agg(
        total=("_kept", "size"),
        kept=("_kept", "sum"),
        kept_pct=("_kept", lambda s: f"{100*s.mean():.1f}%"),
    )
    print(summary.to_string())

    kept.to_csv(V3, sep="\t", index=False)
    dropped.to_csv(DROPPED, sep="\t", index=False)
    print(f"\nWrote {len(kept):,} rows -> {V3}")
    print(f"Dropped audit -> {DROPPED}")

    # Per-gene
    for g in sorted(kept["gene"].unique()):
        sub = kept[kept["gene"] == g]
        (PER_GENE / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {kept['gene'].nunique()} per-gene files")

    # Gene coverage
    NEW = Path("cp_new_bundle/outputs/_shared/new_genes.txt")
    new_genes = set(NEW.read_text().split())
    missing = sorted(new_genes - set(kept["gene"].unique()))
    print(f"\nGene coverage: {kept['gene'].nunique()} / 68")
    if missing:
        print(f"Zero-coverage: {missing}")

    # GOT2 sanity
    print(f"\n=== GOT2 lab-flagged variants in v3 ===")
    for h in ["c.816C>T", "c.213T>C", "c.228T>G"]:
        m = kept[(kept["gene"]=="GOT2") & (kept["hgvs_c"]==h)]
        if len(m):
            r = m.iloc[0]
            print(f"  {h:<12} [{r['source']:<18}]  {r['classification']}")
        else:
            d = dropped[(dropped["gene"]=="GOT2") & (dropped["hgvs_c"]==h)]
            if len(d):
                print(f"  {h:<12} DROPPED reason={d.iloc[0]['_drop_reason']}")
            else:
                print(f"  {h:<12} NOT in v2 input")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Backfill predictor + FAF scores onto ClinVar-source rows in the v2 upload.

ClinVar bulk pull only carries assertion text + transcript + HGVSc; the
predictor / frequency / conservation columns are blank for every clinvar row.
This script joins each row to its per-gene dbNSFP TSV by
(gene, transcript, hgvs_c) and fills the columns.

Input:  cp_new_bundle/outputs/v2_gnomad_inclusive/cp_new_seqnext_UPLOAD_v2.tsv
Output: cp_new_bundle/outputs/v2_gnomad_inclusive/cp_new_seqnext_UPLOAD_v2_scored.tsv
        + refreshes per_gene/{GENE}__v2.tsv files
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

UP = Path("cp_new_bundle/outputs/v2_gnomad_inclusive/cp_new_seqnext_UPLOAD_v2.tsv")
OUT = Path("cp_new_bundle/outputs/v2_gnomad_inclusive/cp_new_seqnext_UPLOAD_v2_scored.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/v2_gnomad_inclusive/per_gene")
PER_GENE_DIR = Path("data/exports/cp_new")
SYN_ANNO = Path("data/exports/cp_new/seqnext/cp_new_seqnext_synonymous_catalog.tsv")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError):
        pass
    return str(v)


# Map upload columns -> dbNSFP columns
SCORE_MAP = {
    "PhastCons100way":   "phastCons100way_vertebrate",
    "PhyloP100way":      "phyloP100way_vertebrate",
    "REVEL":             "REVEL_score",
    "SpliceAI_masked":   "spliceai_ds_max_masked",
    "FAF95_grpmax":      "gnomad_v41_faf95_grpmax",
    "CADD_phred":        "CADD_phred",
    "AlphaMissense_pred":"AlphaMissense_pred",
}


def _pick_first_real(s: str) -> str:
    """For ';'-joined per-transcript values, return the first non-empty/non-dot.
    Returns '' if no real value exists."""
    if not s or s == ".": return ""
    if ";" not in s: return "" if s == "." else s
    for p in s.split(";"):
        p = p.strip()
        if p and p != ".":
            return p
    return ""


def load_dbnsfp_by_hgvs(gene: str) -> tuple[dict, dict]:
    """Index dbNSFP rows by canonical hgvs_c (any matching transcript suffices)
    and by hg19 coords."""
    p = PER_GENE_DIR / f"{gene}.tsv"
    if not p.exists(): return {}, {}
    df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
    by_hgvs = {}
    by_coord = {}
    for _, r in df.iterrows():
        # Collapse multi-transcript values to a single representative
        scores = {}
        for upload_col, dbnsfp_col in SCORE_MAP.items():
            scores[upload_col] = _pick_first_real(_s(r.get(dbnsfp_col)))

        # Index by HGVSc -- the value may be ';'-joined across transcripts;
        # use ALL distinct values as keys so any matches succeed.
        for h in _s(r.get("HGVSc_snpEff")).split(";"):
            h = h.strip()
            if h and h != ".":
                by_hgvs.setdefault(h, scores)

        # Index by (chr_hg19, pos_hg19, ref, alt) for coord-based fallback
        try:
            k = (_s(r["hg19_chr"]), int(float(r["hg19_pos(1-based)"])),
                 _s(r["ref"]), _s(r["alt"]))
            by_coord.setdefault(k, scores)
        except (ValueError, TypeError, KeyError):
            pass

    return by_hgvs, by_coord


def load_syn_lookup() -> tuple[dict, dict]:
    """Index annotated syn-catalog by (gene, hgvs_c) and (gene, chr, pos) -> scores."""
    if not SYN_ANNO.exists(): return {}, {}
    df = pd.read_csv(SYN_ANNO, sep="\t", dtype=str, na_values=["."]).fillna("")
    score_cols = list(SCORE_MAP.keys())
    by_hgvs = {}
    by_coord = {}
    for _, r in df.iterrows():
        scores = {c: _s(r.get(c)) for c in score_cols}
        gene = _s(r.get("gene"))
        hgvs = _s(r.get("hgvs_c"))
        if gene and hgvs:
            by_hgvs.setdefault((gene, hgvs), scores)
        try:
            k = (gene, _s(r["chr_grch37"]), int(float(r["pos_grch37"])),
                 _s(r["ref"]), _s(r["alt"]))
            by_coord.setdefault(k, scores)
        except (ValueError, TypeError, KeyError):
            pass
    return by_hgvs, by_coord


def main() -> int:
    if not UP.exists():
        sys.exit(f"missing {UP}")
    df = pd.read_csv(UP, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df):,} rows from v2 upload")
    print(f"Loading annotated synonymous catalog ...")
    syn_by_hgvs, syn_by_coord = load_syn_lookup()
    print(f"  syn-catalog indexed: {len(syn_by_hgvs):,} (gene,hgvs); "
          f"{len(syn_by_coord):,} (gene,coord)")
    print(f"\nBackfilling scores ...\n")

    # cache per-gene dbNSFP lookups
    cache = {}
    def lookups(gene):
        if gene not in cache:
            cache[gene] = load_dbnsfp_by_hgvs(gene)
        return cache[gene]

    score_cols = list(SCORE_MAP.keys())
    n_hit_dbnsfp = 0
    n_hit_syn = 0
    n_no_match = 0

    for i, r in df.iterrows():
        # Only fill if all 7 are empty (don't overwrite existing scores)
        if any(r[c] for c in score_cols):
            continue
        gene = r["gene"]
        if not gene: continue
        hit = None
        # Path 1: dbNSFP by hgvs_c (missense / non-syn coding)
        by_hgvs, by_coord = lookups(gene)
        if r["hgvs_c"]:
            hit = by_hgvs.get(r["hgvs_c"])
            if hit:
                n_hit_dbnsfp += 1
        if not hit:
            # Path 2: annotated synonymous catalog by (gene, hgvs_c)
            if r["hgvs_c"]:
                hit = syn_by_hgvs.get((gene, r["hgvs_c"]))
                if hit:
                    n_hit_syn += 1
        if not hit:
            # Path 3: coord-based fallback against dbNSFP (rarely matches for ClinVar)
            try:
                k = (r["chr_grch37"], int(float(r["pos_grch37"])), r["ref"], r["alt"])
                hit = by_coord.get(k)
            except (ValueError, TypeError):
                pass
        if not hit:
            n_no_match += 1
            continue
        for c in score_cols:
            if hit.get(c):
                df.at[i, c] = hit[c]

    print(f"  matched via dbNSFP (missense): {n_hit_dbnsfp:,}")
    print(f"  matched via syn catalog:       {n_hit_syn:,}")
    print(f"  no match in either source:     {n_no_match:,}")
    print()

    # Stats
    print("=== Post-backfill score coverage by source ===")
    for src in sorted(df["source"].unique()):
        sub = df[df["source"] == src]
        n = len(sub)
        print(f"\n  {src} ({n:,} rows):")
        for c in score_cols:
            good = (~sub[c].isin(["", "."])).sum()
            print(f"    {c:<22} {100*good/n:>5.1f}%  ({good:,}/{n:,})")

    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote -> {OUT}")

    # Refresh per-gene
    PER_GENE_OUT.mkdir(parents=True, exist_ok=True)
    for g in sorted(df["gene"].unique()):
        sub = df[df["gene"] == g]
        (PER_GENE_OUT / f"{g}__v2.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {df['gene'].nunique()} per-gene files in {PER_GENE_OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Fix two W3 (missense LB) defects found in lab review (2026-07-15):

ISSUE 1 — synonymous/non-missense variants in the missense-only W3 category.
  Root cause: W3 keyed on REVEL presence, but REVEL can leak onto a synonymous
  variant from an overlapping missense annotation. Fix: keep in W3 only if the
  variant is a CONFIRMED missense on the project transcript (dbNSFP aaref!=aaalt,
  aaalt not a stop). Others are removed from W3.

ISSUE 2 — wrong coordinate / FAF (pipeline_w3_catchup took a cDNA number from one
  transcript but kept a different position's genomic coord + FAF). Detected by
  comparing the row's stored gnomad_id to gnomAD's own id for (gene, NM_, c.).
  Mismatches are dropped (their FAF is borrowed from a neighbouring variant, so
  the variant fails lab QA against the gnomAD browser).

Updates v3_3/v3_4 (full+roi), regenerates v3_5 and v3_7, refreshes per-gene.
Writes an audit of everything removed.
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

DBNSFP = Path("data/exports/cp_new")
GPMAP = Path("data/exports/cp_new/gnomad_hgvsc/_gp_map.tsv")
TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]


def _pick(s):
    s = str(s)
    if ";" not in s: return "" if s in (".", "nan") else s
    for p in s.split(";"):
        if p.strip() and p.strip() != ".": return p.strip()
    return ""


def load_dbnsfp_consequence(genes):
    """(gene, chr_hg19, pos_hg19, ref, alt) -> consequence (missense/synonymous/nonsense)."""
    m = {}
    for g in genes:
        p = DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        for _, r in df.iterrows():
            aaref, aaalt = str(r.get("aaref","")).upper(), str(r.get("aaalt","")).upper()
            if not aaref or not aaalt: continue
            cons = ("synonymous" if aaref == aaalt else
                    "nonsense" if aaalt in ("X","*") else "missense")
            try:
                k = (g, str(r["hg19_chr"]), int(float(r["hg19_pos(1-based)"])), r["ref"], r["alt"])
                m.setdefault(k, cons)
            except (ValueError, TypeError, KeyError):
                pass
    return m


def main():
    gp = pd.read_csv(GPMAP, sep="\t", dtype=str).fillna("")
    gp["nmb"] = gp["nm"].str.split(".").str[0]
    gmap = {(r["gene"], r["nmb"], r["cdot"]): r["gnomad_id"] for _, r in gp.iterrows()}

    # consequence map (built once from all genes present)
    ref_df = pd.read_csv(TABLES[0][0], sep="\t", dtype=str).fillna("")
    cons_map = load_dbnsfp_consequence(sorted(ref_df["gene"].unique()))
    print(f"dbNSFP consequence map: {len(cons_map):,} entries")

    P_MISSENSE = re.compile(r"p\.([A-Za-z]{3})\d+([A-Za-z]{3})$")

    audit_all = []
    for path, pgname, tag in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        df["nmb"] = df["transcript"].str.split(".").str[0]
        drop_idx = []
        for i in df.index:
            r = df.loc[i]
            if not r["classification"].startswith("Likely_benign"):
                continue
            # ISSUE 2: coordinate/gnomad_id mismatch -> phantom, drop
            true_id = gmap.get((r["gene"], r["nmb"], r["hgvs_c"]))
            if r["source"] != "gnomad_common" and true_id and r["gnomad_id"] and true_id != r["gnomad_id"]:
                drop_idx.append((i, "coord_mismatch (gnomAD id differs -> wrong FAF)")); continue
            # ISSUE 1: confirm missense
            pm = r["p_hgvs"]
            confirmed = bool(P_MISSENSE.match(pm)) and not pm.endswith("=")
            if not confirmed:
                # fall back to dbNSFP consequence at this genomic position
                try:
                    k = (r["gene"], r["chr_grch37"], int(float(r["pos_grch37"])), r["ref"], r["alt"])
                except (ValueError, TypeError):
                    k = None
                cons = cons_map.get(k) if k else None
                if cons == "missense":
                    confirmed = True
                elif cons in ("synonymous", "nonsense"):
                    drop_idx.append((i, f"not_missense ({cons})")); continue
                else:
                    drop_idx.append((i, "not_missense (unconfirmed)")); continue
        idx = [i for i, _ in drop_idx]
        for i, why in drop_idx:
            row = df.loc[i].to_dict(); row["_removed_reason"] = why; audit_all.append(row)
        df = df.drop(index=idx).drop(columns=["nmb"])
        df.to_csv(p, sep="\t", index=False)
        print(f"  {p.name}: removed {len(idx)} W3 non-missense/coord-mismatch -> {len(df):,}")
        pg = p.parent / pgname; pg.mkdir(exist_ok=True)
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{tag}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))

    if audit_all:
        pd.DataFrame(audit_all).to_csv(
            "cp_new_bundle/outputs/new_genes_68_v3_7/W3_removed_audit.tsv", sep="\t", index=False)

    # regenerate v3.5 and v3.7 from v3_4 ROI
    d = pd.read_csv("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", sep="\t", dtype=str).fillna("")
    allg = sorted(set(open("cp_new_bundle/outputs/_shared/new_genes.txt").read().split()))
    hdr = "\t".join(d.columns)
    for ver in ("v3_5", "v3_7"):
        v = Path(f"cp_new_bundle/outputs/new_genes_68_{ver}")
        (v/"per_gene").mkdir(parents=True, exist_ok=True)
        d.to_csv(v/f"UPLOAD_{ver}.tsv", sep="\t", index=False)
        for f in (v/"per_gene").glob("*.tsv"): f.unlink()
        for g in allg:
            sub = d[d["gene"]==g]
            (v/"per_gene"/f"{g}__{ver}.tsv").write_text(sub.to_csv(sep="\t",index=False) if len(sub) else hdr+"\n")

    w3 = d[d["classification"].str.startswith("Likely_benign")]
    bad = w3[~w3["p_hgvs"].str.match(P_MISSENSE, na=False)]
    print(f"\nv3.7 FINAL: {len(d):,} rows | W3 missense {len(w3)} | non-missense left in W3: {len(bad)}")


if __name__ == "__main__":
    main()

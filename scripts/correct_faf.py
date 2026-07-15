#!/usr/bin/env python3
"""Correct FAF95_grpmax against authoritative gnomAD v4.1 and re-apply the
frequency gates (lab review 2026-07-15).

The ClinVar-row FAF backfill (pass 3) matched by lifted position and sometimes
grabbed a NEIGHBOURING variant's FAF (e.g. RNF43 c.1007G>A stored 0.0026 vs real
9e-05). This resets every row's FAF to gnomAD's own fafmax_faf95_max at that
variant's gnomad_id (from the observed pull), then re-gates:
  W1 (Benign FAF>5%): drop if real FAF <= 0.05
  W3 (LB missense, FAF>0.1%): drop if real FAF <= 0.001
W2 synonymous are frequency-independent (kept; FAF is informational).

gnomad_common rows keep their own pull FAF. Updates all tables, regenerates
v3_5/v3_7, writes an audit.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

OBS = Path("data/exports/cp_new/gnomad_observed/_all_observed.tsv")
TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]
W1_MIN, W3_MIN = 0.05, 0.001


def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def main():
    obs = pd.read_csv(OBS, sep="\t", dtype=str).fillna("")
    faf = {(r["chr_hg38"], r["pos_hg38"], r["ref"], r["alt"]): r["faf95"] for _, r in obs.iterrows()}
    print(f"gnomAD observed FAF index: {len(faf):,}")

    audit = []
    for path, pgname, tag in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        n_corr = n_dropw1 = n_dropw3 = 0
        drop_idx = []
        for i in df.index:
            r = df.loc[i]
            parts = r["gnomad_id"].split("-")
            real = faf.get((parts[0], parts[1], parts[2], parts[3])) if len(parts) == 4 else None
            if real is None:
                # not in observed pull -> treat as ~0 for gating (variant essentially absent)
                real_v = 0.0
            else:
                real_v = _num(real) or 0.0
                if r["FAF95_grpmax"] != real and abs((_num(r["FAF95_grpmax"]) or 0.0) - real_v) > 1e-9:
                    df.at[i, "FAF95_grpmax"] = real  # authoritative gnomAD value
                    n_corr += 1
            c = r["classification"]
            if c.startswith("Benign FAF") and real_v <= W1_MIN:
                drop_idx.append((i, f"W1_faf_now_{real_v:.5f}<=0.05")); n_dropw1 += 1
            elif c.startswith("Likely_benign") and real_v <= W3_MIN:
                drop_idx.append((i, f"W3_faf_now_{real_v:.5f}<=0.001")); n_dropw3 += 1
        for i, why in drop_idx:
            row = df.loc[i].to_dict(); row["_removed_reason"] = why; audit.append(row)
        df = df.drop(index=[i for i, _ in drop_idx])
        df.to_csv(p, sep="\t", index=False)
        print(f"  {p.name}: FAF corrected {n_corr} | dropped W1 {n_dropw1}, W3 {n_dropw3} -> {len(df):,}")
        pg = p.parent / pgname; pg.mkdir(exist_ok=True)
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{tag}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))

    if audit:
        pd.DataFrame(audit).to_csv("cp_new_bundle/outputs/new_genes_68_v3_7/FAF_regate_audit.tsv",
                                   sep="\t", index=False)

    d = pd.read_csv("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", sep="\t", dtype=str).fillna("")
    allg = sorted(set(open("cp_new_bundle/outputs/_shared/new_genes.txt").read().split()))
    hdr = "\t".join(d.columns)
    for ver in ("v3_5", "v3_7"):
        v = Path(f"cp_new_bundle/outputs/new_genes_68_{ver}")
        d.to_csv(v/f"UPLOAD_{ver}.tsv", sep="\t", index=False)
        for f in (v/"per_gene").glob("*.tsv"): f.unlink()
        for g in allg:
            sub = d[d["gene"]==g]
            (v/"per_gene"/f"{g}__{ver}.tsv").write_text(sub.to_csv(sep="\t",index=False) if len(sub) else hdr+"\n")
    print(f"\nv3.7 FINAL: {len(d):,} rows")


if __name__ == "__main__":
    main()

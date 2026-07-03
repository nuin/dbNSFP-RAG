#!/usr/bin/env python3
"""Re-annotate the gnomAD common-variants file with HGVSc from synonymous
catalog + dbNSFP, after building a hg38->hg19 position-level map from all
per-gene dbNSFP TSVs.

The initial pull (scripts/pull_gnomad_common.py) only resolved 0.7% of
variants because the syn-catalog is indexed by hg19 but gnomAD gives hg38.
This script bridges the two coordinate systems.

Output: overwrites cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv
in place (after backing up to .raw.tsv).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

GNOMAD_TSV = Path("cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv")
SYN_CATALOG = Path("data/exports/cp_new/cp_new_synonymous_catalog.tsv")
PER_GENE_DIR = Path("data/exports/cp_new")


def _s(v):
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError):
        pass
    return str(v)


def build_hg38_to_hg19_map(genes: list[str]) -> dict[tuple, str]:
    """Walk all per-gene dbNSFP TSVs once. Index (chr, hg38_pos) -> hg19_pos."""
    m = {}
    for g in genes:
        p = PER_GENE_DIR / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        if "#chr" not in df.columns: continue
        for _, r in df.iterrows():
            try:
                k = (_s(r["#chr"]), int(float(r["pos(1-based)"])))
                h19 = _s(r.get("hg19_pos(1-based)"))
                if h19 and k not in m:
                    m[k] = h19
            except (TypeError, ValueError, KeyError):
                continue
    return m


def load_syn_catalog() -> dict[tuple, dict]:
    """Index syn catalog by (chr_grch37, pos_grch37, ref, alt) -> HGVSc/transcript."""
    if not SYN_CATALOG.exists(): return {}
    df = pd.read_csv(SYN_CATALOG, sep="\t", dtype=str, na_values=["."]).fillna("")
    out = {}
    for _, r in df.iterrows():
        try:
            k = (_s(r["chr_grch37"]), int(float(r["pos_grch37"])), _s(r["ref"]), _s(r["alt"]))
        except (TypeError, ValueError, KeyError):
            continue
        if k in out: continue
        out[k] = {
            "transcript": _s(r.get("transcript")),
            "hgvs_c": _s(r.get("hgvs_c")),
        }
    return out


def load_dbnsfp_lookup(gene: str) -> dict[tuple, dict]:
    p = PER_GENE_DIR / f"{gene}.tsv"
    if not p.exists(): return {}
    df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
    out = {}
    for _, r in df.iterrows():
        try:
            k = (_s(r["#chr"]), int(float(r["pos(1-based)"])), _s(r["ref"]), _s(r["alt"]))
        except (TypeError, ValueError, KeyError):
            continue
        if k in out: continue
        aaref, aaalt = _s(r.get("aaref")), _s(r.get("aaalt"))
        out[k] = {
            "transcript": _s(r.get("Ensembl_transcriptid")).split(";")[0],
            "hgvs_c": _s(r.get("HGVSc_snpEff")).split(";")[0],
            "consequence": "missense" if (aaref and aaalt and aaref != aaalt) else "coding_other",
            "chr_hg19": _s(r.get("hg19_chr")),
            "pos_hg19": _s(r.get("hg19_pos(1-based)")),
        }
    return out


def main() -> int:
    if not GNOMAD_TSV.exists():
        sys.exit(f"missing {GNOMAD_TSV}")
    df = pd.read_csv(GNOMAD_TSV, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df):,} gnomAD common rows")

    # Backup
    backup = GNOMAD_TSV.with_suffix(".raw.tsv")
    if not backup.exists():
        backup.write_text(GNOMAD_TSV.read_text())
        print(f"backup -> {backup}")

    genes = sorted(df["gene"].unique())
    print(f"Building hg38->hg19 position map from {len(genes)} per-gene dbNSFP TSVs ...")
    hg38_to_hg19 = build_hg38_to_hg19_map(genes)
    print(f"  hg38->hg19 entries: {len(hg38_to_hg19):,}")

    print(f"Loading syn catalog ...")
    syn = load_syn_catalog()
    print(f"  syn-catalog entries: {len(syn):,}")

    # cache dbNSFP per gene
    dbnsfp_cache: dict[str, dict] = {}
    def dbn(g):
        if g not in dbnsfp_cache:
            dbnsfp_cache[g] = load_dbnsfp_lookup(g)
        return dbnsfp_cache[g]

    out_rows = []
    n_syn, n_db, n_unann = 0, 0, 0
    for _, r in df.iterrows():
        gene = r["gene"]
        chr38 = r["chr_hg38"]
        try:
            pos38 = int(float(r["pos_hg38"]))
        except (ValueError, TypeError):
            continue
        ref, alt = r["ref"], r["alt"]

        # path 1: dbNSFP hit by (chr_hg38, pos_hg38, ref, alt) -- missense
        k_full = (chr38, pos38, ref, alt)
        d = dbn(gene).get(k_full)
        transcript = ""
        hgvs_c = ""
        chr19 = ""
        pos19 = ""
        consequence = ""
        anno_source = ""
        if d:
            transcript = d["transcript"]
            hgvs_c = d["hgvs_c"]
            consequence = d["consequence"]
            chr19 = d["chr_hg19"]
            pos19 = d["pos_hg19"]
            anno_source = "dbnsfp"
            n_db += 1
        else:
            # path 2: lookup hg19 pos via hg38->hg19 map, then syn-catalog
            hg19 = hg38_to_hg19.get((chr38, pos38))
            if hg19:
                pos19 = hg19
                # syn catalog uses hg19 chr (same as hg38 chr for autosomes)
                chr19 = chr38
                try:
                    syn_hit = syn.get((chr19, int(float(hg19)), ref, alt))
                except (TypeError, ValueError):
                    syn_hit = None
                if syn_hit:
                    transcript = syn_hit["transcript"]
                    hgvs_c = syn_hit["hgvs_c"]
                    consequence = "synonymous"
                    anno_source = "syn_catalog"
                    n_syn += 1
                else:
                    consequence = "coding_unannotated"
                    anno_source = "hg38_to_hg19_only"
                    n_unann += 1
            else:
                # Outside coding; no hg19 mapping from dbNSFP available
                consequence = "noncoding"
                anno_source = "none"
                n_unann += 1

        out_rows.append({
            "gene": gene,
            "transcript": transcript,
            "hgvs_c": hgvs_c,
            "classification": f"Benign FAF >5% (gnomAD v4.1, {consequence or 'unannotated'})",
            "chr_grch37": chr19,
            "pos_grch37": pos19,
            "chr_hg38": chr38,
            "pos_hg38": pos38,
            "ref": ref,
            "alt": alt,
            "FAF95_grpmax": r["FAF95_grpmax"],
            "gnomad_v41_af_joint": r["gnomad_v41_af_joint"],
            "gnomad_v41_grpmax_anc": r["gnomad_v41_grpmax_anc"],
            "gnomad_v41_nhomalt": r["gnomad_v41_nhomalt"],
            "gnomad_v41_filter": r["gnomad_v41_filter"],
            "consequence": consequence,
            "anno_source": anno_source,
            "vv_status": "ok" if hgvs_c else "needs_hgvs",
            "source": "gnomad_common",
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(GNOMAD_TSV, sep="\t", index=False)

    print(f"\n=== Re-annotation summary ===")
    print(f"Total rows: {len(out_df):,}")
    print(f"  dbNSFP-matched (missense/coding_other): {n_db:,}")
    print(f"  syn-catalog-matched (synonymous):        {n_syn:,}")
    print(f"  unannotated:                             {n_unann:,}")
    print(f"  HGVSc resolved: {(out_df['hgvs_c'] != '').sum():,}  "
          f"({100*(out_df['hgvs_c'] != '').sum()/len(out_df):.1f}%)")
    print(f"\nby consequence:")
    print(out_df["consequence"].value_counts().to_string())
    print(f"\nby anno_source:")
    print(out_df["anno_source"].value_counts().to_string())
    print(f"\nFile updated: {GNOMAD_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

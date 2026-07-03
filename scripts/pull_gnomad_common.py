#!/usr/bin/env python3
"""Pull common variants (grpmax FAF95 > 5%) from gnomAD v4.1 joint sites for
each cp_new gene, independent of dbNSFP/SpliceAI coverage.

Why: the W1 (FAF>5%) rule in classify_cp_new.py runs after a dbNSFP join.
dbNSFP 5.3 is non-synonymous-only, so synonymous common variants never reach
W1. Lab reported missed >5% synonymous variants; this fills the gap.

Approach: tabix gnomAD v4.1 over each gene's hg38 range, filter to grpmax
FAF95 >= 0.05, join back to:
  - dbNSFP per-gene TSVs to recover HGVSc on missense
  - synonymous catalog to recover HGVSc on synonymous coding
Unannotated rows (intronic/UTR/etc) are kept with empty HGVSc and flagged
needs_hgvs so VariantValidator can be run later.

Output:
  - data/exports/cp_new/gnomad_common/{GENE}.tsv  -- per-gene rows
  - cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv -- combined
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

GNOMAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"
)
DEFAULT_THRESHOLD = 0.05  # 5%
RANGES_JSON = Path("cp_new_bundle/outputs/new_genes_only/_gene_hg38_ranges.json")
SYN_CATALOG = Path("data/exports/cp_new/cp_new_synonymous_catalog.tsv")
PER_GENE_DIR = Path("data/exports/cp_new")
OUT_PER_GENE = Path("data/exports/cp_new/gnomad_common")
OUT_COMBINED = Path("cp_new_bundle/outputs/new_genes_only/cp_new_gnomad_common.tsv")


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def tabix_region(url: str, chrom: str, start: int, end: int) -> list[str]:
    region = f"chr{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ["tabix", url, region],
            capture_output=True, text=True, timeout=900, check=False,
        )
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {region}", file=sys.stderr)
        return []
    if result.returncode != 0:
        print(f"  tabix err {region}: {result.stderr.strip()[:200]}", file=sys.stderr)
        return []
    return [l for l in result.stdout.splitlines() if l and not l.startswith("#")]


def fetch_gene_common(gene: str, chrom: str, start: int, end: int,
                      threshold: float) -> pd.DataFrame:
    """Tabix gnomAD region, filter to FAF95 grpmax >= threshold."""
    url = GNOMAD_URL.format(chrom=chrom)
    rows = []
    for line in tabix_region(url, chrom, start, end):
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        c, pos, _id, ref, alt, _q, vcf_filter, info = parts[:8]
        ifields = parse_info(info)
        faf = ifields.get("fafmax_faf95_max_joint") or ifields.get("fafmax_faf95_max")
        if not faf or faf == ".":
            continue
        try:
            faf_v = float(faf)
        except ValueError:
            continue
        if faf_v < threshold:
            continue
        rows.append({
            "gene": gene,
            "chr_hg38": c.removeprefix("chr"),
            "pos_hg38": int(pos),
            "ref": ref,
            "alt": alt,
            "gnomad_v41_faf95_grpmax": faf_v,
            "gnomad_v41_af_joint": ifields.get("AF_joint") or ifields.get("AF"),
            "gnomad_v41_grpmax_anc": ifields.get("fafmax_faf95_max_gen_anc_joint")
                or ifields.get("grpmax_joint"),
            "gnomad_v41_nhomalt": ifields.get("nhomalt_joint") or ifields.get("nhomalt"),
            "gnomad_v41_filter": vcf_filter,
        })
    return pd.DataFrame(rows)


def _s(v) -> str:
    """NaN-safe string conversion (pandas NaN -> '')."""
    if v is None: return ""
    try:
        if pd.isna(v): return ""
    except (TypeError, ValueError):
        pass
    return str(v)


def load_dbnsfp_lookup(gene: str) -> dict:
    """index per-gene dbNSFP rows by hg38 coord/ref/alt -> HGVSc + transcript + hg19."""
    p = PER_GENE_DIR / f"{gene}.tsv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
    out = {}
    for _, r in df.iterrows():
        try:
            k = (_s(r["#chr"]), int(float(r["pos(1-based)"])), _s(r["ref"]), _s(r["alt"]))
        except (TypeError, ValueError, KeyError):
            continue
        if k in out:
            continue
        aaref = _s(r.get("aaref"))
        aaalt = _s(r.get("aaalt"))
        out[k] = {
            "transcript": _s(r.get("Ensembl_transcriptid")).split(";")[0],
            "hgvs_c_snpeff": _s(r.get("HGVSc_snpEff")).split(";")[0],
            "hgvs_p_snpeff": _s(r.get("HGVSp_snpEff")).split(";")[0],
            "chr_hg19": _s(r.get("hg19_chr")),
            "pos_hg19": _s(r.get("hg19_pos(1-based)")),
            "consequence": "missense" if (aaref and aaalt and aaref != aaalt) else "other",
        }
    return out


def load_syn_catalog() -> dict:
    """index synonymous catalog by (hg19 chr, hg19 pos, ref, alt) -> HGVSc + NM_."""
    if not SYN_CATALOG.exists():
        return {}
    df = pd.read_csv(SYN_CATALOG, sep="\t", dtype=str, na_values=["."]).fillna("")
    out = {}
    chr_col = "chr_grch37" if "chr_grch37" in df.columns else "chr_hg19"
    pos_col = "pos_grch37" if "pos_grch37" in df.columns else "pos_hg19"
    for _, r in df.iterrows():
        try:
            k = (str(r[chr_col]), int(float(r[pos_col])), r["ref"], r["alt"])
        except (TypeError, ValueError, KeyError):
            continue
        if k in out:
            continue
        out[k] = {
            "transcript": r.get("transcript", ""),
            "hgvs_c": r.get("hgvs_c", ""),
            "gene": r.get("gene", ""),
        }
    return out


def annotate(df_common: pd.DataFrame, gene: str, dbnsfp: dict,
             syn_cat: dict) -> pd.DataFrame:
    """Attach hg19 + HGVSc + transcript from dbNSFP/syn-catalog where possible."""
    out_rows = []
    for _, r in df_common.iterrows():
        k_hg38 = (str(r["chr_hg38"]), int(r["pos_hg38"]), r["ref"], r["alt"])
        dbn = dbnsfp.get(k_hg38, {})
        chr_hg19 = dbn.get("chr_hg19", "")
        pos_hg19 = dbn.get("pos_hg19", "")
        transcript = ""
        hgvs_c = ""
        consequence = ""
        source_anno = ""
        if dbn:
            transcript = dbn.get("transcript", "")
            hgvs_c = dbn.get("hgvs_c_snpeff", "")
            consequence = dbn.get("consequence", "")
            source_anno = "dbnsfp"
        if chr_hg19 and pos_hg19:
            try:
                k_hg19 = (str(chr_hg19), int(float(pos_hg19)), r["ref"], r["alt"])
                syn = syn_cat.get(k_hg19)
                if syn and not hgvs_c:
                    transcript = syn["transcript"]
                    hgvs_c = syn["hgvs_c"]
                    consequence = "synonymous"
                    source_anno = "syn_catalog"
            except (TypeError, ValueError):
                pass
        out_rows.append({
            "gene": gene,
            "transcript": transcript,
            "hgvs_c": hgvs_c,
            "classification": f"Benign FAF >5% (gnomAD v4.1, {consequence or 'unannotated'})",
            "chr_grch37": chr_hg19,
            "pos_grch37": pos_hg19,
            "chr_hg38": r["chr_hg38"],
            "pos_hg38": r["pos_hg38"],
            "ref": r["ref"],
            "alt": r["alt"],
            "FAF95_grpmax": r["gnomad_v41_faf95_grpmax"],
            "gnomad_v41_af_joint": r["gnomad_v41_af_joint"],
            "gnomad_v41_grpmax_anc": r["gnomad_v41_grpmax_anc"],
            "gnomad_v41_nhomalt": r["gnomad_v41_nhomalt"],
            "gnomad_v41_filter": r["gnomad_v41_filter"],
            "consequence": consequence or "unannotated",
            "anno_source": source_anno or "none",
            "vv_status": "ok" if hgvs_c else "needs_hgvs",
            "source": "gnomad_common",
        })
    return pd.DataFrame(out_rows)


def process_gene(gene: str, info: dict, threshold: float, syn_cat: dict,
                 buffer: int) -> tuple[str, pd.DataFrame]:
    chrom = info["chr"]
    start = max(1, info["start"] - buffer)
    end = info["end"] + buffer
    df = fetch_gene_common(gene, chrom, start, end, threshold)
    if df.empty:
        return gene, pd.DataFrame()
    dbnsfp = load_dbnsfp_lookup(gene)
    return gene, annotate(df, gene, dbnsfp, syn_cat)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                    help=f"FAF95 grpmax threshold (default {DEFAULT_THRESHOLD})")
    ap.add_argument("--buffer", type=int, default=50,
                    help="bp to extend each gene range (default 50)")
    ap.add_argument("--genes", nargs="*", help="restrict to a subset of genes")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    if not RANGES_JSON.exists():
        sys.exit(f"missing {RANGES_JSON} -- run gene-range builder first")

    ranges = json.loads(RANGES_JSON.read_text())
    if args.genes:
        ranges = {g: ranges[g] for g in args.genes if g in ranges}
    print(f"Pulling gnomAD v4.1 common variants (FAF95 >= {args.threshold}) "
          f"for {len(ranges)} genes ...")

    OUT_PER_GENE.mkdir(parents=True, exist_ok=True)
    OUT_COMBINED.parent.mkdir(parents=True, exist_ok=True)

    syn_cat = load_syn_catalog()
    print(f"  synonymous catalog index: {len(syn_cat):,} entries")

    all_dfs = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_gene, g, info, args.threshold,
                             syn_cat, args.buffer): g
                   for g, info in ranges.items()}
        for i, fut in enumerate(as_completed(futures), 1):
            gene, df = fut.result()
            if df is not None and not df.empty:
                (OUT_PER_GENE / f"{gene}.tsv").write_text(df.to_csv(sep="\t", index=False))
                all_dfs.append(df)
            print(f"  [{i:>2}/{len(ranges)}] {gene:<10}  {len(df):>5} common variants"
                  if df is not None else f"  [{i:>2}/{len(ranges)}] {gene:<10}  (skip)")

    if not all_dfs:
        print("No common variants found.")
        return 0
    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(OUT_COMBINED, sep="\t", index=False)
    dt = time.time() - t0

    print(f"\n=== Summary  ({dt:.0f}s) ===")
    print(f"Total common variants: {len(combined):,}")
    print(f"  by consequence:")
    print("  " + combined["consequence"].value_counts().to_string().replace("\n", "\n  "))
    print(f"  by anno_source:")
    print("  " + combined["anno_source"].value_counts().to_string().replace("\n", "\n  "))
    print(f"  with HGVSc resolved: {(combined['hgvs_c'] != '').sum():,}  "
          f"({100*(combined['hgvs_c'] != '').sum()/len(combined):.1f}%)")
    print(f"  needs VV HGVSc:      {(combined['vv_status'] == 'needs_hgvs').sum():,}")
    print(f"\nPer-gene files -> {OUT_PER_GENE}/")
    print(f"Combined       -> {OUT_COMBINED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

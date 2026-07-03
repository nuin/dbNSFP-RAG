#!/usr/bin/env python3
"""Backfill BOTH masked (-M 1) and unmasked (-M 0) SpliceAI DS_MAX into the
v3 upload (and refresh per-gene files).

Adds a new column `SpliceAI_unmasked` next to existing `SpliceAI_masked`,
so analysts can compare against spliceai.org (which defaults to unmasked).
Optionally re-validates W2/W3 using the unmasked value (which is what the
project's "<=0.1" criterion was originally specified against).

Usage:
  uv run python scripts/backfill_spliceai_dual.py                  # add column, don't re-filter
  uv run python scripts/backfill_spliceai_dual.py --apply-w2-w3    # also re-drop rows failing unmasked threshold
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

UP = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
PER_GENE_OUT = Path("cp_new_bundle/outputs/new_genes_68_v3/per_gene")
PER_GENE_DBNSFP = Path("data/exports/cp_new")
VCF_DIR = Path("data/exports/cp_new")

MASKED_PATTERNS = ("cp_new.spliceai.vcf", "cp_new.bed.chunk_*.spliceai.vcf",
                   "cp_new_v2_needed.spliceai.vcf")
UNMASKED_PATTERNS = ("cp_new_v2_needed.spliceai_unmasked.vcf",)

W2_SPLICEAI_MAX = 0.1
W3_SPLICEAI_MAX = 0.1


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


def parse_spliceai_vcf(path: Path) -> dict[tuple[str, int, str, str], float]:
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    out = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8: continue
            chrom, pos, _id, ref, alt, _q, _f, info = parts[:8]
            chrom = chrom.removeprefix("chr")
            try: pos_i = int(pos)
            except ValueError: continue
            spliceai = None
            for kv in info.split(";"):
                if kv.startswith("SpliceAI="):
                    spliceai = kv[len("SpliceAI="):]; break
            if not spliceai: continue
            best = 0.0
            for entry in spliceai.split(","):
                fields = entry.split("|")
                if len(fields) < 6: continue
                try:
                    ds = max(float(x) for x in fields[2:6] if x not in ("", "."))
                except ValueError: continue
                if ds > best: best = ds
            key = (chrom, pos_i, ref, alt)
            if key not in out or best > out[key]:
                out[key] = best
    return out


def load_vcfs(patterns) -> dict:
    files = []
    for pat in patterns:
        files.extend(VCF_DIR.glob(pat))
    m = {}
    for vcf in sorted(set(files)):
        before = len(m)
        d = parse_spliceai_vcf(vcf)
        m.update(d)
        print(f"  {vcf.name}: {len(d):,} variants  ({len(m)-before:,} new)")
    return m


def build_hg19_to_hg38(genes: list[str]) -> dict[tuple[str, int], int]:
    out = {}
    for g in genes:
        p = PER_GENE_DBNSFP / f"{g}.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, na_values=["."], low_memory=False).fillna("")
        if "hg19_chr" not in df.columns: continue
        for _, r in df.iterrows():
            try:
                k = (_s(r["hg19_chr"]), int(float(r["hg19_pos(1-based)"])))
                v = int(float(r["pos(1-based)"]))
                out.setdefault(k, v)
            except (TypeError, ValueError, KeyError): continue
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply-w2-w3", action="store_true",
                    help="also drop rows where unmasked SpliceAI > 0.1")
    args = ap.parse_args()

    if not UP.exists(): sys.exit(f"missing {UP}")
    df = pd.read_csv(UP, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df):,} rows from v3 upload")

    print(f"\nLoading masked (-M 1) VCFs:")
    masked = load_vcfs(MASKED_PATTERNS)
    print(f"  total masked entries: {len(masked):,}")

    print(f"\nLoading unmasked (-M 0) VCFs:")
    unmasked = load_vcfs(UNMASKED_PATTERNS)
    print(f"  total unmasked entries: {len(unmasked):,}")

    genes = sorted(df["gene"].unique())
    print(f"\nBuilding hg19->hg38 map from {len(genes)} dbNSFP TSVs ...")
    hg = build_hg19_to_hg38(genes)
    print(f"  {len(hg):,} positions mapped")

    # Add SpliceAI_unmasked column if not present
    if "SpliceAI_unmasked" not in df.columns:
        # insert next to SpliceAI_masked
        cols = list(df.columns)
        ix = cols.index("SpliceAI_masked") + 1
        cols.insert(ix, "SpliceAI_unmasked")
        df["SpliceAI_unmasked"] = ""
        df = df[cols]

    n_masked_set = 0
    n_unmasked_set = 0
    for i, r in df.iterrows():
        if not r["chr_grch37"] or not r["pos_grch37"] or not r["ref"] or not r["alt"]:
            continue
        try:
            hg19_pos = int(float(r["pos_grch37"]))
        except (TypeError, ValueError): continue
        hg38_pos = hg.get((r["chr_grch37"], hg19_pos))
        if not hg38_pos: continue
        k = (r["chr_grch37"], hg38_pos, r["ref"], r["alt"])

        # masked: only fill if currently empty / "0 (assumed)"
        m = masked.get(k)
        if m is not None and r["SpliceAI_masked"] in ("", "0 (assumed)"):
            df.at[i, "SpliceAI_masked"] = f"{m:.4f}"
            n_masked_set += 1
        u = unmasked.get(k)
        if u is not None and not r["SpliceAI_unmasked"]:
            df.at[i, "SpliceAI_unmasked"] = f"{u:.4f}"
            n_unmasked_set += 1

    print(f"\nBackfilled:")
    print(f"  SpliceAI_masked   filled: {n_masked_set:,}")
    print(f"  SpliceAI_unmasked filled: {n_unmasked_set:,}")

    if args.apply_w2_w3:
        n_before = len(df)
        keep_idx = []
        for i, r in df.iterrows():
            c = r["classification"]
            sa = _f(r["SpliceAI_unmasked"]) if r["SpliceAI_unmasked"] else _f(r["SpliceAI_masked"])
            # Only drop if we HAVE a measured value and it exceeds threshold
            if sa is None:
                keep_idx.append(i); continue
            if "synonymous, SpliceAI" in c and sa > W2_SPLICEAI_MAX:
                continue
            if c.startswith("Likely_benign FAF >0.1%") and sa > W3_SPLICEAI_MAX:
                continue
            keep_idx.append(i)
        df = df.loc[keep_idx]
        print(f"  --apply-w2-w3: dropped {n_before - len(df):,} rows where unmasked > 0.1")

    df.to_csv(UP, sep="\t", index=False)
    print(f"\nUpdated -> {UP}")

    for g in genes:
        sub = df[df["gene"] == g]
        (PER_GENE_OUT / f"{g}__v3.tsv").write_text(sub.to_csv(sep="\t", index=False))
    print(f"Refreshed {len(genes)} per-gene files")

    # Coverage summary
    print(f"\n=== Final SpliceAI coverage ===")
    for col in ("SpliceAI_masked", "SpliceAI_unmasked"):
        good = (~df[col].isin(["", ".", "0 (assumed)"])).sum()
        print(f"  {col:<20} {good:>6,} / {len(df):,}  ({100*good/len(df):.1f}%)")

    # Comparison: distribution of (unmasked - masked) where both present
    both = df[(df["SpliceAI_masked"].apply(lambda s: _f(s) is not None)) &
              (df["SpliceAI_unmasked"] != "")].copy()
    if len(both):
        both["_m"] = both["SpliceAI_masked"].apply(_f)
        both["_u"] = both["SpliceAI_unmasked"].apply(_f)
        diff = (both["_u"] - both["_m"])
        print(f"\n=== Masked vs Unmasked (both present, n={len(both):,}) ===")
        print(f"  median  unmasked - masked:  {diff.median():+.4f}")
        print(f"  mean    unmasked - masked:  {diff.mean():+.4f}")
        print(f"  rows where unmasked > 0.1 but masked = 0: "
              f"{((both['_u']>0.1) & (both['_m']==0)).sum():,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

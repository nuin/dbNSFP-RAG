#!/usr/bin/env python3
"""Join SpliceAI masked DS_MAX back into cp_new per-gene TSVs.

Inputs:
  --spliceai-vcf : output of `spliceai -M 1` over the consolidated cp_new VCF
                   (i.e. INFO has SpliceAI=ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|...)
  --tsv-dir      : directory of per-gene TSVs from cp_new.py
                   (each row must have #chr, pos(1-based), ref, alt)

For every variant row in every TSV, sets spliceai_ds_max_masked to the maximum
of the four DS scores across all transcripts SpliceAI annotated. Variants with
no SpliceAI annotation get NA.

Writes in-place: each *.tsv is replaced with the annotated version.

Usage:
  spliceai -I cp_new.vcf -O cp_new.spliceai.vcf -R hg38.fa -A grch38 -M 1
  uv run python scripts/annotate_spliceai.py \
      --spliceai-vcf cp_new.spliceai.vcf \
      --tsv-dir data/exports/cp_new
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def parse_spliceai_vcf(vcf_path: Path) -> dict[tuple[str, int, str, str], float]:
    """Read SpliceAI VCF, return {(chr_no_prefix, pos, ref, alt): ds_max}."""
    opener = gzip.open if vcf_path.suffix in (".gz", ".bgz") else open
    out: dict[tuple[str, int, str, str], float] = {}
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _q, _f, info = parts[:8]
            chrom = chrom.removeprefix("chr")
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            # Find SpliceAI= field in INFO
            spliceai = None
            for kv in info.split(";"):
                if kv.startswith("SpliceAI="):
                    spliceai = kv.split("=", 1)[1]
                    break
            if not spliceai:
                continue
            best = None
            for entry in spliceai.split(","):
                fields = entry.split("|")
                if len(fields) < 6:
                    continue
                try:
                    vals = [float(x) for x in fields[2:6] if x not in ("", ".")]
                except ValueError:
                    continue
                if not vals:
                    continue
                m = max(vals)
                best = m if best is None else max(best, m)
            if best is not None:
                key = (chrom, pos_i, ref, alt)
                prev = out.get(key)
                out[key] = best if prev is None else max(prev, best)
    return out


def annotate_tsv(tsv_path: Path, scores: dict[tuple[str, int, str, str], float]) -> int:
    """Read TSV, add/overwrite spliceai_ds_max_masked column, write in place."""
    first = tsv_path.read_text().splitlines()[:1]
    if not first or not first[0].startswith("#chr"):
        return 0
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return 0
    if df.empty or "#chr" not in df.columns:
        return 0

    matched = 0
    new_col = []
    for _, r in df.iterrows():
        try:
            pos = int(float(r["pos(1-based)"]))
        except (TypeError, ValueError):
            new_col.append("")
            continue
        key = (str(r["#chr"]), pos, str(r["ref"]), str(r["alt"]))
        v = scores.get(key)
        if v is None:
            new_col.append("")
        else:
            new_col.append(f"{v:.4f}")
            matched += 1
    df["spliceai_ds_max_masked"] = new_col
    df.to_csv(tsv_path, sep="\t", index=False)
    return matched


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spliceai-vcf", type=Path, required=True,
                   help="SpliceAI output VCF (run with -M 1)")
    p.add_argument("--tsv-dir", type=Path, required=True,
                   help="Per-gene TSV directory (e.g. data/exports/cp_new)")
    args = p.parse_args()

    if not args.spliceai_vcf.exists():
        sys.exit(f"SpliceAI VCF not found: {args.spliceai_vcf}")
    if not args.tsv_dir.is_dir():
        sys.exit(f"TSV directory not found: {args.tsv_dir}")

    print(f"Loading SpliceAI scores from {args.spliceai_vcf}...")
    scores = parse_spliceai_vcf(args.spliceai_vcf)
    print(f"  {len(scores):,} variants with SpliceAI annotation")

    tsvs = sorted(args.tsv_dir.glob("*.tsv"))
    print(f"Annotating {len(tsvs)} per-gene TSVs in {args.tsv_dir}...")

    total_matched = 0
    for tsv in tqdm(tsvs, desc="annotate"):
        total_matched += annotate_tsv(tsv, scores)

    print(f"\nMatched {total_matched:,} TSV rows to SpliceAI scores")
    return 0


if __name__ == "__main__":
    sys.exit(main())

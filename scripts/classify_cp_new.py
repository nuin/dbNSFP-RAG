#!/usr/bin/env python3
"""Apply the three cp_new benign/LB filtering workflows and emit SeqNext rows.

Inputs:
  * data/exports/cp_new/*.tsv     -- per-gene dbNSFP+gnomAD output from cp_new.py
  * BED file with columns chr,start,end,region,gene,transcript
                                  -- defines the lab's coverage + reporting tx

Workflows (priority order, first match wins):
  1. BENIGN -- gnomad_v41_faf95_grpmax > 0.05
  2. BENIGN -- synonymous AND spliceai_ds_max_masked <= 0.1
              AND phastCons100way_vertebrate < 1.0
              AND (if intronic) phyloP < 0.1
  3. LIKELY_BENIGN -- gnomad_v41_faf95_grpmax > 0.001
              AND REVEL_score < 0.290 AND spliceai_ds_max_masked <= 0.1

Output (SeqNext import format), one TSV per gene under data/exports/cp_new/seqnext/:
  gene  transcript  hgvs_c  classification

Each gene file is named {GENE}_seqnext.tsv. A combined cp_new_seqnext.tsv is
also written at the same directory for convenience.

Variants outside any BED region are silently skipped (the assay doesn't cover
them so they can't be reported anyway). Workflows 2 and 3 require SpliceAI;
when spliceai_ds_max_masked is empty (e.g. before running spliceai -M 1) the
script counts those variants as 'pending' and does NOT emit a classification.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from cp_new_rules import classify, load_bed, lookup_bed_transcript
from tqdm import tqdm

DEFAULT_BED = Path("/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed")
DEFAULT_TSV_DIR = Path("data/exports/cp_new")
DEFAULT_OUT_DIR = DEFAULT_TSV_DIR / "seqnext"
DEFAULT_VV_CACHE = DEFAULT_TSV_DIR / ".vv_cache.json"


@dataclass
class ClassificationRun:
    bed_idx: dict
    gene_default_tx: dict[str, str]
    vv_cache: dict[str, str]
    vv_status_cache: dict[str, str]
    counters: dict[str, int]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--tsv-dir", type=Path, default=DEFAULT_TSV_DIR)
    p.add_argument("--bed", type=Path, default=DEFAULT_BED)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument(
        "--vv-cache",
        type=Path,
        default=DEFAULT_VV_CACHE,
        help="VariantValidator cache for accurate intronic detection. "
        "Built incrementally by resolve_hgvs_vv.py; on the first "
        "pipeline run it's empty so classify uses dbNSFP's c.",
    )
    return p.parse_args()


def load_resolution_cache(cache_path: Path) -> tuple[dict, dict]:
    # Load VV cache for accurate intronic detection + status carry-over
    vv_cache: dict[str, str] = {}  # key -> resolved c. (only "ok" entries)
    vv_status_cache: dict[str, str] = {}  # key -> raw status string (ok/flagged:*)
    if cache_path.exists():
        with open(cache_path) as f:
            raw = json.load(f)
        for k, v in raw.items():
            if not (isinstance(v, list) and len(v) == 2):
                continue
            c_part, status = v
            vv_status_cache[k] = status
            if status == "ok" and c_part:
                vv_cache[k] = "c." + c_part.split("c.", 1)[-1]
        print(f"VV cache: {len(vv_status_cache):,} entries ({len(vv_cache):,} ok)\n")
    else:
        print("VV cache: not found -- using dbNSFP HGVSc for intronic detection\n")

    return vv_cache, vv_status_cache


def initial_counters() -> dict[str, int]:
    return {
        "total_rows": 0,
        "not_in_bed": 0,
        "classified": 0,
        "w1_benign_faf5": 0,
        "w2_benign_synonymous": 0,
        "w3_likely_benign": 0,
        "w2_pending_spliceai": 0,
        "w3_pending_spliceai": 0,
    }


def export_row(
    row, identity: dict, result: tuple[str, str], resolved: str | None, status: str
) -> dict:
    raw = str(row.get("HGVSc_snpEff") or "")
    classification, rule = result
    out_row = {
        **identity,
        "hgvs_c": resolved or (raw.split(";")[0] if raw else ""),
        "classification": f"{classification} {rule}",
    }
    score_columns = {
        "PhastCons100way": "phastCons100way_vertebrate",
        "PhyloP100way": "phyloP100way_vertebrate",
        "REVEL": "REVEL_score",
        "SpliceAI_masked": "spliceai_ds_max_masked",
        "FAF95_grpmax": "gnomad_v41_faf95_grpmax",
        "CADD_phred": "CADD_phred",
        "AlphaMissense_pred": "AlphaMissense_pred",
        "ClinVar_sig": "clinvar_clnsig",
    }
    out_row.update({column: row.get(source) or "" for column, source in score_columns.items()})
    out_row["vv_status"] = status
    # Keep the SeqNext column order stable for existing imports.
    columns = [
        "gene",
        "transcript",
        "hgvs_c",
        "classification",
        "chr_grch37",
        "pos_grch37",
        "ref",
        "alt",
        *score_columns,
        "vv_status",
    ]
    return {column: out_row[column] for column in columns}


def classify_row(row, run: ClassificationRun, seen: set) -> dict | None:
    run.counters["total_rows"] += 1
    # dbNSFP `genename` can be ';'-delimited multi-gene/multi-transcript.
    # Try each part and pick the first one that's actually in the BED
    # (i.e. a known panel gene); fall back to first non-empty otherwise.
    raw = str(row.get("genename") or "").strip()
    parts = [p.strip() for p in raw.split(";") if p.strip()]
    gene = next((p for p in parts if p in run.gene_default_tx), parts[0] if parts else "")
    try:
        hg19_chr = str(row.get("hg19_chr") or "").replace("chr", "")
        hg19_pos = int(float(row.get("hg19_pos(1-based)")))
    except (TypeError, ValueError):
        return None

    transcript = lookup_bed_transcript(run.bed_idx, run.gene_default_tx, hg19_chr, hg19_pos, gene)
    if transcript is None:
        run.counters["not_in_bed"] += 1
        return None

    ref = str(row.get("ref") or "")
    alt = str(row.get("alt") or "")
    key = (hg19_chr, hg19_pos, ref, alt, transcript)
    if key in seen:
        return None
    seen.add(key)

    # Look up VV-resolved c. (anchored to BED's NM_) for intronic detection
    vv_key = f"{hg19_chr}-{hg19_pos}-{ref}-{alt}|{transcript}"
    hgvsc_resolved = run.vv_cache.get(vv_key)

    result = classify(row, hgvsc_resolved=hgvsc_resolved)
    if result is None:
        return None

    _, rule = result
    run.counters["classified"] += 1
    if rule == "FAF >5%":
        run.counters["w1_benign_faf5"] += 1
    elif rule.startswith("synonymous"):
        run.counters["w2_benign_synonymous"] += 1
    elif rule.startswith("FAF >0.1%"):
        run.counters["w3_likely_benign"] += 1

    identity = {
        "gene": gene,
        "transcript": transcript,
        "chr_grch37": hg19_chr,
        "pos_grch37": hg19_pos,
        "ref": ref,
        "alt": alt,
    }
    return export_row(
        row, identity, result, hgvsc_resolved, run.vv_status_cache.get(vv_key, "not_resolved")
    )


def classify_file(tsv: Path, run: ClassificationRun) -> list[dict]:
    first = tsv.read_text().splitlines()[:1]
    if not first or not first[0].startswith("#chr"):
        return []
    df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".", ""], low_memory=False)
    seen: set[tuple[str, int, str, str, str]] = set()
    rows = []
    for _, row in df.iterrows():
        classified = classify_row(row, run, seen)
        if classified is not None:
            rows.append(classified)
    return rows


def main() -> int:
    args = parse_args()
    if not args.bed.exists():
        sys.exit(f"BED not found: {args.bed}")
    if not args.tsv_dir.is_dir():
        sys.exit(f"TSV dir not found: {args.tsv_dir}")

    print(f"BED: {args.bed}")
    print(f"TSVs: {args.tsv_dir}")
    print(f"Out: {args.out_dir}\n")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    bed_idx, gene_default_tx = load_bed(args.bed)
    n_bed_regions = sum(len(v) for v in bed_idx.values())
    print(
        f"BED: {n_bed_regions} regions over {len(bed_idx)} (chr, gene) pairs, "
        f"{len(gene_default_tx)} gene defaults"
    )

    vv_cache, vv_status_cache = load_resolution_cache(args.vv_cache)
    counters = initial_counters()
    run = ClassificationRun(bed_idx, gene_default_tx, vv_cache, vv_status_cache, counters)
    combined_rows = []
    for tsv in tqdm(sorted(args.tsv_dir.glob("*.tsv")), desc="classify"):
        gene_rows = classify_file(tsv, run)
        if gene_rows:
            out_path = args.out_dir / f"{tsv.stem}_seqnext.tsv"
            pd.DataFrame(gene_rows).to_csv(out_path, sep="\t", index=False)
            combined_rows.extend(gene_rows)

    # Combined output
    combined_path = args.out_dir / "cp_new_seqnext.tsv"
    if combined_rows:
        pd.DataFrame(combined_rows).to_csv(combined_path, sep="\t", index=False)

    print("\n=== Summary ===")
    for k, v in counters.items():
        print(f"  {k:30s} {v:>8,}")
    print(f"\nClassified rows -> {combined_path} ({len(combined_rows):,} rows)")
    if counters["w2_pending_spliceai"] or counters["w3_pending_spliceai"]:
        print("\nWorkflows 2 & 3 are partially pending: run spliceai -M 1 on")
        print("  data/exports/cp_new/cp_new.vcf, then annotate_spliceai.py,")
        print("  then re-run this script.")

    return 0


if __name__ == "__main__":
    sys.exit(main())

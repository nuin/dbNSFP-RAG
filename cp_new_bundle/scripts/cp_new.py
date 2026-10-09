#!/usr/bin/env python3
"""Pull dbNSFP + gnomAD v4.1 (+ optional SpliceAI masked) for the cp_new panel.

For each of the 165 cp_new genes, this script:
  1. Filters dbNSFP variant rows to that gene (hg38 primary coords).
  2. Batch-queries gnomAD v4.1 joint sites VCF (remote tabix over HTTPS) for
     the gene's hg38 interval to pick up fafmax_faf95_max ("grpmax FAF").
  3. Optionally tabix-queries Illumina's SpliceAI masked precomputed VCFs
     (local path required; --spliceai-snv-vcf / --spliceai-indel-vcf).
     If not provided, the spliceai_ds_max_masked column is left empty.
  4. Writes one TSV per gene to data/exports/cp_new/{GENE}.tsv. Output rows
     carry BOTH GRCh38 coords (for join keys) AND GRCh37 coords pulled
     directly from dbNSFP's hg19_chr / hg19_pos(1-based) columns -- no
     liftover step needed.

Resumable: skips any gene whose output TSV already exists and is non-empty.
Pass --force to overwrite.

This script does NOT modify the SQLite database. Loading into SQLite happens
in a separate step (see task 5).

Requires: tabix (htslib) on PATH. bcftools optional.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from cp_new_data import (
    GNOMAD_URL_TEMPLATE,
    emit_consolidated_vcf,
    fetch_gnomad_region,
    fetch_spliceai_region,
    find_chr_file,
    load_gene_to_chr,
    read_dbnsfp_chr_for_genes,
    require_tabix,
)

from src.config import DBNSFP_DIR, KEEP_COLUMNS
from src.ingest import filter_columns, get_available_columns, list_chromosome_files
from src.panels import get_panel

PANEL_NAME = "cp_new"
EXTRA_COLS = ["hg19_chr", "hg19_pos(1-based)", "#chr", "pos(1-based)"]


@dataclass(frozen=True)
class AnnotationOptions:
    out_dir: Path
    spliceai_snv: Path | None
    spliceai_indel: Path | None
    skip_gnomad: bool
    skip_spliceai: bool


def annotate_and_write_gene(
    gene: str,
    df: pd.DataFrame,
    options: AnnotationOptions,
) -> int:
    """Annotate gene's pre-filtered rows with gnomAD + SpliceAI, write TSV."""
    out_path = options.out_dir / f"{gene}.tsv"
    if df.empty:
        print(f"  {gene}: no dbNSFP variants")
        out_path.write_text("# no variants\n")
        return 0

    df = df.drop(columns=["_gene_upper"], errors="ignore").copy()
    df["#chr"] = df["#chr"].astype(str)
    df["pos(1-based)"] = pd.to_numeric(df["pos(1-based)"], errors="coerce").astype("Int64")
    if "hg19_pos(1-based)" in df.columns:
        df["hg19_pos(1-based)"] = pd.to_numeric(df["hg19_pos(1-based)"], errors="coerce").astype(
            "Int64"
        )

    chrom = df["#chr"].iloc[0]
    if df["#chr"].nunique() > 1:
        print(f"  {gene}: WARNING multi-contig ({df['#chr'].unique().tolist()})")
    pos_min = int(df["pos(1-based)"].min())
    pos_max = int(df["pos(1-based)"].max())

    if options.skip_gnomad:
        df["gnomad_v41_faf95_grpmax"] = pd.NA
    else:
        gn = fetch_gnomad_region(chrom, pos_min, pos_max)
        if not gn.empty:
            df = df.merge(
                gn,
                left_on=["#chr", "pos(1-based)", "ref", "alt"],
                right_on=["gnomad_chr", "gnomad_pos", "gnomad_ref", "gnomad_alt"],
                how="left",
            ).drop(columns=["gnomad_chr", "gnomad_pos", "gnomad_ref", "gnomad_alt"])
        else:
            for c in (
                "gnomad_v41_faf95_grpmax",
                "gnomad_v41_af_joint",
                "gnomad_v41_popmax_af",
                "gnomad_v41_nhomalt",
                "gnomad_v41_vep_csq",
            ):
                df[c] = pd.NA

    if options.skip_spliceai or (options.spliceai_snv is None and options.spliceai_indel is None):
        df["spliceai_ds_max_masked"] = pd.NA
    else:
        sp = fetch_spliceai_region(
            options.spliceai_snv, options.spliceai_indel, chrom, pos_min, pos_max
        )
        if not sp.empty:
            df = df.merge(
                sp,
                left_on=["#chr", "pos(1-based)", "ref", "alt"],
                right_on=["spliceai_chr", "spliceai_pos", "spliceai_ref", "spliceai_alt"],
                how="left",
            ).drop(columns=["spliceai_chr", "spliceai_pos", "spliceai_ref", "spliceai_alt"])
        else:
            df["spliceai_ds_max_masked"] = pd.NA

    df.to_csv(out_path, sep="\t", index=False)
    print(f"  {gene}: {len(df):,} variants -> {out_path.name}")
    return len(df)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/exports/cp_new"),
        help="Per-gene TSV output directory",
    )
    p.add_argument(
        "--data-dir", type=Path, default=DBNSFP_DIR, help="dbNSFP directory (chr*.gz files)"
    )
    p.add_argument(
        "--spliceai-snv-vcf",
        type=Path,
        default=None,
        help="Illumina spliceai_scores.masked.snv.hg38.vcf.gz path (optional)",
    )
    p.add_argument(
        "--spliceai-indel-vcf",
        type=Path,
        default=None,
        help="Illumina spliceai_scores.masked.indel.hg38.vcf.gz path (optional)",
    )
    p.add_argument(
        "--skip-gnomad",
        action="store_true",
        help="Skip gnomAD network calls (dry pass, dbNSFP only)",
    )
    p.add_argument("--skip-spliceai", action="store_true", help="Leave SpliceAI column empty")
    p.add_argument("--genes", nargs="+", help="Restrict to these gene symbols (debug)")
    p.add_argument("--force", action="store_true", help="Overwrite existing per-gene TSVs")
    p.add_argument(
        "--emit-vcf",
        type=Path,
        default=None,
        help="After per-gene TSVs are written, emit a consolidated VCF "
        "of unique variants here (for running SpliceAI -M 1 externally)",
    )
    return p.parse_args()


def prepare_inputs(args: argparse.Namespace) -> tuple[list[str], list[Path], list[str]]:
    require_tabix()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    genes = sorted(get_panel(PANEL_NAME))
    if args.genes:
        wanted = {g.upper() for g in args.genes}
        genes = [g for g in genes if g.upper() in wanted]
    print(f"Panel: {PANEL_NAME} ({len(genes)} genes)")
    print(f"Output: {args.out_dir}")
    print(f"dbNSFP: {args.data_dir}")
    print(f"gnomAD: {'SKIPPED' if args.skip_gnomad else GNOMAD_URL_TEMPLATE.format(chrom='N')}")
    print(f"SpliceAI SNV : {args.spliceai_snv_vcf or 'NOT PROVIDED'}")
    print(f"SpliceAI INDEL: {args.spliceai_indel_vcf or 'NOT PROVIDED'}")

    if args.spliceai_snv_vcf and not args.spliceai_snv_vcf.exists():
        sys.exit(f"SpliceAI SNV VCF not found: {args.spliceai_snv_vcf}")
    if args.spliceai_indel_vcf and not args.spliceai_indel_vcf.exists():
        sys.exit(f"SpliceAI INDEL VCF not found: {args.spliceai_indel_vcf}")

    chr_files = list_chromosome_files(args.data_dir)
    if not chr_files:
        sys.exit(f"No dbNSFP chr files in {args.data_dir}")
    available = get_available_columns(chr_files[0])
    usecols = filter_columns(available, KEEP_COLUMNS + EXTRA_COLS)

    return genes, chr_files, usecols


def pending_genes(args: argparse.Namespace, genes: list[str]) -> dict[str, list[str]]:
    # Build gene -> chromosome map ONCE
    print("\nLoading gene->chr map from dbNSFP gene file...")
    gene_to_chr = load_gene_to_chr(args.data_dir)
    print(f"Loaded chromosome assignments for {len(gene_to_chr):,} genes")

    # Group requested genes by chromosome (skip unknowns)
    by_chrom: dict[str, list[str]] = {}
    unknown = []
    for g in genes:
        chrom = gene_to_chr.get(g.upper())
        if chrom is None:
            unknown.append(g)
            continue
        by_chrom.setdefault(chrom, []).append(g)
    if unknown:
        print(
            f"WARNING: no chromosome assignment for {len(unknown)} gene(s): "
            f"{', '.join(unknown[:10])}{'...' if len(unknown) > 10 else ''}"
        )

    # Filter to genes that still need work (resume support)
    todo_by_chrom: dict[str, list[str]] = {}
    for chrom, gs in by_chrom.items():
        pending = [
            g
            for g in gs
            if args.force
            or not (args.out_dir / f"{g}.tsv").exists()
            or (args.out_dir / f"{g}.tsv").stat().st_size == 0
        ]
        if pending:
            todo_by_chrom[chrom] = pending

    total_pending = sum(len(v) for v in todo_by_chrom.values())
    print(f"Pending: {total_pending} gene(s) across {len(todo_by_chrom)} chromosome(s)")

    return todo_by_chrom


def process_chromosomes(args, chr_files, usecols, todo_by_chrom) -> int:
    options = AnnotationOptions(
        out_dir=args.out_dir,
        spliceai_snv=args.spliceai_snv_vcf,
        spliceai_indel=args.spliceai_indel_vcf,
        skip_gnomad=args.skip_gnomad,
        skip_spliceai=args.skip_spliceai,
    )
    total = 0
    for chrom in sorted(todo_by_chrom.keys(), key=lambda x: (x.isdigit() is False, x)):
        chr_file = find_chr_file(chr_files, chrom)
        if chr_file is None:
            print(
                f"chr{chrom}: no dbNSFP chr file found, skipping {len(todo_by_chrom[chrom])} gene(s)"
            )
            continue
        pending = todo_by_chrom[chrom]
        print(f"\n=== chr{chrom}: {len(pending)} gene(s) -> {chr_file.name} ===")
        genes_upper = {g.upper() for g in pending}
        chr_df = read_dbnsfp_chr_for_genes(chr_file, genes_upper, usecols)
        if chr_df.empty:
            print(f"  no rows matched on chr{chrom}")
            for g in pending:
                (args.out_dir / f"{g}.tsv").write_text("# no variants\n")
            continue
        for g in tqdm(pending, desc=f"chr{chrom}", leave=False):
            gdf = chr_df[chr_df["_gene_upper"] == g.upper()]
            total += annotate_and_write_gene(
                gene=g,
                df=gdf,
                options=options,
            )

    return total


def main() -> int:
    args = parse_args()
    genes, chr_files, usecols = prepare_inputs(args)
    todo_by_chrom = pending_genes(args, genes)
    total = process_chromosomes(args, chr_files, usecols, todo_by_chrom)
    print(f"\nTotal variants written: {total:,}")
    if args.emit_vcf:
        emit_consolidated_vcf(args.out_dir, args.emit_vcf)
    return 0


if __name__ == "__main__":
    sys.exit(main())

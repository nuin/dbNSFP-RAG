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
import gzip
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import DBNSFP_DIR, KEEP_COLUMNS
from src.ingest import list_chromosome_files, get_available_columns, filter_columns
from src.panels import get_panel

PANEL_NAME = "cp_new"
GNOMAD_URL_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"
)

# dbNSFP columns we always want, even if not in KEEP_COLUMNS
EXTRA_COLS = ["hg19_chr", "hg19_pos(1-based)", "#chr", "pos(1-based)"]

# Gene-level file in dbNSFP carries the chromosome assignment for every gene
DBNSFP_GENE_FILE = "dbNSFP5.3_gene.gz"


@dataclass
class GnomadHit:
    faf95_grpmax: float | None
    af_joint: float | None
    popmax_af: float | None
    nhomalt: int | None


# ---------------------------------------------------------------------------
# Tabix helpers
# ---------------------------------------------------------------------------

def require_tabix() -> str:
    path = shutil.which("tabix")
    if not path:
        sys.exit("ERROR: tabix not found on PATH. brew install htslib  (or apt install tabix)")
    return path


def tabix_region(url_or_path: str, chrom: str, start: int, end: int) -> list[str]:
    """Run tabix on a region and return raw VCF lines."""
    region = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ["tabix", url_or_path, region],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
        )
    except subprocess.TimeoutExpired:
        print(f"  tabix TIMEOUT on {region} -> {url_or_path}", file=sys.stderr)
        return []
    if result.returncode != 0:
        print(f"  tabix error ({result.returncode}) {region}: {result.stderr.strip()[:200]}",
              file=sys.stderr)
        return []
    return [line for line in result.stdout.splitlines() if line and not line.startswith("#")]


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# gnomAD slice -> DataFrame
# ---------------------------------------------------------------------------

def fetch_gnomad_region(chrom: str, start: int, end: int) -> pd.DataFrame:
    """Pull a region from gnomAD v4.1 joint sites and return key fields."""
    url = GNOMAD_URL_TEMPLATE.format(chrom=chrom)
    lines = tabix_region(url, f"chr{chrom}", start, end)
    rows = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        c, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
        # Multiallelic VCF rows are pre-split in gnomAD sites; ALT is single.
        ifields = parse_info(info)
        rows.append({
            "gnomad_chr": c.removeprefix("chr"),
            "gnomad_pos": int(pos),
            "gnomad_ref": ref,
            "gnomad_alt": alt,
            # gnomAD v4.1 joint VCF uses _joint suffix; older releases dropped it
            "gnomad_v41_faf95_grpmax": ifields.get("fafmax_faf95_max_joint") or ifields.get("fafmax_faf95_max"),
            "gnomad_v41_af_joint": ifields.get("AF_joint") or ifields.get("AF"),
            "gnomad_v41_grpmax_anc": ifields.get("fafmax_faf95_max_gen_anc_joint") or ifields.get("grpmax_joint"),
            "gnomad_v41_nhomalt": ifields.get("nhomalt_joint") or ifields.get("nhomalt"),
            "gnomad_v41_filter": parts[6],  # PASS / EXOMES_FILTERED / etc.
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# SpliceAI masked slice -> DataFrame
# ---------------------------------------------------------------------------

def fetch_spliceai_region(
    snv_vcf: Path | None,
    indel_vcf: Path | None,
    chrom: str,
    start: int,
    end: int,
) -> pd.DataFrame:
    """Pull a region from Illumina SpliceAI masked VCFs. Returns DS_MAX per site."""
    rows = []
    for src in (snv_vcf, indel_vcf):
        if src is None:
            continue
        lines = tabix_region(str(src), chrom, start, end)
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            c, pos, _id, ref, alt, _q, _f, info = parts[:8]
            # SpliceAI INFO: SpliceAI=ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            spliceai = parse_info(info).get("SpliceAI")
            if not spliceai:
                continue
            best = 0.0
            for entry in spliceai.split(","):
                fields = entry.split("|")
                if len(fields) < 6:
                    continue
                try:
                    ds = max(float(x) for x in fields[2:6] if x not in ("", "."))
                except ValueError:
                    continue
                best = max(best, ds)
            rows.append({
                "spliceai_chr": c.removeprefix("chr"),
                "spliceai_pos": int(pos),
                "spliceai_ref": ref,
                "spliceai_alt": alt,
                "spliceai_ds_max_masked": best,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Gene -> chromosome map (read once from dbNSFP gene file)
# ---------------------------------------------------------------------------

def load_gene_to_chr(data_dir: Path) -> dict[str, str]:
    """Read dbNSFP gene file and return {GENE_UPPER: chr_str}."""
    gf = data_dir / DBNSFP_GENE_FILE
    if not gf.exists():
        sys.exit(f"Gene file not found: {gf}")
    df = pd.read_csv(gf, sep="\t", usecols=["Gene_name", "chr"], dtype=str, na_values=".")
    df = df.dropna(subset=["Gene_name", "chr"])
    return dict(zip(df["Gene_name"].str.upper(), df["chr"].astype(str)))


def find_chr_file(chr_files: list[Path], chrom: str) -> Path | None:
    """Pick the dbNSFP chr file for a given chromosome label."""
    needle = f"chr{chrom}."
    for f in chr_files:
        if needle in f.name:
            return f
    return None


# ---------------------------------------------------------------------------
# dbNSFP slice for a set of genes on one chromosome
# ---------------------------------------------------------------------------

def read_dbnsfp_chr_for_genes(
    chr_file: Path,
    genes_upper: set[str],
    usecols: list[str],
) -> pd.DataFrame:
    """Single pass through one chr file, keep rows for any gene in genes_upper."""
    keep = []
    with gzip.open(chr_file, "rt") as f:
        reader = pd.read_csv(
            f,
            sep="\t",
            usecols=usecols,
            chunksize=50000,
            na_values=".",
            low_memory=False,
            dtype={
                "#chr": str, "ref": str, "alt": str,
                "aaref": str, "aaalt": str, "genename": str,
                "hg19_chr": str,
            },
        )
        for chunk in reader:
            up = chunk["genename"].fillna("").str.upper()
            mask = up.isin(genes_upper)
            if mask.any():
                hit = chunk[mask].copy()
                hit["_gene_upper"] = up[mask].values
                keep.append(hit)
    if not keep:
        return pd.DataFrame(columns=usecols + ["_gene_upper"])
    return pd.concat(keep, ignore_index=True)


# ---------------------------------------------------------------------------
# Per-gene driver (operates on already-loaded chr DataFrame)
# ---------------------------------------------------------------------------

def annotate_and_write_gene(
    gene: str,
    df: pd.DataFrame,
    out_dir: Path,
    spliceai_snv: Path | None,
    spliceai_indel: Path | None,
    skip_gnomad: bool,
    skip_spliceai: bool,
) -> int:
    """Annotate gene's pre-filtered rows with gnomAD + SpliceAI, write TSV."""
    out_path = out_dir / f"{gene}.tsv"
    if df.empty:
        print(f"  {gene}: no dbNSFP variants")
        out_path.write_text("# no variants\n")
        return 0

    df = df.drop(columns=["_gene_upper"], errors="ignore").copy()
    df["#chr"] = df["#chr"].astype(str)
    df["pos(1-based)"] = pd.to_numeric(df["pos(1-based)"], errors="coerce").astype("Int64")
    if "hg19_pos(1-based)" in df.columns:
        df["hg19_pos(1-based)"] = pd.to_numeric(df["hg19_pos(1-based)"], errors="coerce").astype("Int64")

    chrom = df["#chr"].iloc[0]
    if df["#chr"].nunique() > 1:
        print(f"  {gene}: WARNING multi-contig ({df['#chr'].unique().tolist()})")
    pos_min = int(df["pos(1-based)"].min())
    pos_max = int(df["pos(1-based)"].max())

    if skip_gnomad:
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
            for c in ("gnomad_v41_faf95_grpmax", "gnomad_v41_af_joint",
                      "gnomad_v41_popmax_af", "gnomad_v41_nhomalt", "gnomad_v41_vep_csq"):
                df[c] = pd.NA

    if skip_spliceai or (spliceai_snv is None and spliceai_indel is None):
        df["spliceai_ds_max_masked"] = pd.NA
    else:
        sp = fetch_spliceai_region(spliceai_snv, spliceai_indel, chrom, pos_min, pos_max)
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=Path("data/exports/cp_new"),
                   help="Per-gene TSV output directory")
    p.add_argument("--data-dir", type=Path, default=DBNSFP_DIR,
                   help="dbNSFP directory (chr*.gz files)")
    p.add_argument("--spliceai-snv-vcf", type=Path, default=None,
                   help="Illumina spliceai_scores.masked.snv.hg38.vcf.gz path (optional)")
    p.add_argument("--spliceai-indel-vcf", type=Path, default=None,
                   help="Illumina spliceai_scores.masked.indel.hg38.vcf.gz path (optional)")
    p.add_argument("--skip-gnomad", action="store_true",
                   help="Skip gnomAD network calls (dry pass, dbNSFP only)")
    p.add_argument("--skip-spliceai", action="store_true",
                   help="Leave SpliceAI column empty")
    p.add_argument("--genes", nargs="+",
                   help="Restrict to these gene symbols (debug)")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing per-gene TSVs")
    p.add_argument("--emit-vcf", type=Path, default=None,
                   help="After per-gene TSVs are written, emit a consolidated VCF "
                        "of unique variants here (for running SpliceAI -M 1 externally)")
    args = p.parse_args()

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
        print(f"WARNING: no chromosome assignment for {len(unknown)} gene(s): "
              f"{', '.join(unknown[:10])}{'...' if len(unknown) > 10 else ''}")

    # Filter to genes that still need work (resume support)
    todo_by_chrom: dict[str, list[str]] = {}
    for chrom, gs in by_chrom.items():
        pending = [g for g in gs
                   if args.force or not (args.out_dir / f"{g}.tsv").exists()
                   or (args.out_dir / f"{g}.tsv").stat().st_size == 0]
        if pending:
            todo_by_chrom[chrom] = pending

    total_pending = sum(len(v) for v in todo_by_chrom.values())
    print(f"Pending: {total_pending} gene(s) across {len(todo_by_chrom)} chromosome(s)")

    total = 0
    for chrom in sorted(todo_by_chrom.keys(), key=lambda x: (x.isdigit() is False, x)):
        chr_file = find_chr_file(chr_files, chrom)
        if chr_file is None:
            print(f"chr{chrom}: no dbNSFP chr file found, skipping {len(todo_by_chrom[chrom])} gene(s)")
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
                out_dir=args.out_dir,
                spliceai_snv=args.spliceai_snv_vcf,
                spliceai_indel=args.spliceai_indel_vcf,
                skip_gnomad=args.skip_gnomad,
                skip_spliceai=args.skip_spliceai,
            )

    print(f"\nTotal variants written: {total:,}")

    if args.emit_vcf:
        emit_consolidated_vcf(args.out_dir, args.emit_vcf)
    return 0


def emit_consolidated_vcf(tsv_dir: Path, vcf_path: Path) -> None:
    """Write a minimal VCF of unique (chr, pos, ref, alt) across all per-gene TSVs.

    Output is for feeding to `spliceai -M 1` externally. Variants are written
    in chr,pos sort order. Chromosomes get a `chr` prefix to match hg38 builds.
    """
    seen: set[tuple[str, int, str, str]] = set()
    for tsv in sorted(tsv_dir.glob("*.tsv")):
        # Skip placeholder "# no variants" files; they have no header
        first = tsv.read_text().splitlines()[:1]
        if not first or not first[0].startswith("#chr"):
            continue
        try:
            df = pd.read_csv(tsv, sep="\t",
                             usecols=["#chr", "pos(1-based)", "ref", "alt"],
                             dtype=str)
        except (ValueError, pd.errors.EmptyDataError):
            continue
        for _, r in df.iterrows():
            try:
                pos = int(float(r["pos(1-based)"]))
            except (TypeError, ValueError):
                continue
            seen.add((str(r["#chr"]), pos, str(r["ref"]), str(r["alt"])))

    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    with open(vcf_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=.,Number=0,Type=Flag,Description="placeholder">\n')
        # contig lines for chr1..22, X, Y -- SpliceAI doesn't strictly need them
        # but downstream tools (bcftools) complain without them.
        for c in [str(i) for i in range(1, 23)] + ["X", "Y"]:
            out.write(f"##contig=<ID=chr{c}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in sorted(seen, key=lambda x: (x[0], x[1])):
            chrom_out = chrom if chrom.startswith("chr") else f"chr{chrom}"
            out.write(f"{chrom_out}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"Emitted {len(seen):,} unique variants to {vcf_path}")


if __name__ == "__main__":
    sys.exit(main())

"""Data ingestion module for dbNSFP files."""

import gzip
import io
import re
import zipfile
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from .config import KEEP_COLUMNS, RAW_DATA_DIR, CHUNK_SIZE, DBNSFP_ZIP


def filter_by_genes(df: pd.DataFrame, genes: set[str], gene_column: str = "genename") -> pd.DataFrame:
    """
    Filter DataFrame to only include rows matching target genes.

    Args:
        df: Input DataFrame chunk
        genes: Set of gene symbols to keep (case-insensitive)
        gene_column: Column containing gene names

    Returns:
        Filtered DataFrame with only matching genes
    """
    if gene_column not in df.columns:
        return df

    # Normalize gene names for comparison
    genes_upper = {g.upper() for g in genes}
    mask = df[gene_column].fillna("").str.upper().isin(genes_upper)
    return df[mask]


# Chromosome sort order (natural: 1-22, X, Y, M)
CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "M"]

# Chromosome order by file size (smallest to largest based on dbNSFP5.3.1a)
CHR_SIZE_ORDER = ["M", "Y", "21", "18", "13", "22", "20", "14", "X", "15",
                  "8", "9", "10", "4", "16", "5", "7", "6", "12", "17",
                  "11", "3", "19", "2", "1"]


def sort_chr_files(files: list[str], by_size: bool = True) -> list[str]:
    """Sort chromosome files by size (smallest first) or natural order."""
    order = CHR_SIZE_ORDER if by_size else CHR_ORDER

    def chr_key(f):
        match = re.search(r"chr([0-9XYMN]+)", f)
        if match:
            chr_val = match.group(1)
            try:
                return order.index(chr_val)
            except ValueError:
                return 100
        return 100
    return sorted(files, key=chr_key)


# =============================================================================
# ZIP-based reading (Option 1 - no extraction needed)
# =============================================================================

def list_chromosome_files_in_zip(zip_path: Path | None = None) -> list[str]:
    """List chromosome variant files inside the dbNSFP zip."""
    zip_path = zip_path or DBNSFP_ZIP

    if not zip_path or not zip_path.exists():
        raise FileNotFoundError(f"Zip file not found: {zip_path}")

    with zipfile.ZipFile(zip_path, "r") as zf:
        all_files = zf.namelist()
        chr_files = [f for f in all_files if "_variant.chr" in f and f.endswith(".gz")]

    if not chr_files:
        raise FileNotFoundError(f"No chromosome variant files found in {zip_path}")

    return sort_chr_files(chr_files)


def get_columns_from_zip(zip_path: Path, inner_file: str) -> list[str]:
    """Read header from a gzipped file inside a zip."""
    with zipfile.ZipFile(zip_path, "r") as zf:
        with zf.open(inner_file) as f:
            with gzip.open(f, "rt") as gz:
                header = gz.readline().strip().split("\t")
    return header


def parse_chromosome_from_zip(
    zip_path: Path,
    inner_file: str,
    columns: list[str] | None = None,
    chunksize: int = CHUNK_SIZE,
):
    """
    Stream parse a chromosome file directly from inside a zip.

    The file is gzipped inside the zip, so we decompress on-the-fly.
    """
    columns = columns or KEEP_COLUMNS

    # Get available columns
    available = get_columns_from_zip(zip_path, inner_file)
    usecols = filter_columns(available, columns)

    missing = set(columns) - set(usecols)
    if missing:
        fname = Path(inner_file).name
        print(f"  Note: {len(missing)} columns not in {fname}")

    with zipfile.ZipFile(zip_path, "r") as zf:
        with zf.open(inner_file) as f:
            with gzip.open(f, "rt") as gz:
                reader = pd.read_csv(
                    gz,
                    sep="\t",
                    usecols=usecols,
                    chunksize=chunksize,
                    na_values=".",
                    low_memory=False,
                    dtype={
                        "#chr": str,
                        "ref": str,
                        "alt": str,
                        "aaref": str,
                        "aaalt": str,
                        "genename": str,
                    },
                )
                for chunk in reader:
                    if "#chr" in chunk.columns:
                        chunk["#chr"] = chunk["#chr"].astype(str).str.replace("chr", "")
                    yield chunk


def iter_variants_from_zip(
    zip_path: Path | None = None,
    chromosomes: list[str] | None = None,
    columns: list[str] | None = None,
    show_progress: bool = True,
):
    """
    Iterate through chromosome files in zip, yielding variant chunks.

    Args:
        zip_path: Path to dbNSFP zip file
        chromosomes: Optional list of chromosomes to process (e.g., ["17", "13"])
        columns: Columns to extract
        show_progress: Show progress bar
    """
    zip_path = zip_path or DBNSFP_ZIP
    chr_files = list_chromosome_files_in_zip(zip_path)

    # Filter chromosomes if specified
    if chromosomes:
        chr_set = set(str(c).replace("chr", "") for c in chromosomes)
        chr_files = [f for f in chr_files if any(f"chr{c}.gz" in f or f"chr{c}." in f for c in chr_set)]

    for chr_file in chr_files:
        fname = Path(chr_file).name
        desc = f"Processing {fname}"
        chunks = parse_chromosome_from_zip(zip_path, chr_file, columns)

        if show_progress:
            chunks = tqdm(chunks, desc=desc, unit="batch")

        for chunk in chunks:
            yield chunk


# =============================================================================
# Directory-based reading (extracted files)
# =============================================================================

def list_chromosome_files(data_dir: Path | None = None) -> list[Path]:
    """Find all dbNSFP chromosome variant files in a directory."""
    data_dir = data_dir or RAW_DATA_DIR
    files = list(data_dir.glob("dbNSFP*_variant.chr*.gz"))

    if not files:
        # Check for nested directory
        files = list(data_dir.glob("*/dbNSFP*_variant.chr*.gz"))

    if not files:
        raise FileNotFoundError(
            f"No dbNSFP chromosome files found in {data_dir}. "
            "Expected files like: dbNSFP5.0_variant.chr1.gz"
        )

    return [Path(f) for f in sort_chr_files([str(f) for f in files])]


def get_available_columns(chr_file: Path) -> list[str]:
    """Read header to get available columns."""
    with gzip.open(chr_file, "rt") as f:
        header = f.readline().strip().split("\t")
    return header


def filter_columns(available: list[str], wanted: list[str]) -> list[str]:
    """Return intersection of available and wanted columns, preserving order."""
    available_set = set(available)
    return [col for col in wanted if col in available_set]


def parse_chromosome_file(
    chr_file: Path,
    columns: list[str] | None = None,
    chunksize: int = CHUNK_SIZE,
):
    """
    Stream parse a dbNSFP chromosome file from disk.
    """
    columns = columns or KEEP_COLUMNS

    available = get_available_columns(chr_file)
    usecols = filter_columns(available, columns)

    missing = set(columns) - set(usecols)
    if missing:
        print(f"  Note: {len(missing)} columns not in {chr_file.name}")

    with gzip.open(chr_file, "rt") as f:
        reader = pd.read_csv(
            f,
            sep="\t",
            usecols=usecols,
            chunksize=chunksize,
            na_values=".",
            low_memory=False,
            dtype={
                "#chr": str,
                "ref": str,
                "alt": str,
                "aaref": str,
                "aaalt": str,
                "genename": str,
            },
        )
        for chunk in reader:
            if "#chr" in chunk.columns:
                chunk["#chr"] = chunk["#chr"].astype(str).str.replace("chr", "")
            yield chunk


def iter_all_variants(
    data_dir: Path | None = None,
    columns: list[str] | None = None,
    show_progress: bool = True,
):
    """
    Iterate through all chromosome files in directory.
    """
    files = list_chromosome_files(data_dir)

    for chr_file in files:
        desc = f"Processing {chr_file.name}"
        chunks = parse_chromosome_file(chr_file, columns)

        if show_progress:
            chunks = tqdm(chunks, desc=desc, unit="batch")

        for chunk in chunks:
            yield chunk


# =============================================================================
# Unified interface
# =============================================================================

def iter_variants(
    source: Path | None = None,
    chromosomes: list[str] | None = None,
    columns: list[str] | None = None,
    show_progress: bool = True,
):
    """
    Unified variant iterator - works with both zip files and directories.

    Args:
        source: Path to zip file or directory (auto-detects)
        chromosomes: Optional list of chromosomes to process
        columns: Columns to extract
        show_progress: Show progress bar
    """
    # Auto-detect source
    if source is None:
        if DBNSFP_ZIP and DBNSFP_ZIP.exists():
            source = DBNSFP_ZIP
        else:
            source = RAW_DATA_DIR

    source = Path(source)

    if source.suffix == ".zip" or (source.is_file() and zipfile.is_zipfile(source)):
        # Read from zip
        yield from iter_variants_from_zip(
            source,
            chromosomes=chromosomes,
            columns=columns,
            show_progress=show_progress
        )
    else:
        # Read from directory
        if chromosomes:
            print(f"Warning: chromosome filtering not implemented for directory mode")
        yield from iter_all_variants(source, columns, show_progress)


def main():
    """CLI entry point for testing ingestion."""
    import argparse

    parser = argparse.ArgumentParser(description="Test dbNSFP ingestion")
    parser.add_argument("--zip", type=Path, help="Path to dbNSFP zip file")
    parser.add_argument("--dir", type=Path, help="Path to extracted files")
    parser.add_argument("--chr", type=str, help="Chromosome to test (e.g., 17)")
    parser.add_argument("--sample", type=int, default=5, help="Sample rows to display")
    args = parser.parse_args()

    source = args.zip or args.dir or DBNSFP_ZIP or RAW_DATA_DIR
    chromosomes = [args.chr] if args.chr else None

    print(f"Source: {source}")

    try:
        count = 0
        for chunk in iter_variants(source, chromosomes=chromosomes):
            if count == 0:
                print(f"\nColumns loaded ({len(chunk.columns)}):")
                print(chunk.columns.tolist()[:20], "...")
                print(f"\nSample data ({args.sample} rows):")
                print(chunk.head(args.sample))
            count += len(chunk)
            if count >= 10000:
                break

        print(f"\nProcessed {count:,} variants (stopped at 10k for test)")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

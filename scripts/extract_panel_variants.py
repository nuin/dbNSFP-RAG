#!/usr/bin/env python3
"""Extract all variants for a gene panel to a TSV file."""

import argparse
import gzip
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import DBNSFP_DIR, KEEP_COLUMNS
from src.panels import get_panel
from src.ingest import list_chromosome_files, get_available_columns, filter_columns

import pandas as pd
from tqdm import tqdm


def extract_panel_variants(
    panel_name: str,
    output_file: Path,
    data_dir: Path | None = None,
    limit: int | None = None,
):
    """Extract all variants for genes in a panel."""
    data_dir = data_dir or DBNSFP_DIR
    genes = get_panel(panel_name)
    genes_upper = {g.upper() for g in genes}

    print(f"Panel: {panel_name} ({len(genes)} genes)")
    print(f"Data directory: {data_dir}")
    print(f"Output file: {output_file}")

    chr_files = list_chromosome_files(data_dir)
    print(f"Found {len(chr_files)} chromosome files")

    # Get columns from first file
    available = get_available_columns(chr_files[0])
    usecols = filter_columns(available, KEEP_COLUMNS)

    # Check which columns are missing
    missing = set(KEEP_COLUMNS) - set(usecols)
    if missing:
        print(f"\nMissing columns ({len(missing)}):")
        for col in sorted(missing):
            print(f"  - {col}")

    # Show gnomAD columns available
    gnomad_cols = [c for c in available if 'gnomad' in c.lower()]
    print(f"\nAvailable gnomAD columns ({len(gnomad_cols)}):")
    for col in gnomad_cols[:20]:
        print(f"  - {col}")
    if len(gnomad_cols) > 20:
        print(f"  ... and {len(gnomad_cols) - 20} more")

    total_variants = 0
    first_chunk = True

    with open(output_file, 'w') as out_f:
        for chr_file in tqdm(chr_files, desc="Processing chromosomes"):
            with gzip.open(chr_file, 'rt') as f:
                reader = pd.read_csv(
                    f,
                    sep='\t',
                    usecols=usecols,
                    chunksize=10000,
                    na_values='.',
                    low_memory=False,
                    dtype={
                        '#chr': str,
                        'ref': str,
                        'alt': str,
                        'aaref': str,
                        'aaalt': str,
                        'genename': str,
                    },
                )

                for chunk in reader:
                    # Filter by genes
                    if 'genename' in chunk.columns:
                        mask = chunk['genename'].fillna('').str.upper().isin(genes_upper)
                        filtered = chunk[mask]
                    else:
                        filtered = chunk

                    if len(filtered) > 0:
                        # Write header only once
                        filtered.to_csv(
                            out_f,
                            sep='\t',
                            index=False,
                            header=first_chunk,
                            mode='a' if not first_chunk else 'w'
                        )
                        first_chunk = False
                        total_variants += len(filtered)

                    if limit and total_variants >= limit:
                        print(f"\nReached limit of {limit} variants")
                        break

            if limit and total_variants >= limit:
                break

    print(f"\nExtracted {total_variants:,} variants to {output_file}")
    return total_variants


def main():
    parser = argparse.ArgumentParser(description="Extract panel variants to TSV")
    parser.add_argument("--panel", default="NGSgenes", help="Panel name")
    parser.add_argument("--output", "-o", type=Path, help="Output TSV file")
    parser.add_argument("--data-dir", type=Path, help="dbNSFP data directory")
    parser.add_argument("--limit", type=int, help="Limit number of variants")
    args = parser.parse_args()

    output = args.output or Path(f"data/exports/{args.panel}_variants_raw.tsv")
    output.parent.mkdir(parents=True, exist_ok=True)

    extract_panel_variants(
        panel_name=args.panel,
        output_file=output,
        data_dir=args.data_dir,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()

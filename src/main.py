"""Command-line entry point for the dbNSFP ingestion pipeline."""

import argparse
import logging
import sys
from pathlib import Path

from .build_pipeline import load_checkpoint as load_checkpoint, save_checkpoint as save_checkpoint
from .builders import build_index, build_panel_index, build_sqlite_index
from .config import DBNSFP_DIR, DBNSFP_ZIP
from .ingest import list_chromosome_files_in_zip
from .panels import list_panels


def show_stats():
    """Show index statistics."""
    from .vectorstore import VariantVectorStore

    store = VariantVectorStore()
    count = store.count()

    print("dbNSFP RAG Index Statistics (ClickHouse)")
    print("-" * 40)
    print(f"Total variants indexed: {count:,}")

    if count > 0:
        print("\nSample variants by gene:")
        for gene in ["BRCA1", "BRCA2", "TP53", "EGFR", "KRAS"]:
            results = store.search_by_gene(gene, k=100)
            if results:
                print(f"  {gene}: {len(results)} variants")


def clear_index():
    """Clear the vector index and checkpoint."""
    from .vectorstore import VariantVectorStore

    store = VariantVectorStore()
    count = store.count()
    checkpoint_path = store.db_path / "checkpoint.json"

    if count == 0 and not checkpoint_path.exists():
        print("Index is already empty")
        return 0

    confirm = input(f"Delete {count:,} variants and checkpoint? [y/N]: ")
    if confirm.lower() == "y":
        store.delete_all()
        if checkpoint_path.exists():
            checkpoint_path.unlink()
        print("Index and checkpoint cleared")
    else:
        print("Cancelled")

    return 0


def list_chromosomes(source: Path):
    """List available chromosomes in the source."""
    if source.suffix == ".zip":
        files = list_chromosome_files_in_zip(source)
        print(f"Chromosomes in {source.name}:")
        for f in files:
            print(f"  - {Path(f).name}")
    else:
        print("Listing only supported for zip files")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="dbNSFP RAG Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build full index (auto-resumes if stopped)
  python -m src.main build

  # Start fresh (ignore checkpoint)
  python -m src.main build --fresh

  # Index specific chromosomes only
  python -m src.main build --chr 17 --chr 13

  # Index sample for testing
  python -m src.main build --limit 50000

  # Show index stats
  python -m src.main stats

  # List chromosomes in zip (shows size order)
  python -m src.main list

  # Clear index and checkpoint
  python -m src.main clear
        """,
    )

    # Use directory if available, otherwise zip
    default_source = DBNSFP_DIR if DBNSFP_DIR and DBNSFP_DIR.exists() else DBNSFP_ZIP

    parser.add_argument(
        "--source",
        type=Path,
        default=default_source,
        help="Path to dbNSFP zip or directory",
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    build_parser = subparsers.add_parser("build", help="Build index")
    build_parser.add_argument(
        "--chr",
        action="append",
        dest="chromosomes",
        help="Chromosome(s) to index (can repeat, e.g., --chr 17 --chr 13)",
    )
    build_parser.add_argument(
        "--limit",
        type=int,
        help="Limit number of variants to process",
    )
    build_parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start fresh, ignoring any checkpoint",
    )

    # Stats command
    subparsers.add_parser("stats", help="Show index statistics")

    # List command
    subparsers.add_parser("list", help="List chromosomes in source")

    # Clear command
    subparsers.add_parser("clear", help="Clear the index")

    panel_parser = subparsers.add_parser("build-panel", help="Build panel-specific index")
    panel_parser.add_argument(
        "--panel",
        "-p",
        required=True,
        help="Panel name (e.g., NGSgenes, hereditary_cancer)",
    )
    panel_parser.add_argument(
        "--build",
        "-b",
        choices=["grch37", "grch38"],
        default="grch37",
        help="Genome build (default: grch37)",
    )
    panel_parser.add_argument(
        "--chr",
        action="append",
        dest="chromosomes",
        help="Specific chromosome(s) to process",
    )
    panel_parser.add_argument(
        "--limit",
        type=int,
        help="Limit variants for testing",
    )
    panel_parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start fresh, ignoring checkpoint",
    )

    sqlite_parser = subparsers.add_parser("build-sqlite", help="Build SQLite database (all panels)")
    sqlite_parser.add_argument(
        "--build",
        "-b",
        choices=["grch37", "grch38"],
        default="grch37",
        help="Genome build (default: grch37)",
    )
    sqlite_parser.add_argument(
        "--chr",
        action="append",
        dest="chromosomes",
        help="Specific chromosome(s) to process",
    )
    sqlite_parser.add_argument(
        "--limit",
        type=int,
        help="Limit variants for testing",
    )
    sqlite_parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start fresh, ignoring checkpoint",
    )
    sqlite_parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Override SQLite database path",
    )

    # List-panels command
    subparsers.add_parser("list-panels", help="List available gene panels")

    args = parser.parse_args()

    handlers = {
        "build": lambda: build_index(
            args.source, chromosomes=args.chromosomes, limit=args.limit, resume=not args.fresh
        ),
        "build-panel": lambda: build_panel_index(
            args.source,
            args.panel,
            use_grch37=args.build == "grch37",
            chromosomes=args.chromosomes,
            limit=args.limit,
            resume=not args.fresh,
        ),
        "build-sqlite": lambda: build_sqlite_index(
            args.source,
            use_grch37=args.build == "grch37",
            chromosomes=args.chromosomes,
            limit=args.limit,
            resume=not args.fresh,
            db_path=args.db,
        ),
        "stats": show_stats,
        "list": lambda: list_chromosomes(args.source),
        "clear": clear_index,
        "list-panels": show_panels,
    }
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if args.command is None:
        parser.print_help()
        return 0
    return handlers[args.command]()


def show_panels():
    print("Available gene panels:")
    for name, count in list_panels().items():
        print(f"  {name}: {count} genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())

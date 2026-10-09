"""Database browsing endpoints."""

import logging
import sqlite3

from fastapi import APIRouter, HTTPException, Query

from api.dependencies import get_database

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/variant/{variant_id}")
async def get_variant_info(
    variant_id: str,
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """Get raw variant information from database."""
    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(status_code=404, detail="Variant not found")

    return variant_data


@router.get("/gene/{gene_symbol}")
async def get_gene_variants(
    gene_symbol: str,
    limit: int = 100,
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """Get variants for a specific gene."""
    store = get_database(genome_build)
    results = store.search_by_gene(gene_symbol.upper(), k=limit)

    return {
        "gene": gene_symbol.upper(),
        "genome_build": genome_build,
        "count": len(results),
        "variants": [
            {
                "id": r["id"],
                "pos": r["metadata"].get("pos"),
                "ref": r["metadata"].get("ref"),
                "alt": r["metadata"].get("alt"),
                "consequence": r["metadata"].get("consequence"),
                "clinvar_sig": r["metadata"].get("clinvar_sig"),
                "cadd_phred": r["metadata"].get("cadd_phred"),
            }
            for r in results
        ],
    }


@router.get("/chromosomes")
async def list_chromosomes(
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """List all chromosomes with variant counts, naturally sorted (1-22, X, Y)."""
    store = get_database(genome_build)
    if not hasattr(store, "list_chromosomes"):
        raise HTTPException(status_code=501, detail="Browse requires SQLite backend")
    return {
        "genome_build": genome_build,
        "chromosomes": store.list_chromosomes(),
    }


@router.get("/chromosome/{chrom}/genes")
async def list_genes_on_chromosome(
    chrom: str,
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """List genes on a chromosome with variant counts, sorted by count descending."""
    store = get_database(genome_build)
    if not hasattr(store, "list_genes_by_chromosome"):
        raise HTTPException(status_code=501, detail="Browse requires SQLite backend")
    genes = store.list_genes_by_chromosome(chrom)
    return {
        "chromosome": chrom,
        "genome_build": genome_build,
        "gene_count": len(genes),
        "genes": genes,
    }


@router.get("/panels")
async def list_panels():
    """List available gene panels and variant counts."""
    from src.panels import list_panels as _list_panels

    panels_info = _list_panels()

    # Try to get variant counts from database if available
    panel_counts = {}
    try:
        db = get_database("GRCh37")
        if hasattr(db, "count_by_panel"):
            for name in panels_info:
                panel_counts[name] = db.count_by_panel(name)
    except (OSError, ValueError, RuntimeError, ImportError, sqlite3.Error) as exc:
        logger.warning("Panel variant counts unavailable: %s", exc)
        panel_counts = {}

    return {
        "panels": [
            {
                "name": name,
                "gene_count": gene_count,
                "variant_count": panel_counts.get(name),
            }
            for name, gene_count in panels_info.items()
        ]
    }


@router.get("/panel/{panel_name}/genes")
async def get_panel_genes(panel_name: str):
    """Get the list of genes in a specific panel."""
    from src.panels import get_panel

    try:
        genes = get_panel(panel_name)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "panel": panel_name,
        "gene_count": len(genes),
        "genes": sorted(genes),
    }

"""Gene-level annotation loader for dbNSFP.

This module loads gene-level constraint scores (pLI, LOEUF, mis_z, etc.)
from the dbNSFP gene file and provides a lookup function by gene name.

These annotations are separate from variant-level data and must be joined
during database building.
"""

import gzip
import logging
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

from src.config import DBNSFP_GENE_FILE

logger = logging.getLogger(__name__)

# Columns to extract from gene file
GENE_COLUMNS = [
    "Gene_name",
    "gnomAD_pLI",
    "gnomAD_pRec",
    "gnomAD_pNull",
    "gnomAD_lof.oe",
    "gnomAD_mis.oe",
    "gnomAD_LOEUF",
    "gnomAD_MOEUF",
    "ExAC_pLI",
    "RVIS_percentile_ExAC",
]


def _parse_float(val) -> Optional[float]:
    """Parse float from potentially multi-value or missing field."""
    if pd.isna(val) or val == "." or val == "":
        return None
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, str):
        # Take first value if semicolon-separated
        parts = val.split(";")
        for p in parts:
            p = p.strip()
            if p and p != ".":
                try:
                    return float(p)
                except ValueError:
                    continue
    return None


@lru_cache(maxsize=1)
def load_gene_data(gene_file: Optional[Path] = None) -> dict:
    """
    Load gene-level annotations from dbNSFP gene file.

    Returns a dictionary keyed by gene name with constraint scores.
    Results are cached after first load.
    """
    gene_file = gene_file or DBNSFP_GENE_FILE

    if not gene_file.exists():
        logger.warning(f"Gene file not found: {gene_file}")
        return {}

    logger.info(f"Loading gene annotations from {gene_file}")

    gene_dict = {}

    try:
        # Read the gene file
        with gzip.open(gene_file, "rt") as f:
            # Read header
            header = f.readline().strip().split("\t")

            # Find column indices
            col_indices = {}
            for col in GENE_COLUMNS:
                if col in header:
                    col_indices[col] = header.index(col)
                else:
                    logger.warning(f"Column {col} not found in gene file")

            if "Gene_name" not in col_indices:
                logger.error("Gene_name column not found in gene file")
                return {}

            gene_name_idx = col_indices["Gene_name"]

            # Process rows
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) <= gene_name_idx:
                    continue

                gene_name = fields[gene_name_idx]
                if not gene_name or gene_name == ".":
                    continue

                # Extract constraint scores
                gene_dict[gene_name] = {
                    "gnomad_pli": _parse_float(fields[col_indices["gnomAD_pLI"]])
                        if "gnomAD_pLI" in col_indices else None,
                    "gnomad_prec": _parse_float(fields[col_indices["gnomAD_pRec"]])
                        if "gnomAD_pRec" in col_indices else None,
                    "loeuf": _parse_float(fields[col_indices["gnomAD_LOEUF"]])
                        if "gnomAD_LOEUF" in col_indices else None,
                    "gnomad_lof_oe": _parse_float(fields[col_indices["gnomAD_lof.oe"]])
                        if "gnomAD_lof.oe" in col_indices else None,
                    "gnomad_mis_oe": _parse_float(fields[col_indices["gnomAD_mis.oe"]])
                        if "gnomAD_mis.oe" in col_indices else None,
                    "moeuf": _parse_float(fields[col_indices["gnomAD_MOEUF"]])
                        if "gnomAD_MOEUF" in col_indices else None,
                }

        logger.info(f"Loaded annotations for {len(gene_dict)} genes")
        return gene_dict

    except Exception as e:
        logger.error(f"Error loading gene file: {e}")
        return {}


def get_gene_constraints(gene_name: str, gene_data: Optional[dict] = None) -> dict:
    """
    Get constraint scores for a specific gene.

    Args:
        gene_name: Gene symbol (e.g., "BRCA1")
        gene_data: Pre-loaded gene data dict, or None to load from file

    Returns:
        Dictionary with constraint scores, or empty dict if gene not found
    """
    if gene_data is None:
        gene_data = load_gene_data()

    return gene_data.get(gene_name, {})


def clear_cache():
    """Clear the cached gene data."""
    load_gene_data.cache_clear()

"""
UCSC RepeatMasker API client for BP3 ACMG criterion evaluation.

BP3: In-frame deletions/insertions in a repetitive region without a known function.

This module queries the UCSC Genome Browser API to check if a genomic position
falls within a repetitive region annotated by RepeatMasker.

## Usage:
    from src.repeat_api import is_in_repeat_region, evaluate_bp3

    # Check if a position is in a repeat region
    in_repeat, repeat_class, repeat_family = is_in_repeat_region("17", 41197801, "grch37")

    # Full BP3 evaluation for an in-frame variant
    criterion_met, evidence = evaluate_bp3("17", 41197801, "inframe_deletion", "grch37")

## API Reference:
- UCSC Genome Browser API: https://genome.ucsc.edu/goldenPath/help/api.html
- RepeatMasker track: https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk

## Caching:
- Results are cached using LRU cache (max 10,000 entries) to minimize API calls
- Cache persists for the lifetime of the process

FOR RESEARCH USE ONLY - Not for clinical diagnostic use.
"""

import logging
from functools import lru_cache
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# UCSC API endpoint
UCSC_API_BASE = "https://api.genome.ucsc.edu/getData/track"

# Timeout for API requests (seconds)
REQUEST_TIMEOUT = 10


@lru_cache(maxsize=10000)
def is_in_repeat_region(
    chrom: str,
    pos: int,
    build: str = "grch37",
) -> tuple[bool, Optional[str], Optional[str]]:
    """
    Check if a genomic position is in a repetitive region using UCSC RepeatMasker.

    Args:
        chrom: Chromosome (e.g., "1", "X", "chr1")
        pos: 1-based genomic position
        build: Genome build ("grch37" or "grch38")

    Returns:
        Tuple of (is_in_repeat, repeat_class, repeat_family)
        - is_in_repeat: True if position overlaps a repeat
        - repeat_class: RepeatMasker class (e.g., "LINE", "SINE", "Simple_repeat")
        - repeat_family: Repeat family (e.g., "L1", "Alu", "(CA)n")
    """
    # Map build to UCSC genome name
    genome = "hg19" if build.lower() in ("grch37", "hg19") else "hg38"

    # Normalize chromosome format for UCSC (needs "chr" prefix)
    chrom_clean = str(chrom).replace("chr", "")
    ucsc_chrom = f"chr{chrom_clean}"

    # Query a small window around the position
    start = max(0, pos - 1)  # 0-based start
    end = pos + 1  # 0-based exclusive end

    params = {
        "genome": genome,
        "track": "rmsk",  # RepeatMasker track
        "chrom": ucsc_chrom,
        "start": start,
        "end": end,
    }

    try:
        response = requests.get(UCSC_API_BASE, params=params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        data = response.json()

        # Check if any repeat annotations were returned
        repeats = data.get("rmsk", [])
        if repeats:
            # Get the first overlapping repeat
            repeat = repeats[0]
            repeat_class = repeat.get("repClass", "Unknown")
            repeat_family = repeat.get("repFamily", "Unknown")
            repeat_name = repeat.get("repName", "Unknown")

            logger.debug(
                f"Position {ucsc_chrom}:{pos} in repeat: {repeat_class}/{repeat_family} ({repeat_name})"
            )
            return True, repeat_class, repeat_family

        return False, None, None

    except requests.exceptions.Timeout:
        logger.warning(f"UCSC API timeout for {ucsc_chrom}:{pos}")
        return False, None, None
    except requests.exceptions.RequestException as e:
        logger.warning(f"UCSC API error for {ucsc_chrom}:{pos}: {e}")
        return False, None, None
    except (KeyError, ValueError) as e:
        logger.warning(f"UCSC API parse error for {ucsc_chrom}:{pos}: {e}")
        return False, None, None


def evaluate_bp3(
    chrom: str,
    pos: int,
    consequence: str,
    build: str = "grch37",
) -> tuple[bool, str]:
    """
    Evaluate BP3 criterion: In-frame indel in repetitive region.

    Args:
        chrom: Chromosome
        pos: Genomic position
        consequence: Variant consequence (e.g., "inframe_deletion")
        build: Genome build

    Returns:
        Tuple of (criterion_met, evidence_string)
    """
    # Check if variant is in-frame
    cons_lower = consequence.lower() if consequence else ""
    is_inframe = "inframe" in cons_lower or "in_frame" in cons_lower

    if not is_inframe:
        return False, f"Not an in-frame variant (consequence: {consequence})"

    # Query repeat status
    in_repeat, repeat_class, repeat_family = is_in_repeat_region(chrom, pos, build)

    if in_repeat:
        return True, f"In-frame indel in {repeat_class}/{repeat_family} repetitive region"
    else:
        return False, "In-frame indel but not in repetitive region"


def get_repeat_info(
    chrom: str,
    pos: int,
    build: str = "grch37",
) -> Optional[dict]:
    """
    Get detailed repeat information for a position.

    Args:
        chrom: Chromosome
        pos: Genomic position
        build: Genome build

    Returns:
        Dictionary with repeat details or None if not in repeat
    """
    in_repeat, repeat_class, repeat_family = is_in_repeat_region(chrom, pos, build)

    if in_repeat:
        return {
            "in_repeat": True,
            "repeat_class": repeat_class,
            "repeat_family": repeat_family,
        }
    return None

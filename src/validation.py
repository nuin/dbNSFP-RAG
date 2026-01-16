"""Input validation module for ACMG Variant Classification System.

This module provides validation functions for variant input to prevent
silent failures and incorrect lookups.

Risk Mitigations:
- FM-01: Invalid variant accepted
- FM-02: Malformed chromosome
- FM-03: Negative/zero position
- FM-04: Non-ACGT allele
- FM-05: Multi-allelic variant
- FM-06: Reference allele mismatch
"""

import re
from dataclasses import dataclass
from typing import Optional


@dataclass
class ValidationResult:
    """Result of a validation operation."""

    valid: bool
    error: Optional[str] = None
    warning: Optional[str] = None
    normalized: Optional[str] = None
    variant_id: Optional[str] = None


# Valid chromosomes
VALID_CHROMOSOMES = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "X", "Y", "M", "MT"
}

# Chromosome aliases
CHROMOSOME_ALIASES = {
    "MT": "M",
    "23": "X",
    "24": "Y",
    "25": "M",
}

# Valid nucleotides
VALID_NUCLEOTIDES = set("ACGT")

# ACMG criteria codes
VALID_PATHOGENIC_CODES = {
    "PVS1",
    "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
}

VALID_BENIGN_CODES = {
    "BA1",
    "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
}

ALL_VALID_ACMG_CODES = VALID_PATHOGENIC_CODES | VALID_BENIGN_CODES


def validate_chromosome(chrom: Optional[str]) -> ValidationResult:
    """Validate and normalize chromosome identifier.

    Args:
        chrom: Chromosome identifier (e.g., "17", "chr17", "X")

    Returns:
        ValidationResult with normalized chromosome or error
    """
    if chrom is None:
        return ValidationResult(
            valid=False,
            error="Chromosome cannot be None"
        )

    # Convert to string and strip whitespace
    chrom_str = str(chrom).strip().upper()

    # Remove 'CHR' prefix if present
    if chrom_str.startswith("CHR"):
        chrom_str = chrom_str[3:]

    # Apply aliases
    if chrom_str in CHROMOSOME_ALIASES:
        chrom_str = CHROMOSOME_ALIASES[chrom_str]

    # Validate
    if chrom_str not in VALID_CHROMOSOMES:
        return ValidationResult(
            valid=False,
            error=f"Invalid chromosome '{chrom}'. Valid values: 1-22, X, Y, M/MT"
        )

    return ValidationResult(
        valid=True,
        normalized=chrom_str
    )


def validate_position(pos: Optional[int | str]) -> ValidationResult:
    """Validate genomic position.

    Args:
        pos: 1-based genomic position

    Returns:
        ValidationResult with validated position or error
    """
    if pos is None:
        return ValidationResult(
            valid=False,
            error="Position cannot be None"
        )

    # Try to convert to integer
    try:
        pos_int = int(pos)
    except (ValueError, TypeError):
        return ValidationResult(
            valid=False,
            error=f"Position must be a valid integer, got '{pos}'"
        )

    # Check range
    if pos_int <= 0:
        return ValidationResult(
            valid=False,
            error=f"Position must be positive (1-based), got {pos_int}"
        )

    # Reasonable upper bound (longest human chromosome ~250Mb)
    if pos_int > 300_000_000:
        return ValidationResult(
            valid=False,
            error=f"Position {pos_int} exceeds maximum chromosome length",
            warning="Position may be from wrong genome build"
        )

    return ValidationResult(
        valid=True,
        normalized=str(pos_int)
    )


def validate_allele(allele: Optional[str]) -> ValidationResult:
    """Validate nucleotide allele.

    Args:
        allele: Reference or alternate allele (e.g., "A", "ACGT")

    Returns:
        ValidationResult with validated allele or error
    """
    if allele is None:
        return ValidationResult(
            valid=False,
            error="Allele cannot be None"
        )

    # Convert to string and strip/uppercase
    allele_str = str(allele).strip().upper()

    # Check for empty
    if not allele_str:
        return ValidationResult(
            valid=False,
            error="Allele cannot be empty"
        )

    # Check for multi-allelic (comma-separated)
    if "," in allele_str:
        return ValidationResult(
            valid=False,
            error=f"Multi-allelic variants not supported. Please decompose '{allele}' into separate variants"
        )

    # Check for deletion marker (might need special handling)
    if allele_str == "-":
        return ValidationResult(
            valid=False,
            error="Deletion marker '-' not supported. Use VCF-style representation (e.g., ref='AC', alt='A')"
        )

    # Validate nucleotides
    invalid_chars = set(allele_str) - VALID_NUCLEOTIDES
    if invalid_chars:
        return ValidationResult(
            valid=False,
            error=f"Invalid characters in allele '{allele}': {invalid_chars}. Only A, C, G, T allowed"
        )

    return ValidationResult(
        valid=True,
        normalized=allele_str
    )


def validate_variant(
    chrom: Optional[str],
    pos: Optional[int | str],
    ref: Optional[str],
    alt: Optional[str]
) -> ValidationResult:
    """Validate complete variant specification.

    Args:
        chrom: Chromosome identifier
        pos: 1-based position
        ref: Reference allele
        alt: Alternate allele

    Returns:
        ValidationResult with variant_id if valid, or error details
    """
    errors = []
    warnings = []

    # Validate each component
    chrom_result = validate_chromosome(chrom)
    if not chrom_result.valid:
        errors.append(chrom_result.error)

    pos_result = validate_position(pos)
    if not pos_result.valid:
        errors.append(pos_result.error)
    if pos_result.warning:
        warnings.append(pos_result.warning)

    ref_result = validate_allele(ref)
    if not ref_result.valid:
        errors.append(f"Reference allele: {ref_result.error}")

    alt_result = validate_allele(alt)
    if not alt_result.valid:
        errors.append(f"Alternate allele: {alt_result.error}")

    # If any validation failed, return combined errors
    if errors:
        return ValidationResult(
            valid=False,
            error="; ".join(errors),
            warning="; ".join(warnings) if warnings else None
        )

    # Additional checks with valid components
    norm_ref = ref_result.normalized
    norm_alt = alt_result.normalized

    # Check ref != alt
    if norm_ref == norm_alt:
        return ValidationResult(
            valid=False,
            error=f"Reference and alternate alleles are identical: {norm_ref}"
        )

    # Build variant ID
    norm_chrom = chrom_result.normalized
    norm_pos = pos_result.normalized
    variant_id = f"{norm_chrom}_{norm_pos}_{norm_ref}_{norm_alt}"

    return ValidationResult(
        valid=True,
        variant_id=variant_id,
        normalized=f"chr{norm_chrom}:{norm_pos} {norm_ref}>{norm_alt}",
        warning="; ".join(warnings) if warnings else None
    )


def validate_reference_against_database(
    variant_id: str,
    user_ref: str,
    stored_ref: str
) -> ValidationResult:
    """Validate user-provided reference allele against database.

    This catches coordinate/genome build mismatches where the user provides
    coordinates from the wrong build.

    Args:
        variant_id: The variant being queried
        user_ref: Reference allele provided by user
        stored_ref: Reference allele stored in database

    Returns:
        ValidationResult with error if mismatch
    """
    if not stored_ref:
        # No stored ref to compare against
        return ValidationResult(valid=True)

    user_ref_norm = user_ref.strip().upper()
    stored_ref_norm = stored_ref.strip().upper()

    if user_ref_norm != stored_ref_norm:
        return ValidationResult(
            valid=False,
            error=(
                f"Reference allele mismatch for {variant_id}: "
                f"you provided '{user_ref}', but database has '{stored_ref}'. "
                f"This may indicate wrong genome build or incorrect coordinates."
            )
        )

    return ValidationResult(valid=True)


def validate_acmg_code(code: Optional[str]) -> ValidationResult:
    """Validate ACMG/AMP evidence code.

    Args:
        code: ACMG evidence code (e.g., "PVS1", "PM2", "BP4")

    Returns:
        ValidationResult indicating if code is valid
    """
    if code is None:
        return ValidationResult(
            valid=False,
            error="ACMG code cannot be None"
        )

    code_upper = str(code).strip().upper()

    if code_upper in ALL_VALID_ACMG_CODES:
        return ValidationResult(
            valid=True,
            normalized=code_upper
        )

    # Check if it's a modified code (e.g., "PM2_Supporting")
    base_code = code_upper.split("_")[0]
    if base_code in ALL_VALID_ACMG_CODES:
        return ValidationResult(
            valid=True,
            normalized=code_upper,
            warning=f"Modified strength code detected: {code}"
        )

    return ValidationResult(
        valid=False,
        error=f"Invalid ACMG code '{code}'. Valid codes: PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7"
    )


def validate_hgvs(hgvs: str) -> ValidationResult:
    """Validate HGVS notation and extract components.

    Args:
        hgvs: HGVS string (e.g., "17:g.41197801T>A", "NC_000017.10:g.41197801T>A")

    Returns:
        ValidationResult with parsed components or error
    """
    if not hgvs:
        return ValidationResult(
            valid=False,
            error="HGVS string cannot be empty"
        )

    # Pattern for genomic HGVS: chr:g.posRef>Alt or NC_xxx:g.posRef>Alt
    patterns = [
        # Simple: 17:g.41197801T>A
        r"^(\d+|[XYM]):g\.(\d+)([ACGT]+)>([ACGT]+)$",
        # With chr: chr17:g.41197801T>A
        r"^chr(\d+|[XYM]):g\.(\d+)([ACGT]+)>([ACGT]+)$",
        # RefSeq: NC_000017.10:g.41197801T>A
        r"^NC_0000(\d+)\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)$",
    ]

    hgvs_upper = hgvs.strip().upper()

    for pattern in patterns:
        match = re.match(pattern, hgvs_upper, re.IGNORECASE)
        if match:
            chrom, pos, ref, alt = match.groups()
            return validate_variant(chrom, pos, ref, alt)

    return ValidationResult(
        valid=False,
        error=f"Invalid HGVS format: '{hgvs}'. Expected format: '17:g.41197801T>A'"
    )


def validate_gene_symbol(symbol: Optional[str]) -> ValidationResult:
    """Validate gene symbol format.

    Args:
        symbol: Gene symbol (e.g., "BRCA1", "TP53")

    Returns:
        ValidationResult with normalized symbol or error
    """
    if symbol is None:
        return ValidationResult(
            valid=False,
            error="Gene symbol cannot be None"
        )

    symbol_str = str(symbol).strip().upper()

    if not symbol_str:
        return ValidationResult(
            valid=False,
            error="Gene symbol cannot be empty"
        )

    # Basic format check (letters, numbers, hyphens)
    if not re.match(r"^[A-Z][A-Z0-9\-]*$", symbol_str):
        return ValidationResult(
            valid=False,
            error=f"Invalid gene symbol format: '{symbol}'"
        )

    return ValidationResult(
        valid=True,
        normalized=symbol_str
    )


def validate_genome_build(build: Optional[str]) -> ValidationResult:
    """Validate genome build identifier.

    Args:
        build: Genome build (e.g., "GRCh37", "hg19", "GRCh38", "hg38")

    Returns:
        ValidationResult with normalized build or error
    """
    if build is None:
        return ValidationResult(
            valid=True,  # Optional, defaults to database build
            normalized=None
        )

    build_upper = str(build).strip().upper()

    # Normalize aliases
    build_map = {
        "GRCH37": "GRCh37",
        "HG19": "GRCh37",
        "GRCH38": "GRCh38",
        "HG38": "GRCh38",
    }

    if build_upper in build_map:
        return ValidationResult(
            valid=True,
            normalized=build_map[build_upper]
        )

    return ValidationResult(
        valid=False,
        error=f"Invalid genome build '{build}'. Valid values: GRCh37/hg19, GRCh38/hg38"
    )

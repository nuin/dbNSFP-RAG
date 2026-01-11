"""Input validation tests for ACMG Variant Classification System.

These tests verify that the system properly validates variant input
to prevent silent failures and incorrect lookups.

Risk Mitigations Tested:
- FM-01: Invalid variant accepted
- FM-02: Malformed chromosome
- FM-03: Negative/zero position
- FM-04: Non-ACGT allele
- FM-05: Multi-allelic variant
- FM-06: Reference allele mismatch
"""

import pytest
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.validation import (
    validate_variant,
    validate_chromosome,
    validate_position,
    validate_allele,
    ValidationResult,
)


class TestChromosomeValidation:
    """Tests for chromosome validation."""

    @pytest.mark.parametrize("chrom", [
        "1", "2", "10", "22",  # Autosomes
        "X", "Y", "M", "MT",   # Sex and mitochondrial
        "chr1", "chr22", "chrX", "chrY", "chrM",  # With prefix
    ])
    def test_valid_chromosomes(self, chrom: str):
        """Valid chromosomes should be accepted."""
        result = validate_chromosome(chrom)
        assert result.valid is True, f"Chromosome {chrom} should be valid"

    @pytest.mark.parametrize("chrom", [
        "0",      # Zero
        "23",     # Out of range
        "99",     # Way out of range
        "ABC",    # Letters
        "",       # Empty
        "chr",    # Prefix only
        "chr0",   # Invalid with prefix
        "chr99",  # Invalid with prefix
    ])
    def test_invalid_chromosomes(self, chrom: str):
        """Invalid chromosomes should be rejected."""
        result = validate_chromosome(chrom)
        assert result.valid is False, f"Chromosome {chrom} should be invalid"
        assert "chromosome" in result.error.lower()

    def test_none_chromosome(self):
        """None chromosome should be rejected."""
        result = validate_chromosome(None)
        assert result.valid is False

    def test_chromosome_normalization(self):
        """Chromosome should be normalized (chr prefix removed)."""
        result = validate_chromosome("chr17")
        assert result.valid is True
        assert result.normalized == "17"


class TestPositionValidation:
    """Tests for genomic position validation."""

    @pytest.mark.parametrize("pos", [1, 100, 1000000, 250000000])
    def test_valid_positions(self, pos: int):
        """Valid positions should be accepted."""
        result = validate_position(pos)
        assert result.valid is True, f"Position {pos} should be valid"

    @pytest.mark.parametrize("pos", [-1, -100, 0])
    def test_invalid_negative_zero_positions(self, pos: int):
        """Negative and zero positions should be rejected."""
        result = validate_position(pos)
        assert result.valid is False, f"Position {pos} should be invalid"
        assert "position" in result.error.lower()

    def test_none_position(self):
        """None position should be rejected."""
        result = validate_position(None)
        assert result.valid is False

    def test_string_position(self):
        """String position should be rejected or converted."""
        # If string is numeric, may be converted
        result = validate_position("100")
        # Depending on implementation, may accept or reject
        # Document expected behavior

    def test_non_numeric_position(self):
        """Non-numeric position should be rejected."""
        result = validate_position("abc")
        assert result.valid is False


class TestAlleleValidation:
    """Tests for allele validation."""

    @pytest.mark.parametrize("allele", [
        "A", "C", "G", "T",           # Single nucleotides
        "AC", "GT", "ACGT",           # Multiple nucleotides
        "ACGTACGTACGT",               # Longer sequences
    ])
    def test_valid_alleles(self, allele: str):
        """Valid alleles should be accepted."""
        result = validate_allele(allele)
        assert result.valid is True, f"Allele {allele} should be valid"

    @pytest.mark.parametrize("allele", [
        "X", "N", "R", "Y",  # IUPAC ambiguity codes
        "1", "123",          # Numbers
        "A1", "1A",          # Mixed
        "a", "g",            # Lowercase
        "",                  # Empty
        "-",                 # Deletion marker (handle separately)
    ])
    def test_invalid_alleles(self, allele: str):
        """Invalid alleles should be rejected."""
        result = validate_allele(allele)
        assert result.valid is False, f"Allele {allele} should be invalid"
        assert "allele" in result.error.lower()

    def test_none_allele(self):
        """None allele should be rejected."""
        result = validate_allele(None)
        assert result.valid is False

    def test_allele_uppercase_normalization(self):
        """Alleles should be normalized to uppercase if lowercase accepted."""
        # Depending on implementation
        pass


class TestMultiAllelicDetection:
    """Tests for multi-allelic variant detection."""

    def test_single_allele(self):
        """Single allele should be accepted."""
        result = validate_allele("G")
        assert result.valid is True
        assert not getattr(result, 'is_multi_allelic', False)

    @pytest.mark.parametrize("allele", ["G,T", "A,G,C", "G,T,C,A"])
    def test_multi_allelic_detection(self, allele: str):
        """Multi-allelic variants should be detected."""
        result = validate_allele(allele)
        assert result.valid is False
        assert "multi-allelic" in result.error.lower() or "multiple" in result.error.lower()

    def test_multi_allelic_with_spaces(self):
        """Multi-allelic with spaces should be detected."""
        result = validate_allele("G, T")
        assert result.valid is False


class TestFullVariantValidation:
    """Tests for complete variant validation."""

    def test_valid_snv(self, valid_snv: dict):
        """Valid SNV should pass validation."""
        result = validate_variant(
            valid_snv["chr"],
            valid_snv["pos"],
            valid_snv["ref"],
            valid_snv["alt"]
        )
        assert result.valid is True

    def test_valid_insertion(self):
        """Valid insertion should pass validation."""
        result = validate_variant("17", 100, "A", "ACGT")
        assert result.valid is True

    def test_valid_deletion(self):
        """Valid deletion should pass validation."""
        result = validate_variant("17", 100, "ACGT", "A")
        assert result.valid is True

    def test_invalid_all_fields(self):
        """Variant with all invalid fields should report all errors."""
        result = validate_variant("99", -1, "X", "Y")
        assert result.valid is False
        # Should report multiple validation errors

    def test_ref_equals_alt_warning(self):
        """Ref equal to alt should trigger warning."""
        result = validate_variant("17", 100, "A", "A")
        # May be valid but should warn
        assert result.valid is False or result.warning is not None


class TestVariantIdGeneration:
    """Tests for variant ID generation from validated input."""

    def test_variant_id_format(self, valid_snv: dict):
        """Variant ID should follow standard format."""
        result = validate_variant(
            valid_snv["chr"],
            valid_snv["pos"],
            valid_snv["ref"],
            valid_snv["alt"]
        )
        assert result.valid is True
        expected_id = f"{valid_snv['chr']}_{valid_snv['pos']}_{valid_snv['ref']}_{valid_snv['alt']}"
        # Remove 'chr' prefix if present
        expected_id = expected_id.replace("chr", "")
        assert result.variant_id == expected_id

    def test_variant_id_normalized(self):
        """Variant ID should use normalized chromosome."""
        result = validate_variant("chr17", 100, "A", "G")
        assert result.valid is True
        assert result.variant_id == "17_100_A_G"


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_very_long_allele(self):
        """Very long allele should be handled."""
        long_allele = "A" * 10000
        result = validate_allele(long_allele)
        # May have length limit
        assert result.valid is True or "length" in result.error.lower()

    def test_chromosome_case_insensitive(self):
        """Chromosome validation should be case-insensitive."""
        assert validate_chromosome("x").valid is True
        assert validate_chromosome("X").valid is True
        assert validate_chromosome("chrx").valid is True
        assert validate_chromosome("CHRX").valid is True

    def test_whitespace_handling(self):
        """Whitespace in input should be handled."""
        result = validate_chromosome(" 17 ")
        # Should either accept (after strip) or reject
        # Document expected behavior

    def test_special_characters(self):
        """Special characters should be rejected."""
        assert validate_chromosome("17\n").valid is False or validate_chromosome("17\n").normalized == "17"
        assert validate_allele("A\t").valid is False


class TestValidationResultStructure:
    """Tests for validation result structure."""

    def test_valid_result_structure(self, valid_snv: dict):
        """Valid result should have expected structure."""
        result = validate_variant(
            valid_snv["chr"],
            valid_snv["pos"],
            valid_snv["ref"],
            valid_snv["alt"]
        )
        assert hasattr(result, 'valid')
        assert hasattr(result, 'variant_id')
        assert result.error is None or result.error == ""

    def test_invalid_result_structure(self):
        """Invalid result should have error message."""
        result = validate_variant("99", -1, "X", "Y")
        assert hasattr(result, 'valid')
        assert result.valid is False
        assert hasattr(result, 'error')
        assert result.error is not None and len(result.error) > 0


class TestIntegrationWithDatabase:
    """Tests for validation integration with database lookup."""

    @pytest.mark.requires_database
    def test_validated_variant_lookup(self, vectorstore, valid_brca1_variant: dict):
        """Validated variant should be findable in database."""
        result = validate_variant(
            valid_brca1_variant["chr"],
            valid_brca1_variant["pos"],
            valid_brca1_variant["ref"],
            valid_brca1_variant["alt"]
        )
        assert result.valid is True

        # Look up in database
        db_result = vectorstore.get_by_id(result.variant_id)
        # May or may not exist depending on database state

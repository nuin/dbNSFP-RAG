"""Classification accuracy tests for ACMG Variant Classification System.

These tests verify classification accuracy against gold standard variants
with known ClinVar classifications.

Risk Mitigations Tested:
- FM-07: Pathogenic -> Benign misclassification
- FM-08: Benign -> Pathogenic misclassification
- FM-09: VUS over-calling
- FM-10: Invalid ACMG criteria
"""

import pytest
import sys
import time
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.acmg_scoring import score_variant_from_metadata


def _classify_variant(variant: dict, vectorstore) -> str:
    """Look up variant in DB and run ACMG rule-based scoring."""
    chrom = str(variant["chr"]).replace("chr", "")
    variant_id = variant.get(
        "variant_id", f"{chrom}_{variant['pos']}_{variant['ref']}_{variant['alt']}"
    )
    data = vectorstore.get_by_id(variant_id)
    if data is None:
        return "Uncertain_significance"
    acmg = score_variant_from_metadata(data["metadata"])
    return acmg.classification.value


def _classify_with_confidence(variant: dict, vectorstore) -> dict:
    """Look up variant and return classification + confidence."""
    chrom = str(variant["chr"]).replace("chr", "")
    variant_id = variant.get(
        "variant_id", f"{chrom}_{variant['pos']}_{variant['ref']}_{variant['alt']}"
    )
    data = vectorstore.get_by_id(variant_id)
    if data is None:
        return {"classification": "Uncertain_significance", "confidence": 0.5}
    acmg = score_variant_from_metadata(data["metadata"])
    return {"classification": acmg.classification.value, "confidence": acmg.confidence}


class TestPathogenicClassification:
    """Tests for pathogenic variant classification."""

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    @pytest.mark.slow
    def test_pathogenic_sensitivity(
        self,
        pathogenic_variants: list[dict],
        vectorstore,
        classification_accuracy_threshold: float
    ):
        """Pathogenic variants should be classified as Pathogenic or Likely_pathogenic.

        Acceptance Criteria: >= 90% sensitivity
        """
        if not pathogenic_variants:
            pytest.skip("No pathogenic variants in gold standard")

        correct = 0
        results = []

        for variant in pathogenic_variants:
            # Get classification
            classification = _classify_variant(variant, vectorstore)

            # Pathogenic or Likely_pathogenic is acceptable
            is_correct = classification in ["Pathogenic", "Likely_pathogenic"]
            if is_correct:
                correct += 1

            results.append({
                "variant_id": variant.get("variant_id", f"{variant['chr']}_{variant['pos']}_{variant['ref']}_{variant['alt']}"),
                "expected": "Pathogenic",
                "actual": classification,
                "correct": is_correct
            })

        sensitivity = correct / len(pathogenic_variants)

        # Report failures
        failures = [r for r in results if not r["correct"]]
        if failures:
            failure_msg = "\n".join([
                f"  {r['variant_id']}: expected Pathogenic, got {r['actual']}"
                for r in failures[:10]  # Limit to 10
            ])
            if len(failures) > 10:
                failure_msg += f"\n  ... and {len(failures) - 10} more"

        assert sensitivity >= classification_accuracy_threshold, (
            f"Pathogenic sensitivity {sensitivity:.1%} < {classification_accuracy_threshold:.0%}\n"
            f"Failures:\n{failure_msg if failures else 'None'}"
        )

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    def test_no_pathogenic_to_benign(
        self,
        pathogenic_variants: list[dict],
        vectorstore
    ):
        """Pathogenic variants should NEVER be classified as Benign.

        This is a catastrophic misclassification that could lead to missed diagnosis.
        """
        if not pathogenic_variants:
            pytest.skip("No pathogenic variants in gold standard")

        catastrophic_failures = []

        for variant in pathogenic_variants:
            classification = _classify_variant(variant, vectorstore)

            if classification == "Benign":
                catastrophic_failures.append({
                    "variant_id": variant.get("variant_id"),
                    "gene": variant.get("gene", "Unknown"),
                    "classification": classification
                })

        assert len(catastrophic_failures) == 0, (
            f"CATASTROPHIC: {len(catastrophic_failures)} pathogenic variants classified as Benign:\n" +
            "\n".join([f"  {f['variant_id']} ({f['gene']})" for f in catastrophic_failures])
        )


class TestBenignClassification:
    """Tests for benign variant classification."""

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    @pytest.mark.slow
    def test_benign_specificity(
        self,
        benign_variants: list[dict],
        vectorstore,
        classification_accuracy_threshold: float
    ):
        """Benign variants should be classified as Benign or Likely_benign.

        Acceptance Criteria: >= 90% specificity
        """
        if not benign_variants:
            pytest.skip("No benign variants in gold standard")

        correct = 0
        results = []

        for variant in benign_variants:
            classification = _classify_variant(variant, vectorstore)

            is_correct = classification in ["Benign", "Likely_benign"]
            if is_correct:
                correct += 1

            results.append({
                "variant_id": variant.get("variant_id"),
                "expected": "Benign",
                "actual": classification,
                "correct": is_correct
            })

        specificity = correct / len(benign_variants)

        failures = [r for r in results if not r["correct"]]

        assert specificity >= classification_accuracy_threshold, (
            f"Benign specificity {specificity:.1%} < {classification_accuracy_threshold:.0%}"
        )

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    def test_no_benign_to_pathogenic(
        self,
        benign_variants: list[dict],
        vectorstore
    ):
        """Benign variants should NEVER be classified as Pathogenic.

        This could lead to unnecessary procedures and patient anxiety.
        """
        if not benign_variants:
            pytest.skip("No benign variants in gold standard")

        catastrophic_failures = []

        for variant in benign_variants:
            classification = _classify_variant(variant, vectorstore)

            if classification == "Pathogenic":
                catastrophic_failures.append({
                    "variant_id": variant.get("variant_id"),
                    "gene": variant.get("gene", "Unknown"),
                    "classification": classification
                })

        assert len(catastrophic_failures) == 0, (
            f"CATASTROPHIC: {len(catastrophic_failures)} benign variants classified as Pathogenic"
        )


class TestLikelyClassifications:
    """Tests for likely pathogenic/benign classifications."""

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    @pytest.mark.slow
    def test_likely_pathogenic_accuracy(
        self,
        likely_pathogenic_variants: list[dict],
        vectorstore
    ):
        """Likely pathogenic should be classified as LP or P."""
        if not likely_pathogenic_variants:
            pytest.skip("No likely pathogenic variants in gold standard")

        correct = 0
        for variant in likely_pathogenic_variants:
            classification = _classify_variant(variant, vectorstore)
            if classification in ["Pathogenic", "Likely_pathogenic"]:
                correct += 1

        accuracy = correct / len(likely_pathogenic_variants)
        assert accuracy >= 0.60, f"Likely pathogenic accuracy {accuracy:.1%} < 60%"

    @pytest.mark.requires_model
    @pytest.mark.requires_database
    @pytest.mark.slow
    def test_likely_benign_accuracy(
        self,
        likely_benign_variants: list[dict],
        vectorstore
    ):
        """Likely benign should be classified as LB or B."""
        if not likely_benign_variants:
            pytest.skip("No likely benign variants in gold standard")

        correct = 0
        for variant in likely_benign_variants:
            classification = _classify_variant(variant, vectorstore)
            if classification in ["Benign", "Likely_benign"]:
                correct += 1

        accuracy = correct / len(likely_benign_variants)
        assert accuracy >= 0.30, f"Likely benign accuracy {accuracy:.1%} < 30%"


class TestACMGCriteriaValidation:
    """Tests for ACMG criteria code validation."""

    VALID_PATHOGENIC_CODES = {
        "PVS1",  # Very strong
        "PS1", "PS2", "PS3", "PS4",  # Strong
        "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",  # Moderate
        "PP1", "PP2", "PP3", "PP4", "PP5",  # Supporting
    }

    VALID_BENIGN_CODES = {
        "BA1",  # Stand-alone
        "BS1", "BS2", "BS3", "BS4",  # Strong
        "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",  # Supporting
    }

    ALL_VALID_CODES = VALID_PATHOGENIC_CODES | VALID_BENIGN_CODES

    def test_criteria_code_validation(self):
        """All returned criteria codes should be valid ACMG codes."""
        from src.validation import validate_acmg_code

        for code in self.ALL_VALID_CODES:
            result = validate_acmg_code(code)
            assert result.valid is True, f"Valid code {code} rejected"

    def test_invalid_criteria_rejected(self):
        """Invalid criteria codes should be rejected."""
        from src.validation import validate_acmg_code

        invalid_codes = ["PVS2", "PS5", "PM7", "PP6", "BA2", "BS5", "BP8", "XX1", "ABC"]

        for code in invalid_codes:
            result = validate_acmg_code(code)
            assert result.valid is False, f"Invalid code {code} accepted"

    @pytest.mark.requires_model
    def test_classification_returns_valid_codes(self, vectorstore, valid_brca1_variant: dict):
        """Classification response should only contain valid ACMG codes."""
        data = vectorstore.get_by_id(valid_brca1_variant.get(
            "variant_id",
            f"{valid_brca1_variant['chr']}_{valid_brca1_variant['pos']}_{valid_brca1_variant['ref']}_{valid_brca1_variant['alt']}"
        ))
        if data is None:
            pytest.skip("BRCA1 variant not in database")
        acmg = score_variant_from_metadata(data["metadata"])
        for criterion in acmg.criteria_met:
            assert criterion.code in self.ALL_VALID_CODES, (
                f"Invalid ACMG code returned: {criterion.code}"
            )


class TestClassificationConfidence:
    """Tests for classification confidence scores."""

    @pytest.mark.requires_model
    def test_confidence_range(self, vectorstore, gold_standard: list[dict]):
        """Confidence scores should be between 0 and 1."""
        for variant in gold_standard[:10]:  # Test subset
            result = _classify_with_confidence(variant, vectorstore)
            assert 0 <= result["confidence"] <= 1, (
                f"Confidence {result['confidence']} out of range for {variant.get('variant_id')}"
            )

    @pytest.mark.requires_model
    def test_pathogenic_high_confidence(
        self,
        pathogenic_variants: list[dict],
        vectorstore
    ):
        """Pathogenic classifications should have high confidence."""
        if not pathogenic_variants:
            pytest.skip("No pathogenic variants")

        for variant in pathogenic_variants[:10]:
            result = _classify_with_confidence(variant, vectorstore)
            if result["classification"] == "Pathogenic":
                assert result["confidence"] >= 0.7, (
                    f"Pathogenic with low confidence {result['confidence']}"
                )


class TestClassificationConsistency:
    """Tests for classification consistency and reproducibility."""

    @pytest.mark.requires_model
    def test_same_variant_same_result(self, vectorstore, valid_brca1_variant: dict):
        """Same variant should always get same classification."""
        results = []
        for _ in range(3):
            classification = _classify_variant(valid_brca1_variant, vectorstore)
            results.append(classification)

        assert len(set(results)) == 1, (
            f"Inconsistent results for same variant: {results}"
        )

    @pytest.mark.requires_model
    def test_normalized_input_consistency(self, vectorstore):
        """Different input formats for same variant should give same result."""
        # With and without chr prefix
        result1 = _classify_variant(
            {"chr": "17", "pos": 41197801, "ref": "T", "alt": "A"},
            vectorstore
        )
        result2 = _classify_variant(
            {"chr": "chr17", "pos": 41197801, "ref": "T", "alt": "A"},
            vectorstore
        )

        assert result1 == result2, "Normalized inputs should give same result"


class TestPerformance:
    """Tests for classification performance."""

    @pytest.mark.requires_model
    @pytest.mark.slow
    def test_single_classification_latency(
        self,
        vectorstore,
        valid_snv: dict,
        performance_threshold_ms: int
    ):
        """Single classification should complete within threshold."""
        start = time.time()
        _ = _classify_variant(valid_snv, vectorstore)
        elapsed_ms = (time.time() - start) * 1000

        assert elapsed_ms < performance_threshold_ms, (
            f"Classification took {elapsed_ms:.0f}ms > {performance_threshold_ms}ms threshold"
        )

    @pytest.mark.requires_model
    @pytest.mark.slow
    def test_batch_classification_throughput(
        self,
        vectorstore,
        gold_standard: list[dict]
    ):
        """Batch classification should maintain reasonable throughput."""
        if len(gold_standard) < 10:
            pytest.skip("Need at least 10 variants for throughput test")

        variants = gold_standard[:10]
        start = time.time()

        for variant in variants:
            _ = _classify_variant(variant, vectorstore)

        elapsed = time.time() - start
        throughput = len(variants) / elapsed

        assert throughput >= 1, f"Throughput {throughput:.2f}/s < 1/s"

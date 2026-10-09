"""ACMG/AMP scoring public interface. For research use only."""

from src.acmg_automated import (
    evaluate_clinvar_criteria,
    evaluate_computational_criteria,
    evaluate_population_criteria,
    evaluate_variant_type_criteria,
)
from src.acmg_combination import (
    combine_criteria,
)
from src.acmg_definitions import (
    ACMG_CRITERIA_DEFINITIONS,
    create_criterion,
)
from src.acmg_evaluation import (
    evaluate_all_criteria,
)
from src.acmg_types import (
    ACMGCriterion,
    ACMGScore,
    Classification,
    CriterionStatus,
    EvidenceStrength,
    EvidenceType,
    Thresholds,
)

__all__ = [
    "ACMG_CRITERIA_DEFINITIONS",
    "ACMGCriterion",
    "ACMGScore",
    "Classification",
    "CriterionStatus",
    "EvidenceStrength",
    "EvidenceType",
    "Thresholds",
    "combine_criteria",
    "create_criterion",
    "evaluate_all_criteria",
    "evaluate_clinvar_criteria",
    "evaluate_computational_criteria",
    "evaluate_population_criteria",
    "evaluate_variant_type_criteria",
    "score_variant",
    "score_variant_from_metadata",
]


def score_variant(
    gnomad_af: float | None = None,
    gnomad_af_popmax: float | None = None,
    sift_pred: str | None = None,
    polyphen_pred: str | None = None,
    cadd_phred: float | None = None,
    revel_score: float | None = None,
    alphamissense_pred: str | None = None,
    mutation_taster: str | None = None,
    clinvar_sig: str | None = None,
    clinvar_review: str | None = None,
    consequence: str | None = None,
    variant_type: str | None = None,
    is_lof_gene: bool = True,
) -> ACMGScore:
    """
    Score a variant according to ACMG/AMP guidelines.

    This is the main entry point for variant classification.

    Args:
        gnomad_af: gnomAD allele frequency (global)
        gnomad_af_popmax: gnomAD maximum population allele frequency
        sift_pred: SIFT prediction (D=deleterious, T=tolerated)
        polyphen_pred: PolyPhen-2 prediction (D=damaging, P=possibly, B=benign)
        cadd_phred: CADD phred-scaled score
        revel_score: REVEL score (0-1)
        alphamissense_pred: AlphaMissense prediction
        mutation_taster: MutationTaster prediction
        clinvar_sig: ClinVar clinical significance
        clinvar_review: ClinVar review status
        consequence: Variant consequence (e.g., "missense_variant", "frameshift")
        variant_type: Variant type (SNV, indel, etc.)
        is_lof_gene: Whether LOF is a known disease mechanism for this gene

    Returns:
        ACMGScore with classification and criteria details
    """
    all_criteria = []

    # Population criteria
    all_criteria.extend(evaluate_population_criteria(gnomad_af, gnomad_af_popmax))

    # Computational criteria
    all_criteria.extend(
        evaluate_computational_criteria(
            sift_pred, polyphen_pred, cadd_phred, revel_score, alphamissense_pred, mutation_taster
        )
    )

    # ClinVar criteria
    all_criteria.extend(evaluate_clinvar_criteria(clinvar_sig, clinvar_review))

    # Variant type criteria
    all_criteria.extend(evaluate_variant_type_criteria(variant_type, consequence, is_lof_gene))

    # Combine and classify
    return combine_criteria(all_criteria)


def score_variant_from_metadata(metadata: dict) -> ACMGScore:
    """
    Score a variant using metadata dictionary from the vector store.

    This is a convenience function that extracts the relevant fields
    from the metadata and calls score_variant().
    """
    return score_variant(
        gnomad_af=metadata.get("gnomad_af"),
        gnomad_af_popmax=metadata.get("gnomad_af_popmax"),
        sift_pred=metadata.get("sift_pred"),
        polyphen_pred=metadata.get("polyphen_pred"),
        cadd_phred=metadata.get("cadd_phred"),
        revel_score=metadata.get("revel_score"),
        alphamissense_pred=metadata.get("alphamissense_pred"),
        mutation_taster=metadata.get("mutationtaster_pred") or metadata.get("mutation_taster"),
        clinvar_sig=metadata.get("clinvar_sig"),
        clinvar_review=metadata.get("clinvar_review"),
        consequence=metadata.get("consequence"),
        variant_type=metadata.get("variant_type"),
    )

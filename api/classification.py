"""Build classification response payloads from validated metadata."""

from fastapi import HTTPException

from api.interpretation import generate_acmg_interpretation
from api.models import ACMGCriteriaDetail, ClassificationResponse, ClinVarAnnotation, VariantRequest
from src.acmg_scoring import evaluate_all_criteria, score_variant_from_metadata
from src.evidence_links import generate_evidence_links
from src.validation import validate_reference_against_database


def _criterion_detail(criterion) -> ACMGCriteriaDetail:
    return ACMGCriteriaDetail(
        code=criterion.code,
        name=criterion.name,
        short_name=criterion.short_name,
        description=criterion.description,
        strength=criterion.strength.value,
        evidence_type=criterion.evidence_type.value,
        category=criterion.category,
        status=criterion.status.value,
        evidence=criterion.evidence,
    )


def build_classification_response(
    meta: dict, variant_id: str, request: VariantRequest, clinvar_client=None
):
    genome_build = request.genome_build

    # Validate reference allele against database
    stored_ref = meta.get("ref")
    if stored_ref:
        ref_check = validate_reference_against_database(variant_id, request.ref, stored_ref)
        if not ref_check.valid:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": ref_check.error,
                    "expected_ref": stored_ref,
                    "provided_ref": request.ref,
                    "suggestion": "Verify your coordinates match the genome build",
                },
            )

    gene = meta.get("gene")
    scores = {
        "cadd_phred": meta.get("cadd_phred"),
        "revel_score": meta.get("revel_score"),
        "gnomad_af": meta.get("gnomad_af"),
        "sift_pred": meta.get("sift_pred"),
        "polyphen_pred": meta.get("polyphen_pred"),
    }
    scores = {k: v for k, v in scores.items() if v is not None}

    # Extract ClinVar data if available
    clinvar_sig = meta.get("clinvar_sig")
    if clinvar_sig:
        clinvar_data = ClinVarAnnotation(
            id=meta.get("clinvar_id") or None,
            significance=clinvar_sig or None,
            review_status=meta.get("clinvar_review") or None,
            trait=meta.get("clinvar_trait") or None,
        )
    else:
        clinvar_data = None

    # Use rule-based ACMG scoring as PRIMARY classification
    acmg_score = score_variant_from_metadata(meta)

    criteria_met_details = [_criterion_detail(c) for c in acmg_score.criteria_met]
    all_criteria_details = [
        _criterion_detail(c) for c in evaluate_all_criteria(meta, clinvar_client=clinvar_client)
    ]

    # Generate interpretation from rule-based scoring
    clinvar_sig = meta.get("clinvar_sig")
    interpretation = generate_acmg_interpretation(
        acmg_score=acmg_score,
        gene=gene,
        variant_id=variant_id,
        scores=scores,
        clinvar_sig=clinvar_sig,
    )

    # Generate external database links for evidence verification
    evidence_links = generate_evidence_links(meta, build=genome_build.lower())

    return ClassificationResponse(
        variant_id=variant_id,
        gene=gene,
        acmg_classification=acmg_score.classification.value,
        criteria=[],  # Rule-based scoring uses criteria_met instead
        criteria_met=criteria_met_details,
        all_criteria=all_criteria_details,
        rule_applied=acmg_score.rule_applied,
        scoring_method="rule_based",
        confidence=acmg_score.confidence,
        interpretation=interpretation,
        scores=scores,
        clinvar=clinvar_data,
        evidence_links=evidence_links,
        genome_build=genome_build,
        scoring_notes=acmg_score.notes,
    )


def missing_variant_response(variant_id: str, genome_build: str) -> ClassificationResponse:
    return ClassificationResponse(
        variant_id=variant_id,
        gene=None,
        acmg_classification="Not_in_database",
        criteria=[],
        criteria_met=[],
        all_criteria=[],
        rule_applied="",
        scoring_method="rule_based",
        confidence=0.0,
        interpretation=f"Variant {variant_id} was not found in the database. Only variants within the NGSgenes panel ({326} genes) are available for classification.",
        scores={},
        clinvar=None,
        evidence_links={},
        genome_build=genome_build,
        scoring_notes=["Variant not found in database"],
    )

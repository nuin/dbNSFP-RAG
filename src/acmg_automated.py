"""Automated criteria used by the classification scoring interface."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, Thresholds


def evaluate_population_criteria(
    gnomad_af: float | None,
    gnomad_af_popmax: float | None = None,
) -> list[ACMGCriterion]:
    """
    Evaluate population-based ACMG criteria.

    BA1: AF ≥5%
    BS1: AF too high for disorder (>1%)
    PM2: Absent or extremely rare (<0.01%)
    """
    criteria = []

    # Use popmax if available, otherwise use global AF
    af = gnomad_af_popmax if gnomad_af_popmax is not None else gnomad_af

    if af is None:
        # Absent from database - meets PM2
        criteria.append(
            create_criterion(
                "PM2", met=True, evidence="Variant absent from gnomAD population database"
            )
        )
    elif af >= Thresholds.BA1_AF:
        # Very common - BA1
        criteria.append(
            create_criterion("BA1", met=True, evidence=f"Allele frequency {af:.4%} ≥ 5% in gnomAD")
        )
    elif af > Thresholds.BS1_AF:
        # Too common for most Mendelian diseases - BS1
        criteria.append(
            create_criterion("BS1", met=True, evidence=f"Allele frequency {af:.4%} > 1% in gnomAD")
        )
    elif af < Thresholds.PM2_AF:
        # Rare - PM2
        criteria.append(
            create_criterion(
                "PM2", met=True, evidence=f"Allele frequency {af:.6%} is extremely rare in gnomAD"
            )
        )

    return criteria


def evaluate_computational_criteria(
    sift_pred: str | None,
    polyphen_pred: str | None,
    cadd_phred: float | None,
    revel_score: float | None,
    alphamissense_pred: str | None = None,
    mutation_taster: str | None = None,
) -> list[ACMGCriterion]:
    """
    Evaluate computational prediction-based ACMG criteria.

    PP3: Multiple computational evidence supports deleterious
    BP4: Multiple computational evidence suggests no impact
    """
    damaging_count = 0
    benign_count = 0
    evidence_details = []

    # SIFT
    if sift_pred:
        sift_lower = sift_pred.lower()
        if "deleterious" in sift_lower or sift_lower == "d":
            damaging_count += 1
            evidence_details.append("SIFT: Deleterious")
        elif "tolerated" in sift_lower or sift_lower == "t":
            benign_count += 1
            evidence_details.append("SIFT: Tolerated")

    # PolyPhen-2
    if polyphen_pred:
        pp_lower = polyphen_pred.lower()
        if "damaging" in pp_lower or pp_lower in ("d", "p"):
            damaging_count += 1
            evidence_details.append(f"PolyPhen: {polyphen_pred}")
        elif "benign" in pp_lower or pp_lower == "b":
            benign_count += 1
            evidence_details.append("PolyPhen: Benign")

    # CADD
    if cadd_phred is not None:
        if cadd_phred >= Thresholds.CADD_DAMAGING:
            damaging_count += 1
            evidence_details.append(f"CADD phred: {cadd_phred:.1f} (≥20)")
        elif cadd_phred < Thresholds.CADD_BENIGN:
            benign_count += 1
            evidence_details.append(f"CADD phred: {cadd_phred:.1f} (<10)")

    # REVEL
    if revel_score is not None:
        if revel_score >= Thresholds.REVEL_DAMAGING:
            damaging_count += 1
            evidence_details.append(f"REVEL: {revel_score:.3f} (≥0.5)")
        elif revel_score < Thresholds.REVEL_BENIGN:
            benign_count += 1
            evidence_details.append(f"REVEL: {revel_score:.3f} (<0.25)")

    # AlphaMissense
    if alphamissense_pred:
        am_lower = alphamissense_pred.lower()
        if "pathogenic" in am_lower or am_lower == "p":
            damaging_count += 1
            evidence_details.append("AlphaMissense: Pathogenic")
        elif "benign" in am_lower or am_lower == "b":
            benign_count += 1
            evidence_details.append("AlphaMissense: Benign")

    # MutationTaster
    if mutation_taster:
        mt_lower = mutation_taster.lower()
        if "disease" in mt_lower or mt_lower in ("d", "a"):
            damaging_count += 1
            evidence_details.append("MutationTaster: Disease causing")
        elif "polymorphism" in mt_lower or mt_lower in ("n", "p"):
            benign_count += 1
            evidence_details.append("MutationTaster: Polymorphism")

    criteria = []

    # PP3: Multiple lines support damaging
    if damaging_count >= Thresholds.MIN_PREDICTORS_AGREE:
        criteria.append(
            create_criterion(
                "PP3",
                met=True,
                evidence=f"{damaging_count} computational predictors support deleterious effect: "
                + "; ".join(evidence_details),
            )
        )

    # BP4: Multiple lines support benign
    elif benign_count >= Thresholds.MIN_PREDICTORS_AGREE:
        criteria.append(
            create_criterion(
                "BP4",
                met=True,
                evidence=f"{benign_count} computational predictors suggest no impact: "
                + "; ".join(evidence_details),
            )
        )

    return criteria


def evaluate_clinvar_criteria(
    clinvar_sig: str | None,
    clinvar_review: str | None = None,
) -> list[ACMGCriterion]:
    """
    Evaluate ClinVar-based criteria (PP5/BP6).

    Note: PP5 and BP6 are controversial and should be used with caution.
    We only apply them for high-confidence ClinVar submissions.
    """
    criteria = []

    if not clinvar_sig:
        return criteria

    sig_lower = clinvar_sig.lower()

    # Check review status for confidence
    high_confidence = False
    if clinvar_review:
        review_lower = clinvar_review.lower()
        if (
            "expert" in review_lower
            or "practice" in review_lower
            or "multiple" in review_lower
            and "no_conflicts" in review_lower
        ):
            high_confidence = True

    # Only apply PP5/BP6 for high-confidence submissions
    if high_confidence:
        if "pathogenic" in sig_lower and "benign" not in sig_lower:
            if "likely" in sig_lower:
                evidence = f"ClinVar: Likely pathogenic ({clinvar_review})"
            else:
                evidence = f"ClinVar: Pathogenic ({clinvar_review})"
            criteria.append(create_criterion("PP5", met=True, evidence=evidence))

        elif "benign" in sig_lower and "pathogenic" not in sig_lower:
            if "likely" in sig_lower:
                evidence = f"ClinVar: Likely benign ({clinvar_review})"
            else:
                evidence = f"ClinVar: Benign ({clinvar_review})"
            criteria.append(create_criterion("BP6", met=True, evidence=evidence))

    return criteria


def evaluate_variant_type_criteria(
    variant_type: str | None,
    consequence: str | None,
    is_lof_gene: bool = True,  # Assume LOF is mechanism unless known otherwise
) -> list[ACMGCriterion]:
    """
    Evaluate variant type/consequence criteria.

    PVS1: LOF variant in gene where LOF is known mechanism
    BP7: Synonymous with no predicted splice impact
    """
    criteria = []

    if not consequence:
        return criteria

    cons_lower = consequence.lower()

    # PVS1: Loss-of-function variants
    lof_consequences = [
        "frameshift",
        "nonsense",
        "stop_gained",
        "stop_lost",
        "splice_donor",
        "splice_acceptor",
        "start_lost",
        "transcript_ablation",
    ]

    if any(lof in cons_lower for lof in lof_consequences) and is_lof_gene:
        criteria.append(
            create_criterion(
                "PVS1",
                met=True,
                evidence=f"Loss-of-function variant ({consequence}) in gene where LOF is a known disease mechanism",
            )
        )

    # BP7: Synonymous variants
    if "synonymous" in cons_lower:
        criteria.append(
            create_criterion(
                "BP7",
                met=True,
                evidence=f"Synonymous variant ({consequence}) - no predicted splice impact",
            )
        )

    return criteria

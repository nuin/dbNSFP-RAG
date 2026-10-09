"""Detailed computation evidence evaluation."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, CriterionStatus, Thresholds


def _predictor_counts(metadata: dict) -> tuple[int, int, list[str]]:
    sift_pred = metadata.get("sift_pred")
    polyphen_pred = metadata.get("polyphen_pred")
    cadd_phred = metadata.get("cadd_phred")
    revel_score = metadata.get("revel_score")
    alphamissense_pred = metadata.get("alphamissense_pred")
    mutation_taster = metadata.get("mutationtaster_pred") or metadata.get("mutation_taster")
    # PP3: Computational evidence supports damaging
    damaging_count = 0
    benign_count = 0
    predictor_details = []

    if sift_pred:
        if sift_pred.lower() in ("d", "deleterious"):
            damaging_count += 1
            predictor_details.append("SIFT: Deleterious")
        else:
            benign_count += 1
            predictor_details.append("SIFT: Tolerated")

    if polyphen_pred:
        if polyphen_pred.lower() in ("d", "p", "probably_damaging", "possibly_damaging"):
            damaging_count += 1
            predictor_details.append(f"PolyPhen: {polyphen_pred}")
        else:
            benign_count += 1
            predictor_details.append("PolyPhen: Benign")

    if cadd_phred is not None:
        if cadd_phred >= Thresholds.CADD_DAMAGING:
            damaging_count += 1
            predictor_details.append(f"CADD: {cadd_phred:.1f} (≥20, damaging)")
        elif cadd_phred < Thresholds.CADD_BENIGN:
            benign_count += 1
            predictor_details.append(f"CADD: {cadd_phred:.1f} (<10, benign)")
        else:
            predictor_details.append(f"CADD: {cadd_phred:.1f} (intermediate)")

    if revel_score is not None:
        if revel_score >= Thresholds.REVEL_DAMAGING:
            damaging_count += 1
            predictor_details.append(f"REVEL: {revel_score:.3f} (≥0.5, damaging)")
        elif revel_score < Thresholds.REVEL_BENIGN:
            benign_count += 1
            predictor_details.append(f"REVEL: {revel_score:.3f} (<0.25, benign)")
        else:
            predictor_details.append(f"REVEL: {revel_score:.3f} (intermediate)")

    if alphamissense_pred:
        if alphamissense_pred.lower() in ("pathogenic", "p"):
            damaging_count += 1
            predictor_details.append("AlphaMissense: Pathogenic")
        elif alphamissense_pred.lower() in ("benign", "b"):
            benign_count += 1
            predictor_details.append("AlphaMissense: Benign")
        else:
            predictor_details.append(f"AlphaMissense: {alphamissense_pred}")

    if mutation_taster:
        mt_lower = mutation_taster.lower()
        if "disease" in mt_lower or mt_lower in ("d", "a"):
            damaging_count += 1
            predictor_details.append("MutationTaster: Disease-causing")
        elif "polymorphism" in mt_lower or mt_lower in ("n", "p"):
            benign_count += 1
            predictor_details.append("MutationTaster: Polymorphism")

    return damaging_count, benign_count, predictor_details


def evaluate_computation(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    damaging_count, benign_count, predictor_details = _predictor_counts(metadata)
    if predictor_details:
        if damaging_count >= Thresholds.MIN_PREDICTORS_AGREE:
            all_criteria.append(
                create_criterion(
                    "PP3",
                    CriterionStatus.MET,
                    f"{damaging_count} of {len(predictor_details)} predictors support damaging: {'; '.join(predictor_details)}",
                )
            )
        else:
            all_criteria.append(
                create_criterion(
                    "PP3",
                    CriterionStatus.NOT_MET,
                    f"Only {damaging_count} of {len(predictor_details)} predictors support damaging (need ≥{Thresholds.MIN_PREDICTORS_AGREE}): {'; '.join(predictor_details)}",
                )
            )
    else:
        all_criteria.append(
            create_criterion(
                "PP3", CriterionStatus.NOT_EVALUATED, "No computational predictions available"
            )
        )

    # BP4: Computational evidence supports benign
    if predictor_details:
        if benign_count >= Thresholds.MIN_PREDICTORS_AGREE:
            all_criteria.append(
                create_criterion(
                    "BP4",
                    CriterionStatus.MET,
                    f"{benign_count} of {len(predictor_details)} predictors support benign: {'; '.join(predictor_details)}",
                )
            )
        else:
            all_criteria.append(
                create_criterion(
                    "BP4",
                    CriterionStatus.NOT_MET,
                    f"Only {benign_count} of {len(predictor_details)} predictors support benign (need ≥{Thresholds.MIN_PREDICTORS_AGREE}): {'; '.join(predictor_details)}",
                )
            )
    else:
        all_criteria.append(
            create_criterion(
                "BP4", CriterionStatus.NOT_EVALUATED, "No computational predictions available"
            )
        )

    return all_criteria

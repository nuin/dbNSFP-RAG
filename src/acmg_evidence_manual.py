"""Detailed manual evidence evaluation."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, CriterionStatus


def evaluate_manual_ps() -> list[ACMGCriterion]:
    all_criteria = []
    # PS2: De novo confirmed
    all_criteria.append(
        create_criterion(
            "PS2",
            CriterionStatus.NOT_EVALUATED,
            "Requires family data to confirm de novo status (not available from dbNSFP)",
        )
    )

    # PS3: Functional studies
    all_criteria.append(
        create_criterion(
            "PS3",
            CriterionStatus.NOT_EVALUATED,
            "Requires review of published functional studies (not automated)",
        )
    )

    # PS4: High prevalence in affected
    all_criteria.append(
        create_criterion(
            "PS4",
            CriterionStatus.NOT_EVALUATED,
            "Requires case-control study data (not available from dbNSFP)",
        )
    )

    return all_criteria


def evaluate_manual_pm3() -> list[ACMGCriterion]:
    all_criteria = []
    # PM3: In trans with known pathogenic
    all_criteria.append(
        create_criterion(
            "PM3",
            CriterionStatus.NOT_EVALUATED,
            "Requires phasing data to determine cis/trans status (not available)",
        )
    )

    return all_criteria


def evaluate_manual_pm6_pp1() -> list[ACMGCriterion]:
    all_criteria = []
    # PM6: Assumed de novo (without confirmation)
    all_criteria.append(
        create_criterion(
            "PM6",
            CriterionStatus.NOT_EVALUATED,
            "Requires family history/segregation data (not available from dbNSFP)",
        )
    )

    # PP1: Cosegregation with disease
    all_criteria.append(
        create_criterion(
            "PP1",
            CriterionStatus.NOT_EVALUATED,
            "Requires family segregation data (not available from dbNSFP)",
        )
    )

    return all_criteria


def evaluate_manual_pp4() -> list[ACMGCriterion]:
    all_criteria = []
    # PP4: Patient phenotype highly specific
    all_criteria.append(
        create_criterion(
            "PP4",
            CriterionStatus.NOT_EVALUATED,
            "Requires patient phenotype data (not available from dbNSFP)",
        )
    )

    return all_criteria


def evaluate_manual_bs() -> list[ACMGCriterion]:
    all_criteria = []
    # BS3: Functional studies show no damaging effect
    all_criteria.append(
        create_criterion(
            "BS3",
            CriterionStatus.NOT_EVALUATED,
            "Requires review of published functional studies (not automated)",
        )
    )

    # BS4: Lack of segregation in affected family members
    all_criteria.append(
        create_criterion(
            "BS4",
            CriterionStatus.NOT_EVALUATED,
            "Requires family segregation data (not available from dbNSFP)",
        )
    )

    return all_criteria


def evaluate_manual_bp2() -> list[ACMGCriterion]:
    all_criteria = []
    # BP2: Observed in cis/trans with pathogenic
    all_criteria.append(
        create_criterion(
            "BP2",
            CriterionStatus.NOT_EVALUATED,
            "Requires phasing data to determine cis/trans status (not available)",
        )
    )

    return all_criteria


def evaluate_manual_bp5() -> list[ACMGCriterion]:
    all_criteria = []
    # BP5: Variant found in case with alternate cause
    all_criteria.append(
        create_criterion(
            "BP5",
            CriterionStatus.NOT_EVALUATED,
            "Requires case-level data showing alternate molecular cause (not available)",
        )
    )

    return all_criteria

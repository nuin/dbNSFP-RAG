"""Evaluate and order detailed evidence for every ACMG criterion."""

from src.acmg_evidence_clinical import (
    evaluate_bp3,
    evaluate_bp6,
    evaluate_pp5,
    evaluate_ps1_pm5,
)
from src.acmg_evidence_computation import (
    evaluate_computation,
)
from src.acmg_evidence_intrinsic import (
    evaluate_bp1,
    evaluate_bp7,
    evaluate_pm1,
    evaluate_pm4,
    evaluate_pp2,
    evaluate_pvs1,
)
from src.acmg_evidence_manual import (
    evaluate_manual_bp2,
    evaluate_manual_bp5,
    evaluate_manual_bs,
    evaluate_manual_pm3,
    evaluate_manual_pm6_pp1,
    evaluate_manual_pp4,
    evaluate_manual_ps,
)
from src.acmg_evidence_population import (
    evaluate_ba1,
    evaluate_bs1,
    evaluate_bs2,
    evaluate_pm2,
)
from src.acmg_types import ACMGCriterion

CRITERIA_DISPLAY_ORDER = (
    "PVS1",
    "PS2",
    "PS3",
    "PS4",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PS1",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
)


def evaluate_all_criteria(metadata: dict, clinvar_client=None) -> list[ACMGCriterion]:
    """Evaluate all 28 criteria with status and supporting evidence."""
    criteria = []
    criteria.extend(evaluate_pvs1(metadata))
    criteria.extend(evaluate_pm1(metadata))
    criteria.extend(evaluate_pm4(metadata))
    criteria.extend(evaluate_pp2(metadata))
    criteria.extend(evaluate_bp1(metadata))
    criteria.extend(evaluate_bp7(metadata))
    criteria.extend(evaluate_pm2(metadata))
    criteria.extend(evaluate_ba1(metadata))
    criteria.extend(evaluate_bs1(metadata))
    criteria.extend(evaluate_bs2(metadata))
    criteria.extend(evaluate_ps1_pm5(metadata, clinvar_client))
    criteria.extend(evaluate_pp5(metadata))
    criteria.extend(evaluate_bp6(metadata))
    criteria.extend(evaluate_bp3(metadata))
    criteria.extend(evaluate_computation(metadata))
    criteria.extend(evaluate_manual_ps())
    criteria.extend(evaluate_manual_pm3())
    criteria.extend(evaluate_manual_pm6_pp1())
    criteria.extend(evaluate_manual_pp4())
    criteria.extend(evaluate_manual_bs())
    criteria.extend(evaluate_manual_bp2())
    criteria.extend(evaluate_manual_bp5())
    by_code = {criterion.code: criterion for criterion in criteria}
    return [by_code[code] for code in CRITERIA_DISPLAY_ORDER]

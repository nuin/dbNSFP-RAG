"""Detailed population evidence evaluation."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, CriterionStatus, Thresholds


def evaluate_pm2(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    gnomad_af = metadata.get("gnomad_af")
    # PM2: Absent/rare in population
    if gnomad_af is None:
        all_criteria.append(
            create_criterion(
                "PM2", CriterionStatus.MET, "Variant absent from gnomAD population database"
            )
        )
    elif gnomad_af < Thresholds.PM2_AF:
        all_criteria.append(
            create_criterion(
                "PM2",
                CriterionStatus.MET,
                f"Allele frequency {gnomad_af:.6%} is extremely rare (< 0.01%) in gnomAD",
            )
        )
    else:
        all_criteria.append(
            create_criterion(
                "PM2",
                CriterionStatus.NOT_MET,
                f"Allele frequency {gnomad_af:.4%} exceeds rare threshold (0.01%)",
            )
        )

    return all_criteria


def evaluate_ba1(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    gnomad_af = metadata.get("gnomad_af")
    # BA1: Very high population frequency (≥5%)
    if gnomad_af is not None and gnomad_af >= Thresholds.BA1_AF:
        all_criteria.append(
            create_criterion(
                "BA1",
                CriterionStatus.MET,
                f"Allele frequency {gnomad_af:.4%} is ≥5% in gnomAD (stand-alone benign)",
            )
        )
    elif gnomad_af is not None:
        all_criteria.append(
            create_criterion(
                "BA1",
                CriterionStatus.NOT_MET,
                f"Allele frequency {gnomad_af:.4%} is below 5% threshold",
            )
        )
    else:
        all_criteria.append(
            create_criterion(
                "BA1", CriterionStatus.NOT_MET, "Variant not found in gnomAD (cannot be common)"
            )
        )

    return all_criteria


def evaluate_bs1(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    gnomad_af = metadata.get("gnomad_af")
    # BS1: Too common for disorder (>1%)
    if gnomad_af is not None and gnomad_af > Thresholds.BS1_AF:
        all_criteria.append(
            create_criterion(
                "BS1",
                CriterionStatus.MET,
                f"Allele frequency {gnomad_af:.4%} exceeds 1% (too common for most Mendelian disorders)",
            )
        )
    elif gnomad_af is not None:
        all_criteria.append(
            create_criterion(
                "BS1",
                CriterionStatus.NOT_MET,
                f"Allele frequency {gnomad_af:.4%} is below 1% threshold",
            )
        )
    else:
        all_criteria.append(
            create_criterion("BS1", CriterionStatus.NOT_MET, "Variant not found in gnomAD")
        )

    return all_criteria


def evaluate_bs2(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    # BS2: Observed in healthy adult (for penetrant disorder)
    # Partially evaluable using gnomAD homozygote count (for recessive disorders)
    gnomad_hom = metadata.get("gnomad_hom")
    if gnomad_hom is not None and gnomad_hom > 0:
        all_criteria.append(
            create_criterion(
                "BS2",
                CriterionStatus.MET,
                f"Observed as homozygous in {int(gnomad_hom)} healthy individual(s) in gnomAD (supports benign for recessive disorders)",
            )
        )
    else:
        all_criteria.append(
            create_criterion(
                "BS2", CriterionStatus.NOT_MET, "No homozygous individuals observed in gnomAD"
            )
        )

    return all_criteria

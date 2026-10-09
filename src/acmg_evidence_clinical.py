"""Detailed clinical evidence evaluation."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, CriterionStatus


def evaluate_ps1_pm5(metadata: dict, clinvar_client=None) -> list[ACMGCriterion]:
    all_criteria = []
    # PS1 and PM5: Evaluate using ClinVar API if available
    gene = metadata.get("gene", "")
    aa_ref = metadata.get("aa_ref", "")  # May need to extract from variant data
    aa_alt = metadata.get("aa_alt", "")
    aa_pos = metadata.get("aa_pos")

    if clinvar_client and gene:
        try:
            from src.clinvar_api import evaluate_ps1_pm5

            # Convert aa_pos to int if it's a string
            if isinstance(aa_pos, str) and aa_pos.isdigit():
                aa_pos = int(aa_pos)

            if aa_pos and aa_ref and aa_alt:
                ps1_met, pm5_met, evidence = evaluate_ps1_pm5(
                    clinvar_client, gene, aa_pos, aa_ref, aa_alt
                )

                if ps1_met:
                    all_criteria.append(create_criterion("PS1", CriterionStatus.MET, evidence))
                else:
                    all_criteria.append(
                        create_criterion(
                            "PS1",
                            CriterionStatus.NOT_MET,
                            "No pathogenic variant with same amino acid change found in ClinVar",
                        )
                    )

                if pm5_met:
                    all_criteria.append(create_criterion("PM5", CriterionStatus.MET, evidence))
                else:
                    all_criteria.append(
                        create_criterion(
                            "PM5",
                            CriterionStatus.NOT_MET,
                            "No pathogenic variant at same position with different AA change in ClinVar",
                        )
                    )
            else:
                all_criteria.append(
                    create_criterion(
                        "PS1",
                        CriterionStatus.NOT_EVALUATED,
                        "Amino acid change data incomplete (missing position or ref/alt AA)",
                    )
                )
                all_criteria.append(
                    create_criterion(
                        "PM5",
                        CriterionStatus.NOT_EVALUATED,
                        "Amino acid change data incomplete (missing position or ref/alt AA)",
                    )
                )
        except (ImportError, ValueError, TypeError, KeyError, OSError) as e:
            all_criteria.append(
                create_criterion(
                    "PS1", CriterionStatus.NOT_EVALUATED, f"ClinVar lookup failed: {e}"
                )
            )
            all_criteria.append(
                create_criterion(
                    "PM5", CriterionStatus.NOT_EVALUATED, f"ClinVar lookup failed: {e}"
                )
            )
    else:
        all_criteria.append(
            create_criterion(
                "PS1",
                CriterionStatus.NOT_EVALUATED,
                "Requires ClinVar API client for automated evaluation (not configured)",
            )
        )
        all_criteria.append(
            create_criterion(
                "PM5",
                CriterionStatus.NOT_EVALUATED,
                "Requires ClinVar API client for automated evaluation (not configured)",
            )
        )

    return all_criteria


def evaluate_pp5(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    clinvar_sig = metadata.get("clinvar_sig")
    clinvar_review = metadata.get("clinvar_review")
    # PP5: Reputable source reports pathogenic
    if clinvar_sig:
        sig_lower = clinvar_sig.lower()
        has_high_review = clinvar_review and (
            "expert" in clinvar_review.lower() or "multiple" in clinvar_review.lower()
        )
        if "pathogenic" in sig_lower and "benign" not in sig_lower:
            if has_high_review:
                all_criteria.append(
                    create_criterion(
                        "PP5",
                        CriterionStatus.MET,
                        f"ClinVar reports {clinvar_sig} with review status: {clinvar_review}",
                    )
                )
            else:
                all_criteria.append(
                    create_criterion(
                        "PP5",
                        CriterionStatus.NOT_MET,
                        f"ClinVar reports {clinvar_sig} but review status ({clinvar_review}) not sufficient for PP5",
                    )
                )
        else:
            all_criteria.append(
                create_criterion(
                    "PP5",
                    CriterionStatus.NOT_MET,
                    f"ClinVar significance ({clinvar_sig}) is not pathogenic",
                )
            )
    else:
        all_criteria.append(
            create_criterion(
                "PP5", CriterionStatus.NOT_EVALUATED, "No ClinVar annotation available"
            )
        )

    return all_criteria


def evaluate_bp6(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    clinvar_sig = metadata.get("clinvar_sig")
    clinvar_review = metadata.get("clinvar_review")
    # BP6: Reputable source reports benign
    if clinvar_sig:
        sig_lower = clinvar_sig.lower()
        has_high_review = clinvar_review and (
            "expert" in clinvar_review.lower() or "multiple" in clinvar_review.lower()
        )
        if "benign" in sig_lower and "pathogenic" not in sig_lower:
            if has_high_review:
                all_criteria.append(
                    create_criterion(
                        "BP6",
                        CriterionStatus.MET,
                        f"ClinVar reports {clinvar_sig} with review status: {clinvar_review}",
                    )
                )
            else:
                all_criteria.append(
                    create_criterion(
                        "BP6",
                        CriterionStatus.NOT_MET,
                        f"ClinVar reports {clinvar_sig} but review status ({clinvar_review}) not sufficient for BP6",
                    )
                )
        else:
            all_criteria.append(
                create_criterion(
                    "BP6",
                    CriterionStatus.NOT_MET,
                    f"ClinVar significance ({clinvar_sig}) is not benign",
                )
            )
    else:
        all_criteria.append(
            create_criterion(
                "BP6", CriterionStatus.NOT_EVALUATED, "No ClinVar annotation available"
            )
        )

    return all_criteria


def evaluate_bp3(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    # BP3: In-frame indel in repetitive region
    is_inframe = "inframe" in cons_lower or "in_frame" in cons_lower
    if is_inframe:
        # Check if in repeat region using UCSC API
        try:
            from src.repeat_api import is_in_repeat_region

            chrom = metadata.get("chr", "")
            pos = metadata.get("pos", 0)
            if chrom and pos:
                in_repeat, repeat_class, repeat_family = is_in_repeat_region(chrom, pos)
                if in_repeat:
                    all_criteria.append(
                        create_criterion(
                            "BP3",
                            CriterionStatus.MET,
                            f"In-frame indel in {repeat_class}/{repeat_family} repetitive region",
                        )
                    )
                else:
                    all_criteria.append(
                        create_criterion(
                            "BP3",
                            CriterionStatus.NOT_MET,
                            "In-frame indel but not in a repetitive region",
                        )
                    )
            else:
                all_criteria.append(
                    create_criterion(
                        "BP3",
                        CriterionStatus.NOT_EVALUATED,
                        "Position data not available for repeat region lookup",
                    )
                )
        except (ImportError, ValueError, TypeError, KeyError, OSError) as e:
            all_criteria.append(
                create_criterion(
                    "BP3", CriterionStatus.NOT_EVALUATED, f"Repeat region lookup failed: {e}"
                )
            )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "BP3",
                CriterionStatus.NOT_MET,
                f"Not an in-frame variant (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion("BP3", CriterionStatus.NOT_EVALUATED, "Consequence type not available")
        )

    return all_criteria

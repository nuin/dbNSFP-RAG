"""Detailed intrinsic evidence evaluation."""

from src.acmg_definitions import create_criterion
from src.acmg_types import ACMGCriterion, CriterionStatus


def evaluate_pvs1(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    # PVS1: Null variant in LOF gene (enhanced with gene constraint)
    lof_consequences = [
        "frameshift",
        "nonsense",
        "stop_gained",
        "stop_lost",
        "splice_donor",
        "splice_acceptor",
        "start_lost",
    ]
    is_lof = any(lof in cons_lower for lof in lof_consequences)
    loeuf = metadata.get("loeuf")
    gnomad_pli = metadata.get("gnomad_pli")

    if is_lof:
        # Check gene constraint to strengthen/weaken evidence
        constraint_info = ""
        if loeuf is not None:
            if loeuf < 0.35:
                constraint_info = f" Gene is highly LOF-intolerant (LOEUF={loeuf:.2f} < 0.35)."
            elif loeuf < 0.6:
                constraint_info = f" Gene is moderately LOF-constrained (LOEUF={loeuf:.2f})."
            else:
                constraint_info = f" Gene tolerates some LOF variation (LOEUF={loeuf:.2f})."
        elif gnomad_pli is not None:
            if gnomad_pli >= 0.9:
                constraint_info = f" Gene is LOF-intolerant (pLI={gnomad_pli:.2f} ≥ 0.9)."
            else:
                constraint_info = f" Gene has moderate LOF tolerance (pLI={gnomad_pli:.2f})."

        all_criteria.append(
            create_criterion(
                "PVS1",
                CriterionStatus.MET,
                f"Loss-of-function variant ({consequence}) in gene where LOF is a known disease mechanism.{constraint_info}",
            )
        )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "PVS1",
                CriterionStatus.NOT_MET,
                f"Not a loss-of-function variant (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion(
                "PVS1", CriterionStatus.NOT_EVALUATED, "Consequence type not available"
            )
        )

    # PS1: Same amino acid change as known pathogenic
    # (Evaluated together with PM5 below after PM4)

    return all_criteria


def evaluate_pm1(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    # PM1: Mutational hotspot / functional domain
    interpro_domain = metadata.get("interpro_domain", "")
    if interpro_domain and interpro_domain not in ("", ".", "N/A"):
        all_criteria.append(
            create_criterion(
                "PM1", CriterionStatus.MET, f"Located in functional domain: {interpro_domain}"
            )
        )
    else:
        all_criteria.append(
            create_criterion(
                "PM1", CriterionStatus.NOT_MET, "Variant not located in annotated functional domain"
            )
        )

    return all_criteria


def evaluate_pm4(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    # PM4: Protein length change (in-frame indel, stop-loss)
    if "inframe" in cons_lower or "stop_lost" in cons_lower:
        all_criteria.append(
            create_criterion(
                "PM4", CriterionStatus.MET, f"In-frame protein length change ({consequence})"
            )
        )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "PM4",
                CriterionStatus.NOT_MET,
                f"Not an in-frame protein length change (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion("PM4", CriterionStatus.NOT_EVALUATED, "Consequence type not available")
        )

    return all_criteria


def evaluate_pp2(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    # PP2: Missense in gene with low benign missense rate
    # Using gnomad_mis_oe (missense observed/expected ratio) - lower is more constrained
    gnomad_mis_oe = metadata.get("gnomad_mis_oe")
    if "missense" in cons_lower:
        if gnomad_mis_oe is not None:
            # Low missense o/e ratio (<0.6) indicates gene intolerant to missense
            if gnomad_mis_oe < 0.6:
                all_criteria.append(
                    create_criterion(
                        "PP2",
                        CriterionStatus.MET,
                        f"Missense variant in gene with high missense constraint (gnomAD mis o/e = {gnomad_mis_oe:.2f} < 0.6)",
                    )
                )
            else:
                all_criteria.append(
                    create_criterion(
                        "PP2",
                        CriterionStatus.NOT_MET,
                        f"Gene not highly missense-constrained (gnomAD mis o/e = {gnomad_mis_oe:.2f} ≥ 0.6)",
                    )
                )
        else:
            all_criteria.append(
                create_criterion(
                    "PP2",
                    CriterionStatus.NOT_EVALUATED,
                    "Missense variant - gene constraint score not available",
                )
            )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "PP2",
                CriterionStatus.NOT_MET,
                f"Not a missense variant (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion("PP2", CriterionStatus.NOT_EVALUATED, "Consequence type not available")
        )

    return all_criteria


def evaluate_bp1(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    gnomad_mis_oe = metadata.get("gnomad_mis_oe")
    # BP1: Missense in gene where only truncating causes disease
    gnomad_pli = metadata.get("gnomad_pli")
    if "missense" in cons_lower:
        if gnomad_pli is not None:
            # High pLI (>0.9) + high missense o/e (>0.8) suggests truncating is primary mechanism
            # (LOF intolerant but missense tolerated)
            if gnomad_pli >= 0.9 and (gnomad_mis_oe is None or gnomad_mis_oe > 0.8):
                mis_oe_str = f"{gnomad_mis_oe:.2f}" if gnomad_mis_oe else "N/A"
                all_criteria.append(
                    create_criterion(
                        "BP1",
                        CriterionStatus.MET,
                        f"Missense in gene where truncating variants are primary disease mechanism (pLI = {gnomad_pli:.2f}, mis o/e = {mis_oe_str})",
                    )
                )
            else:
                all_criteria.append(
                    create_criterion(
                        "BP1",
                        CriterionStatus.NOT_MET,
                        f"Gene constraint does not suggest truncating-only mechanism (pLI = {gnomad_pli:.2f})",
                    )
                )
        else:
            all_criteria.append(
                create_criterion(
                    "BP1", CriterionStatus.NOT_EVALUATED, "Gene constraint scores not available"
                )
            )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "BP1",
                CriterionStatus.NOT_MET,
                f"Not a missense variant (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion("BP1", CriterionStatus.NOT_EVALUATED, "Consequence type not available")
        )

    return all_criteria


def evaluate_bp7(metadata: dict) -> list[ACMGCriterion]:
    all_criteria = []
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""
    # BP7: Synonymous with no splice impact (enhanced with SpliceAI)
    if "synonymous" in cons_lower:
        # Get SpliceAI scores and compute max delta
        spliceai_ag = metadata.get("spliceai_ag")
        spliceai_al = metadata.get("spliceai_al")
        spliceai_dg = metadata.get("spliceai_dg")
        spliceai_dl = metadata.get("spliceai_dl")

        spliceai_scores = [
            s for s in [spliceai_ag, spliceai_al, spliceai_dg, spliceai_dl] if s is not None
        ]
        spliceai_max = max(spliceai_scores) if spliceai_scores else None

        if spliceai_max is not None:
            if spliceai_max < 0.2:
                all_criteria.append(
                    create_criterion(
                        "BP7",
                        CriterionStatus.MET,
                        f"Synonymous variant ({consequence}) with no predicted splice impact (SpliceAI max={spliceai_max:.2f} < 0.2)",
                    )
                )
            else:
                # Synonymous but potential splice impact - DO NOT apply BP7
                all_criteria.append(
                    create_criterion(
                        "BP7",
                        CriterionStatus.NOT_MET,
                        f"Synonymous variant but potential splice impact (SpliceAI max={spliceai_max:.2f} ≥ 0.2) - BP7 not applied",
                    )
                )
        else:
            # No SpliceAI data available - apply with caution
            all_criteria.append(
                create_criterion(
                    "BP7",
                    CriterionStatus.MET,
                    f"Synonymous variant ({consequence}) - no amino acid change (SpliceAI not available)",
                )
            )
    elif consequence:
        all_criteria.append(
            create_criterion(
                "BP7",
                CriterionStatus.NOT_MET,
                f"Not a synonymous variant (consequence: {consequence})",
            )
        )
    else:
        all_criteria.append(
            create_criterion("BP7", CriterionStatus.NOT_EVALUATED, "Consequence type not available")
        )

    return all_criteria

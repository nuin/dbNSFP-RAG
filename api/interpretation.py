"""Format classifications as human-readable evidence summaries."""


def generate_acmg_interpretation(
    acmg_score, gene: str, variant_id: str, scores: dict, clinvar_sig: str | None = None
) -> str:
    """Generate a clinical interpretation from ACMG rule-based scoring.

    Creates a clear, evidence-based interpretation that aligns with the
    ACMG/AMP guidelines and the criteria actually applied.
    """
    classification = acmg_score.classification.value.replace("_", " ")
    criteria_met = acmg_score.criteria_met
    rule = acmg_score.rule_applied

    lines = []

    # Classification summary
    lines.append(f"**ACMG Classification: {classification}**")
    lines.append("")

    if rule:
        lines.append(f"*Classification Rule: {rule}*")
        lines.append("")

    # Gene context
    if gene:
        lines.append(
            f"This variant in **{gene}** has been classified as **{classification}** based on the following evidence:"
        )
    else:
        lines.append(
            f"This variant has been classified as **{classification}** based on the following evidence:"
        )
    lines.append("")

    # Evidence summary by category
    if criteria_met:
        pathogenic_criteria = [c for c in criteria_met if c.evidence_type.value == "pathogenic"]
        benign_criteria = [c for c in criteria_met if c.evidence_type.value == "benign"]

        if pathogenic_criteria:
            lines.append("**Pathogenic Evidence:**")
            for c in pathogenic_criteria:
                strength_label = c.strength.value.replace("_", " ").title()
                lines.append(f"- **{c.code}** ({strength_label}): {c.name}")
                if c.evidence:
                    lines.append(f"  - {c.evidence}")
            lines.append("")

        if benign_criteria:
            lines.append("**Benign Evidence:**")
            for c in benign_criteria:
                strength_label = c.strength.value.replace("_", " ").title()
                lines.append(f"- **{c.code}** ({strength_label}): {c.name}")
                if c.evidence:
                    lines.append(f"  - {c.evidence}")
            lines.append("")
    else:
        lines.append("No specific ACMG criteria were met based on available evidence.")
        lines.append("")

    # ClinVar concordance
    if clinvar_sig:
        lines.append("**ClinVar Annotation:**")
        lines.append(f"ClinVar reports this variant as: {clinvar_sig}")
        clinvar_lower = clinvar_sig.lower()
        classification_lower = classification.lower()
        if clinvar_lower in classification_lower or classification_lower in clinvar_lower:
            lines.append("✓ Rule-based classification is concordant with ClinVar.")
        elif (
            "pathogenic" in clinvar_lower
            and "benign" in classification_lower
            or "benign" in clinvar_lower
            and "pathogenic" in classification_lower
        ):
            lines.append(
                "⚠ Discordance: Rule-based classification differs from ClinVar. Expert review recommended."
            )
        lines.append("")

    # Key scores summary
    if scores:
        lines.append("**Key Scores:**")
        score_items = []
        if "gnomad_af" in scores and scores["gnomad_af"] is not None:
            af = scores["gnomad_af"]
            if af >= 0.01:
                score_items.append(f"gnomAD AF: {af * 100:.2f}% (common)")
            elif af >= 0.001:
                score_items.append(f"gnomAD AF: {af * 100:.3f}% (low frequency)")
            else:
                score_items.append(f"gnomAD AF: {af:.2e} (rare)")
        if "cadd_phred" in scores and scores["cadd_phred"] is not None:
            cadd = scores["cadd_phred"]
            label = "high" if cadd >= 20 else "moderate" if cadd >= 10 else "low"
            score_items.append(f"CADD: {cadd:.1f} ({label})")
        if "revel_score" in scores and scores["revel_score"] is not None:
            revel = scores["revel_score"]
            label = "pathogenic" if revel >= 0.5 else "uncertain" if revel >= 0.25 else "benign"
            score_items.append(f"REVEL: {revel:.3f} ({label})")
        if score_items:
            lines.append(", ".join(score_items))
        lines.append("")

    # Disclaimer
    lines.append("---")
    lines.append(
        "*This classification is based on automated ACMG/AMP rule application. Expert review is required before clinical use.*"
    )

    return "\n".join(lines)

"""Combine ACMG evidence using the existing ordered classification rules."""

from collections import Counter

from src.acmg_types import ACMGCriterion, ACMGScore, Classification


def _benign_result(counts: Counter, has_pathogenic: bool):
    rules = (
        (
            counts["BA"] >= 1,
            Classification.BENIGN,
            0.95,
            "BA1: Allele frequency ≥5%",
            "Warning: Contradictory pathogenic evidence present but BA1 overrides",
        ),
        (
            counts["BS"] >= 2,
            Classification.BENIGN,
            0.90,
            "≥2 BS criteria met",
            "Warning: Contradictory pathogenic evidence present",
        ),
        (
            counts["BS"] >= 1 and counts["BP"] >= 1,
            Classification.LIKELY_BENIGN,
            0.80,
            "1 BS + 1 BP",
            "Contradictory pathogenic evidence - review recommended",
        ),
        (
            counts["BP"] >= 2,
            Classification.LIKELY_BENIGN,
            0.75,
            "≥2 BP criteria met",
            "Contradictory pathogenic evidence - review recommended",
        ),
    )
    for matches, classification, confidence, rule, warning in rules:
        if matches:
            return classification, confidence, rule, [warning] if has_pathogenic else []
    return None


def _pathogenic_result(counts: Counter):
    pvs, ps, pm, pp = (counts[prefix] for prefix in ("PVS", "PS", "PM", "PP"))
    pathogenic = Classification.PATHOGENIC
    likely = Classification.LIKELY_PATHOGENIC
    rules = (
        (pvs >= 1 and ps >= 1, pathogenic, 0.95, "PVS1 + ≥1 PS"),
        (pvs >= 1 and pm >= 2, pathogenic, 0.90, "PVS1 + ≥2 PM"),
        (pvs >= 1 and pm >= 1 and pp >= 1, pathogenic, 0.85, "PVS1 + 1 PM + 1 PP"),
        (pvs >= 1 and pp >= 2, pathogenic, 0.80, "PVS1 + ≥2 PP"),
        (pvs >= 1 and pm >= 1, likely, 0.80, "PVS1 + 1 PM"),
        (ps >= 2, pathogenic, 0.90, "≥2 PS"),
        (ps >= 1 and pm >= 3, pathogenic, 0.85, "1 PS + ≥3 PM"),
        (ps >= 1 and pm >= 2 and pp >= 2, pathogenic, 0.80, "1 PS + 2 PM + ≥2 PP"),
        (ps >= 1 and pm >= 1 and pp >= 4, pathogenic, 0.80, "1 PS + 1 PM + ≥4 PP"),
        (ps >= 1 and pm >= 1, likely, 0.75, "1 PS + 1-2 PM"),
        (ps >= 1 and pp >= 2, likely, 0.70, "1 PS + ≥2 PP"),
        (pm >= 3, likely, 0.75, "≥3 PM"),
        (pm >= 2 and pp >= 2, likely, 0.70, "2 PM + ≥2 PP"),
        (pm >= 1 and pp >= 4, likely, 0.70, "1 PM + ≥4 PP"),
    )
    for matches, classification, confidence, rule in rules:
        if matches:
            return classification, confidence, rule, []
    return Classification.UNCERTAIN_SIGNIFICANCE, 0.5, "", []


def combine_criteria(criteria: list[ACMGCriterion]) -> ACMGScore:
    """Apply benign rules first, then the ordered pathogenic combinations."""
    met_criteria = [criterion for criterion in criteria if criterion.met]
    prefixes = ("PVS", "PS", "PM", "PP", "BA", "BS", "BP")
    counts = Counter(
        {
            prefix: sum(criterion.code.startswith(prefix) for criterion in met_criteria)
            for prefix in prefixes
        }
    )
    has_pathogenic = any(counts[prefix] for prefix in ("PVS", "PS", "PM", "PP"))
    has_benign = any(counts[prefix] for prefix in ("BA", "BS", "BP"))
    result = _benign_result(counts, has_pathogenic)
    if result is None:
        result = _pathogenic_result(counts)
    classification, confidence, rule, notes = result
    if has_pathogenic and has_benign and classification == Classification.UNCERTAIN_SIGNIFICANCE:
        notes.append("Contradictory pathogenic and benign evidence - classified as VUS")
    return ACMGScore(
        classification=classification,
        criteria_met=met_criteria,
        criteria_not_met=[criterion for criterion in criteria if not criterion.met],
        confidence=confidence,
        rule_applied=rule,
        notes=notes,
    )

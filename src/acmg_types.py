"""Shared ACMG evidence and classification types."""

from dataclasses import dataclass, field
from enum import Enum


class Classification(Enum):
    """ACMG 5-tier classification."""

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain_significance"
    LIKELY_BENIGN = "Likely_benign"
    BENIGN = "Benign"


class EvidenceStrength(Enum):
    """Evidence strength levels."""

    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    STAND_ALONE = "stand_alone"


class EvidenceType(Enum):
    """Evidence type (pathogenic or benign)."""

    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class CriterionStatus(Enum):
    """Evaluation status for a criterion."""

    MET = "met"  # Criterion was met
    NOT_MET = "not_met"  # Criterion was evaluated but not met
    NOT_EVALUATED = "not_evaluated"  # Cannot be evaluated from available data


@dataclass
class ACMGCriterion:
    """Represents a single ACMG evidence criterion."""

    code: str
    name: str
    description: str
    strength: EvidenceStrength
    evidence_type: EvidenceType
    status: CriterionStatus = CriterionStatus.NOT_EVALUATED
    met: bool = False  # Convenience - True if status == MET
    evidence: str | None = None  # Explanation of why criterion was/wasn't met

    # For UI display
    short_name: str = ""
    category: str = ""  # population, computation, intrinsic, clinical, literature


@dataclass
class ACMGScore:
    """Result of ACMG classification scoring."""

    classification: Classification
    criteria_met: list[ACMGCriterion] = field(default_factory=list)
    criteria_not_met: list[ACMGCriterion] = field(default_factory=list)
    confidence: float = 0.5
    rule_applied: str = ""
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "classification": self.classification.value,
            "criteria_met": [
                {
                    "code": c.code,
                    "name": c.name,
                    "short_name": c.short_name,
                    "description": c.description,
                    "strength": c.strength.value,
                    "evidence_type": c.evidence_type.value,
                    "category": c.category,
                    "evidence": c.evidence,
                }
                for c in self.criteria_met
            ],
            "confidence": self.confidence,
            "rule_applied": self.rule_applied,
            "notes": self.notes,
        }


class Thresholds:
    """Thresholds for automated ACMG criteria evaluation."""

    # Population frequency thresholds
    BA1_AF = 0.05  # ≥5% = Benign Stand-Alone
    BS1_AF = 0.01  # >1% = too common (disease-dependent, using conservative)
    PM2_AF = 0.0001  # <0.01% = rare (absent from controls)

    # Computational prediction thresholds
    CADD_DAMAGING = 20  # CADD phred ≥20 suggests damaging
    CADD_BENIGN = 10  # CADD phred <10 suggests benign

    REVEL_DAMAGING = 0.5  # REVEL ≥0.5 suggests damaging
    REVEL_BENIGN = 0.25  # REVEL <0.25 suggests benign

    # Number of computational predictors needed for PP3/BP4
    MIN_PREDICTORS_AGREE = 3

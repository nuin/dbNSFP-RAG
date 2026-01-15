"""
ACMG/AMP Variant Classification Scoring Module

Implements the official ACMG/AMP 2015 guidelines for variant classification.
Reference: Richards et al. Genetics in Medicine (2015) 17:405-424

This module evaluates evidence criteria and combines them according to the
official scoring rules to produce 5-tier classifications.

## Criteria Implementation Status

### Fully Automated (from dbNSFP data):
- PVS1: LOF variant detection with LOEUF/pLI gene constraint scoring
- PM1: Functional domain detection via InterPro annotations
- PM2: Rare variant detection using gnomAD allele frequency
- PM4: In-frame protein length change detection
- PP2: Missense in constrained gene using gnomAD mis_z score
- PP3: Computational predictor consensus (SIFT, PolyPhen, CADD, REVEL, AlphaMissense, MutationTaster)
- PP5: ClinVar pathogenic assertion with review status check
- BA1: Common variant (≥5% AF) detection
- BS1: Too common for disorder (>1% AF)
- BS2: Homozygotes in gnomAD (for recessive disorders)
- BP1: Missense in truncating-variant gene using pLI/mis_z
- BP3: In-frame indel in repeat region (via UCSC RepeatMasker API)
- BP4: Computational predictor consensus for benign
- BP6: ClinVar benign assertion with review status check
- BP7: Synonymous variant with SpliceAI splice impact check

### ClinVar API-Based (when configured):
- PS1: Same amino acid change as known pathogenic variant
- PM5: Different amino acid change at same position as pathogenic

### Manual Input Required (family/clinical data):
- PS2: De novo (confirmed) - requires trio sequencing
- PS3: Functional studies (damaging) - requires literature review
- PS4: Prevalence in affected - requires case-control data
- PM3: In trans with pathogenic - requires phasing
- PM6: De novo (assumed) - requires family history
- PP1: Cosegregation - requires segregation data
- PP4: Specific phenotype - requires patient phenotype
- BS3: Functional studies (benign) - requires literature review
- BS4: Lack of segregation - requires family data
- BP2: In cis with pathogenic - requires phasing
- BP5: Alternate molecular cause - requires case review

FOR RESEARCH USE ONLY - Not for clinical diagnostic use.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from src.clinvar_api import ClinVarClient


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
    evidence: Optional[str] = None  # Explanation of why criterion was/wasn't met

    # For UI display
    short_name: str = ""
    category: str = ""  # population, computation, intrinsic, clinical, literature


# Define all ACMG criteria with their properties
ACMG_CRITERIA_DEFINITIONS = {
    # Pathogenic - Very Strong
    "PVS1": {
        "name": "Null Known",
        "short_name": "LOF",
        "description": "Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a gene where LOF is a known mechanism of disease",
        "strength": EvidenceStrength.VERY_STRONG,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    # Pathogenic - Strong
    "PS1": {
        "name": "Other Nucleotide",
        "short_name": "Same AA",
        "description": "Same amino acid change as a previously established pathogenic variant regardless of nucleotide change",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    "PS2": {
        "name": "De Novo - Confirmed",
        "short_name": "De Novo",
        "description": "De novo (both maternity and paternity confirmed) in a patient with the disease and no family history",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "clinical",
    },
    "PS3": {
        "name": "Functional Consequence",
        "short_name": "Functional",
        "description": "Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "literature",
    },
    "PS4": {
        "name": "High Prevalence",
        "short_name": "Prevalence",
        "description": "The prevalence of the variant in affected individuals is significantly increased compared to the prevalence in controls",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "population",
    },
    # Pathogenic - Moderate
    "PM1": {
        "name": "Hotspot",
        "short_name": "Hotspot",
        "description": "Located in a mutational hot spot and/or critical and well-established functional domain without benign variation",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    "PM2": {
        "name": "Rare",
        "short_name": "Rare",
        "description": "Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or gnomAD",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "population",
    },
    "PM3": {
        "name": "Trans",
        "short_name": "Trans",
        "description": "For recessive disorders, detected in trans with a pathogenic variant",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "clinical",
    },
    "PM4": {
        "name": "In-frame Non-Repetitive",
        "short_name": "In-frame",
        "description": "Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    "PM5": {
        "name": "Different Missense",
        "short_name": "Diff AA",
        "description": "Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    "PM6": {
        "name": "De Novo - Not Confirmed",
        "short_name": "De Novo?",
        "description": "Assumed de novo, but without confirmation of paternity and maternity",
        "strength": EvidenceStrength.MODERATE,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "clinical",
    },
    # Pathogenic - Supporting
    "PP1": {
        "name": "Segregates - One Family",
        "short_name": "Segregates",
        "description": "Co-segregation with disease in multiple affected family members in a gene definitively known to cause the disease",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "clinical",
    },
    "PP2": {
        "name": "Typically Missense",
        "short_name": "Missense",
        "description": "Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "intrinsic",
    },
    "PP3": {
        "name": "Predicted Damaging",
        "short_name": "Computational",
        "description": "Multiple lines of computational evidence support a deleterious effect on the gene or gene product",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "computation",
    },
    "PP4": {
        "name": "Family History",
        "short_name": "Phenotype",
        "description": "Patient's phenotype or family history is highly specific for a disease with a single genetic etiology",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "clinical",
    },
    "PP5": {
        "name": "Reputable Source - Pathogenic",
        "short_name": "ClinVar P",
        "description": "Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.PATHOGENIC,
        "category": "literature",
    },
    # Benign - Stand-Alone
    "BA1": {
        "name": "Very Common",
        "short_name": "AF ≥5%",
        "description": "Allele frequency is ≥5% in Exome Sequencing Project, 1000 Genomes Project, or gnomAD",
        "strength": EvidenceStrength.STAND_ALONE,
        "evidence_type": EvidenceType.BENIGN,
        "category": "population",
    },
    # Benign - Strong
    "BS1": {
        "name": "Too Common",
        "short_name": "Common",
        "description": "Allele frequency is greater than expected for disorder",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.BENIGN,
        "category": "population",
    },
    "BS2": {
        "name": "Healthy Mutant",
        "short_name": "Healthy",
        "description": "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.BENIGN,
        "category": "population",
    },
    "BS3": {
        "name": "No Functional Consequence",
        "short_name": "Functional OK",
        "description": "Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.BENIGN,
        "category": "literature",
    },
    "BS4": {
        "name": "Segregation - None",
        "short_name": "No Seg",
        "description": "Lack of segregation in affected members of a family",
        "strength": EvidenceStrength.STRONG,
        "evidence_type": EvidenceType.BENIGN,
        "category": "clinical",
    },
    # Benign - Supporting
    "BP1": {
        "name": "Not LOF",
        "short_name": "Missense",
        "description": "Missense variant in a gene for which primarily truncating variants are known to cause disease",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "intrinsic",
    },
    "BP2": {
        "name": "With cis Pathogenic",
        "short_name": "Cis/Trans",
        "description": "Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "clinical",
    },
    "BP3": {
        "name": "In-frame Repetitive",
        "short_name": "Repeat",
        "description": "In-frame deletions/insertions in a repetitive region without a known function",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "intrinsic",
    },
    "BP4": {
        "name": "Computationally Inert",
        "short_name": "Computational",
        "description": "Multiple lines of computational evidence suggest no impact on gene or gene product",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "computation",
    },
    "BP5": {
        "name": "Other Cause",
        "short_name": "Alt Cause",
        "description": "Variant found in a case with an alternate molecular basis for disease",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "clinical",
    },
    "BP6": {
        "name": "Reputable Source - Benign",
        "short_name": "ClinVar B",
        "description": "Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "literature",
    },
    "BP7": {
        "name": "Synonymous",
        "short_name": "Silent",
        "description": "A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved",
        "strength": EvidenceStrength.SUPPORTING,
        "evidence_type": EvidenceType.BENIGN,
        "category": "intrinsic",
    },
}


def create_criterion(
    code: str,
    status: CriterionStatus = CriterionStatus.NOT_EVALUATED,
    evidence: str = None,
    met: bool = None,  # Deprecated, use status instead
) -> ACMGCriterion:
    """Create an ACMGCriterion from a code."""
    if code not in ACMG_CRITERIA_DEFINITIONS:
        raise ValueError(f"Unknown ACMG criterion: {code}")

    defn = ACMG_CRITERIA_DEFINITIONS[code]

    # Handle legacy 'met' parameter
    if met is not None:
        status = CriterionStatus.MET if met else CriterionStatus.NOT_MET
    is_met = status == CriterionStatus.MET

    return ACMGCriterion(
        code=code,
        name=defn["name"],
        short_name=defn["short_name"],
        description=defn["description"],
        strength=defn["strength"],
        evidence_type=defn["evidence_type"],
        category=defn["category"],
        status=status,
        met=is_met,
        evidence=evidence,
    )


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


# Thresholds for automated criteria evaluation
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


def evaluate_population_criteria(
    gnomad_af: Optional[float],
    gnomad_af_popmax: Optional[float] = None,
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
        criteria.append(create_criterion(
            "PM2",
            met=True,
            evidence="Variant absent from gnomAD population database"
        ))
    elif af >= Thresholds.BA1_AF:
        # Very common - BA1
        criteria.append(create_criterion(
            "BA1",
            met=True,
            evidence=f"Allele frequency {af:.4%} ≥ 5% in gnomAD"
        ))
    elif af > Thresholds.BS1_AF:
        # Too common for most Mendelian diseases - BS1
        criteria.append(create_criterion(
            "BS1",
            met=True,
            evidence=f"Allele frequency {af:.4%} > 1% in gnomAD"
        ))
    elif af < Thresholds.PM2_AF:
        # Rare - PM2
        criteria.append(create_criterion(
            "PM2",
            met=True,
            evidence=f"Allele frequency {af:.6%} is extremely rare in gnomAD"
        ))

    return criteria


def evaluate_computational_criteria(
    sift_pred: Optional[str],
    polyphen_pred: Optional[str],
    cadd_phred: Optional[float],
    revel_score: Optional[float],
    alphamissense_pred: Optional[str] = None,
    mutation_taster: Optional[str] = None,
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
        criteria.append(create_criterion(
            "PP3",
            met=True,
            evidence=f"{damaging_count} computational predictors support deleterious effect: " +
                     "; ".join(evidence_details)
        ))

    # BP4: Multiple lines support benign
    elif benign_count >= Thresholds.MIN_PREDICTORS_AGREE:
        criteria.append(create_criterion(
            "BP4",
            met=True,
            evidence=f"{benign_count} computational predictors suggest no impact: " +
                     "; ".join(evidence_details)
        ))

    return criteria


def evaluate_clinvar_criteria(
    clinvar_sig: Optional[str],
    clinvar_review: Optional[str] = None,
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
        if "expert" in review_lower or "practice" in review_lower:
            high_confidence = True
        elif "multiple" in review_lower and "no_conflicts" in review_lower:
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
    variant_type: Optional[str],
    consequence: Optional[str],
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
        "frameshift", "nonsense", "stop_gained", "stop_lost",
        "splice_donor", "splice_acceptor", "start_lost",
        "transcript_ablation"
    ]

    if any(lof in cons_lower for lof in lof_consequences):
        if is_lof_gene:
            criteria.append(create_criterion(
                "PVS1",
                met=True,
                evidence=f"Loss-of-function variant ({consequence}) in gene where LOF is a known disease mechanism"
            ))

    # BP7: Synonymous variants
    if "synonymous" in cons_lower:
        criteria.append(create_criterion(
            "BP7",
            met=True,
            evidence=f"Synonymous variant ({consequence}) - no predicted splice impact"
        ))

    return criteria


def combine_criteria(criteria: list[ACMGCriterion]) -> ACMGScore:
    """
    Combine ACMG criteria according to official scoring rules.

    Returns the final classification based on met criteria.

    Pathogenic:
    1. 1 PVS1 + ≥1 PS
    2. 1 PVS1 + ≥2 PM
    3. 1 PVS1 + 1 PM + 1 PP
    4. 1 PVS1 + ≥2 PP
    5. ≥2 PS
    6. 1 PS + ≥3 PM
    7. 1 PS + 2 PM + ≥2 PP
    8. 1 PS + 1 PM + ≥4 PP

    Likely Pathogenic:
    1. 1 PVS1 + 1 PM
    2. 1 PS + 1-2 PM
    3. 1 PS + ≥2 PP
    4. ≥3 PM
    5. 2 PM + ≥2 PP
    6. 1 PM + ≥4 PP

    Benign:
    1. 1 BA1
    2. ≥2 BS

    Likely Benign:
    1. 1 BS + 1 BP
    2. ≥2 BP
    """
    met_criteria = [c for c in criteria if c.met]

    # Count by strength and type
    pvs = [c for c in met_criteria if c.code.startswith("PVS")]
    ps = [c for c in met_criteria if c.code.startswith("PS")]
    pm = [c for c in met_criteria if c.code.startswith("PM")]
    pp = [c for c in met_criteria if c.code.startswith("PP")]

    ba = [c for c in met_criteria if c.code.startswith("BA")]
    bs = [c for c in met_criteria if c.code.startswith("BS")]
    bp = [c for c in met_criteria if c.code.startswith("BP")]

    n_pvs = len(pvs)
    n_ps = len(ps)
    n_pm = len(pm)
    n_pp = len(pp)
    n_ba = len(ba)
    n_bs = len(bs)
    n_bp = len(bp)

    # Check for contradictory evidence
    has_pathogenic_evidence = n_pvs > 0 or n_ps > 0 or n_pm > 0 or n_pp > 0
    has_benign_evidence = n_ba > 0 or n_bs > 0 or n_bp > 0

    classification = Classification.UNCERTAIN_SIGNIFICANCE
    confidence = 0.5
    rule = ""
    notes = []

    # === Check Benign first (BA1 is stand-alone) ===
    if n_ba >= 1:
        classification = Classification.BENIGN
        confidence = 0.95
        rule = "BA1: Allele frequency ≥5%"
        if has_pathogenic_evidence:
            notes.append("Warning: Contradictory pathogenic evidence present but BA1 overrides")

    elif n_bs >= 2:
        classification = Classification.BENIGN
        confidence = 0.90
        rule = "≥2 BS criteria met"
        if has_pathogenic_evidence:
            notes.append("Warning: Contradictory pathogenic evidence present")

    # === Check Likely Benign ===
    elif n_bs >= 1 and n_bp >= 1:
        classification = Classification.LIKELY_BENIGN
        confidence = 0.80
        rule = "1 BS + 1 BP"
        if has_pathogenic_evidence:
            notes.append("Contradictory pathogenic evidence - review recommended")

    elif n_bp >= 2:
        classification = Classification.LIKELY_BENIGN
        confidence = 0.75
        rule = "≥2 BP criteria met"
        if has_pathogenic_evidence:
            notes.append("Contradictory pathogenic evidence - review recommended")

    # === Check Pathogenic (only if no benign classification) ===
    elif classification == Classification.UNCERTAIN_SIGNIFICANCE:

        # PVS1 combinations
        if n_pvs >= 1:
            if n_ps >= 1:
                classification = Classification.PATHOGENIC
                confidence = 0.95
                rule = "PVS1 + ≥1 PS"
            elif n_pm >= 2:
                classification = Classification.PATHOGENIC
                confidence = 0.90
                rule = "PVS1 + ≥2 PM"
            elif n_pm >= 1 and n_pp >= 1:
                classification = Classification.PATHOGENIC
                confidence = 0.85
                rule = "PVS1 + 1 PM + 1 PP"
            elif n_pp >= 2:
                classification = Classification.PATHOGENIC
                confidence = 0.80
                rule = "PVS1 + ≥2 PP"
            elif n_pm >= 1:
                classification = Classification.LIKELY_PATHOGENIC
                confidence = 0.80
                rule = "PVS1 + 1 PM"

        # PS combinations
        if classification == Classification.UNCERTAIN_SIGNIFICANCE:
            if n_ps >= 2:
                classification = Classification.PATHOGENIC
                confidence = 0.90
                rule = "≥2 PS"
            elif n_ps >= 1:
                if n_pm >= 3:
                    classification = Classification.PATHOGENIC
                    confidence = 0.85
                    rule = "1 PS + ≥3 PM"
                elif n_pm >= 2 and n_pp >= 2:
                    classification = Classification.PATHOGENIC
                    confidence = 0.80
                    rule = "1 PS + 2 PM + ≥2 PP"
                elif n_pm >= 1 and n_pp >= 4:
                    classification = Classification.PATHOGENIC
                    confidence = 0.80
                    rule = "1 PS + 1 PM + ≥4 PP"
                elif n_pm >= 1:
                    classification = Classification.LIKELY_PATHOGENIC
                    confidence = 0.75
                    rule = "1 PS + 1-2 PM"
                elif n_pp >= 2:
                    classification = Classification.LIKELY_PATHOGENIC
                    confidence = 0.70
                    rule = "1 PS + ≥2 PP"

        # PM-only combinations
        if classification == Classification.UNCERTAIN_SIGNIFICANCE:
            if n_pm >= 3:
                classification = Classification.LIKELY_PATHOGENIC
                confidence = 0.75
                rule = "≥3 PM"
            elif n_pm >= 2 and n_pp >= 2:
                classification = Classification.LIKELY_PATHOGENIC
                confidence = 0.70
                rule = "2 PM + ≥2 PP"
            elif n_pm >= 1 and n_pp >= 4:
                classification = Classification.LIKELY_PATHOGENIC
                confidence = 0.70
                rule = "1 PM + ≥4 PP"

    # Check for contradictions
    if has_pathogenic_evidence and has_benign_evidence:
        if classification == Classification.UNCERTAIN_SIGNIFICANCE:
            notes.append("Contradictory pathogenic and benign evidence - classified as VUS")

    return ACMGScore(
        classification=classification,
        criteria_met=met_criteria,
        criteria_not_met=[c for c in criteria if not c.met],
        confidence=confidence,
        rule_applied=rule,
        notes=notes,
    )


def score_variant(
    gnomad_af: Optional[float] = None,
    gnomad_af_popmax: Optional[float] = None,
    sift_pred: Optional[str] = None,
    polyphen_pred: Optional[str] = None,
    cadd_phred: Optional[float] = None,
    revel_score: Optional[float] = None,
    alphamissense_pred: Optional[str] = None,
    mutation_taster: Optional[str] = None,
    clinvar_sig: Optional[str] = None,
    clinvar_review: Optional[str] = None,
    consequence: Optional[str] = None,
    variant_type: Optional[str] = None,
    is_lof_gene: bool = True,
) -> ACMGScore:
    """
    Score a variant according to ACMG/AMP guidelines.

    This is the main entry point for variant classification.

    Args:
        gnomad_af: gnomAD allele frequency (global)
        gnomad_af_popmax: gnomAD maximum population allele frequency
        sift_pred: SIFT prediction (D=deleterious, T=tolerated)
        polyphen_pred: PolyPhen-2 prediction (D=damaging, P=possibly, B=benign)
        cadd_phred: CADD phred-scaled score
        revel_score: REVEL score (0-1)
        alphamissense_pred: AlphaMissense prediction
        mutation_taster: MutationTaster prediction
        clinvar_sig: ClinVar clinical significance
        clinvar_review: ClinVar review status
        consequence: Variant consequence (e.g., "missense_variant", "frameshift")
        variant_type: Variant type (SNV, indel, etc.)
        is_lof_gene: Whether LOF is a known disease mechanism for this gene

    Returns:
        ACMGScore with classification and criteria details
    """
    all_criteria = []

    # Population criteria
    all_criteria.extend(evaluate_population_criteria(gnomad_af, gnomad_af_popmax))

    # Computational criteria
    all_criteria.extend(evaluate_computational_criteria(
        sift_pred, polyphen_pred, cadd_phred, revel_score,
        alphamissense_pred, mutation_taster
    ))

    # ClinVar criteria
    all_criteria.extend(evaluate_clinvar_criteria(clinvar_sig, clinvar_review))

    # Variant type criteria
    all_criteria.extend(evaluate_variant_type_criteria(
        variant_type, consequence, is_lof_gene
    ))

    # Combine and classify
    return combine_criteria(all_criteria)


def score_variant_from_metadata(metadata: dict) -> ACMGScore:
    """
    Score a variant using metadata dictionary from the vector store.

    This is a convenience function that extracts the relevant fields
    from the metadata and calls score_variant().
    """
    return score_variant(
        gnomad_af=metadata.get("gnomad_af"),
        gnomad_af_popmax=metadata.get("gnomad_af_popmax"),
        sift_pred=metadata.get("sift_pred"),
        polyphen_pred=metadata.get("polyphen_pred"),
        cadd_phred=metadata.get("cadd_phred"),
        revel_score=metadata.get("revel_score"),
        alphamissense_pred=metadata.get("alphamissense_pred"),
        mutation_taster=metadata.get("mutationtaster_pred") or metadata.get("mutation_taster"),
        clinvar_sig=metadata.get("clinvar_sig"),
        clinvar_review=metadata.get("clinvar_review"),
        consequence=metadata.get("consequence"),
        variant_type=metadata.get("variant_type"),
    )


def evaluate_all_criteria(
    metadata: dict, clinvar_client: Optional["ClinVarClient"] = None
) -> list[ACMGCriterion]:
    """
    Evaluate ALL 28 ACMG criteria with explanations for each.

    Args:
        metadata: Variant metadata dictionary
        clinvar_client: Optional ClinVarClient for PS1/PM5 evaluation

    Returns a list of all criteria with their status (MET, NOT_MET, NOT_EVALUATED)
    and evidence explaining the evaluation.

    This function is designed to provide full transparency into the ACMG
    classification process, showing exactly why each criterion was or wasn't applied.
    """
    all_criteria = []

    # Extract values
    gnomad_af = metadata.get("gnomad_af")
    sift_pred = metadata.get("sift_pred")
    polyphen_pred = metadata.get("polyphen_pred")
    cadd_phred = metadata.get("cadd_phred")
    revel_score = metadata.get("revel_score")
    alphamissense_pred = metadata.get("alphamissense_pred")
    mutation_taster = metadata.get("mutationtaster_pred") or metadata.get("mutation_taster")
    clinvar_sig = metadata.get("clinvar_sig")
    clinvar_review = metadata.get("clinvar_review")
    consequence = metadata.get("consequence", "")
    cons_lower = consequence.lower() if consequence else ""

    # ========== PATHOGENIC CRITERIA ==========

    # PVS1: Null variant in LOF gene (enhanced with gene constraint)
    lof_consequences = ["frameshift", "nonsense", "stop_gained", "stop_lost",
                        "splice_donor", "splice_acceptor", "start_lost"]
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

        all_criteria.append(create_criterion(
            "PVS1", CriterionStatus.MET,
            f"Loss-of-function variant ({consequence}) in gene where LOF is a known disease mechanism.{constraint_info}"
        ))
    elif consequence:
        all_criteria.append(create_criterion(
            "PVS1", CriterionStatus.NOT_MET,
            f"Not a loss-of-function variant (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "PVS1", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

    # PS1: Same amino acid change as known pathogenic
    # (Evaluated together with PM5 below after PM4)

    # PS2: De novo confirmed
    all_criteria.append(create_criterion(
        "PS2", CriterionStatus.NOT_EVALUATED,
        "Requires family data to confirm de novo status (not available from dbNSFP)"
    ))

    # PS3: Functional studies
    all_criteria.append(create_criterion(
        "PS3", CriterionStatus.NOT_EVALUATED,
        "Requires review of published functional studies (not automated)"
    ))

    # PS4: High prevalence in affected
    all_criteria.append(create_criterion(
        "PS4", CriterionStatus.NOT_EVALUATED,
        "Requires case-control study data (not available from dbNSFP)"
    ))

    # PM1: Mutational hotspot / functional domain
    interpro_domain = metadata.get("interpro_domain", "")
    if interpro_domain and interpro_domain not in ("", ".", "N/A"):
        all_criteria.append(create_criterion(
            "PM1", CriterionStatus.MET,
            f"Located in functional domain: {interpro_domain}"
        ))
    else:
        all_criteria.append(create_criterion(
            "PM1", CriterionStatus.NOT_EVALUATED,
            "No domain annotation available (InterPro domain not found)"
        ))

    # PM2: Absent/rare in population
    if gnomad_af is None:
        all_criteria.append(create_criterion(
            "PM2", CriterionStatus.MET,
            "Variant absent from gnomAD population database"
        ))
    elif gnomad_af < Thresholds.PM2_AF:
        all_criteria.append(create_criterion(
            "PM2", CriterionStatus.MET,
            f"Allele frequency {gnomad_af:.6%} is extremely rare (< 0.01%) in gnomAD"
        ))
    else:
        all_criteria.append(create_criterion(
            "PM2", CriterionStatus.NOT_MET,
            f"Allele frequency {gnomad_af:.4%} exceeds rare threshold (0.01%)"
        ))

    # PM3: In trans with known pathogenic
    all_criteria.append(create_criterion(
        "PM3", CriterionStatus.NOT_EVALUATED,
        "Requires phasing data to determine cis/trans status (not available)"
    ))

    # PM4: Protein length change (in-frame indel, stop-loss)
    if "inframe" in cons_lower or "stop_lost" in cons_lower:
        all_criteria.append(create_criterion(
            "PM4", CriterionStatus.MET,
            f"In-frame protein length change ({consequence})"
        ))
    elif consequence:
        all_criteria.append(create_criterion(
            "PM4", CriterionStatus.NOT_MET,
            f"Not an in-frame protein length change (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "PM4", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

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
                    all_criteria.append(create_criterion(
                        "PS1", CriterionStatus.MET, evidence
                    ))
                else:
                    all_criteria.append(create_criterion(
                        "PS1", CriterionStatus.NOT_MET,
                        "No pathogenic variant with same amino acid change found in ClinVar"
                    ))

                if pm5_met:
                    all_criteria.append(create_criterion(
                        "PM5", CriterionStatus.MET, evidence
                    ))
                else:
                    all_criteria.append(create_criterion(
                        "PM5", CriterionStatus.NOT_MET,
                        "No pathogenic variant at same position with different AA change in ClinVar"
                    ))
            else:
                all_criteria.append(create_criterion(
                    "PS1", CriterionStatus.NOT_EVALUATED,
                    "Amino acid change data incomplete (missing position or ref/alt AA)"
                ))
                all_criteria.append(create_criterion(
                    "PM5", CriterionStatus.NOT_EVALUATED,
                    "Amino acid change data incomplete (missing position or ref/alt AA)"
                ))
        except Exception as e:
            all_criteria.append(create_criterion(
                "PS1", CriterionStatus.NOT_EVALUATED,
                f"ClinVar lookup failed: {e}"
            ))
            all_criteria.append(create_criterion(
                "PM5", CriterionStatus.NOT_EVALUATED,
                f"ClinVar lookup failed: {e}"
            ))
    else:
        all_criteria.append(create_criterion(
            "PS1", CriterionStatus.NOT_EVALUATED,
            "Requires ClinVar API client for automated evaluation (not configured)"
        ))
        all_criteria.append(create_criterion(
            "PM5", CriterionStatus.NOT_EVALUATED,
            "Requires ClinVar API client for automated evaluation (not configured)"
        ))

    # PM6: Assumed de novo (without confirmation)
    all_criteria.append(create_criterion(
        "PM6", CriterionStatus.NOT_EVALUATED,
        "Requires family history/segregation data (not available from dbNSFP)"
    ))

    # PP1: Cosegregation with disease
    all_criteria.append(create_criterion(
        "PP1", CriterionStatus.NOT_EVALUATED,
        "Requires family segregation data (not available from dbNSFP)"
    ))

    # PP2: Missense in gene with low benign missense rate
    # Using gnomad_mis_oe (missense observed/expected ratio) - lower is more constrained
    gnomad_mis_oe = metadata.get("gnomad_mis_oe")
    if "missense" in cons_lower:
        if gnomad_mis_oe is not None:
            # Low missense o/e ratio (<0.6) indicates gene intolerant to missense
            if gnomad_mis_oe < 0.6:
                all_criteria.append(create_criterion(
                    "PP2", CriterionStatus.MET,
                    f"Missense variant in gene with high missense constraint (gnomAD mis o/e = {gnomad_mis_oe:.2f} < 0.6)"
                ))
            else:
                all_criteria.append(create_criterion(
                    "PP2", CriterionStatus.NOT_MET,
                    f"Gene not highly missense-constrained (gnomAD mis o/e = {gnomad_mis_oe:.2f} ≥ 0.6)"
                ))
        else:
            all_criteria.append(create_criterion(
                "PP2", CriterionStatus.NOT_EVALUATED,
                "Missense variant - gene constraint score not available"
            ))
    elif consequence:
        all_criteria.append(create_criterion(
            "PP2", CriterionStatus.NOT_MET,
            f"Not a missense variant (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "PP2", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

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

    if predictor_details:
        if damaging_count >= Thresholds.MIN_PREDICTORS_AGREE:
            all_criteria.append(create_criterion(
                "PP3", CriterionStatus.MET,
                f"{damaging_count} of {len(predictor_details)} predictors support damaging: {'; '.join(predictor_details)}"
            ))
        else:
            all_criteria.append(create_criterion(
                "PP3", CriterionStatus.NOT_MET,
                f"Only {damaging_count} of {len(predictor_details)} predictors support damaging (need ≥{Thresholds.MIN_PREDICTORS_AGREE}): {'; '.join(predictor_details)}"
            ))
    else:
        all_criteria.append(create_criterion(
            "PP3", CriterionStatus.NOT_EVALUATED,
            "No computational predictions available"
        ))

    # PP4: Patient phenotype highly specific
    all_criteria.append(create_criterion(
        "PP4", CriterionStatus.NOT_EVALUATED,
        "Requires patient phenotype data (not available from dbNSFP)"
    ))

    # PP5: Reputable source reports pathogenic
    if clinvar_sig:
        sig_lower = clinvar_sig.lower()
        has_high_review = clinvar_review and ("expert" in clinvar_review.lower() or "multiple" in clinvar_review.lower())
        if "pathogenic" in sig_lower and "benign" not in sig_lower:
            if has_high_review:
                all_criteria.append(create_criterion(
                    "PP5", CriterionStatus.MET,
                    f"ClinVar reports {clinvar_sig} with review status: {clinvar_review}"
                ))
            else:
                all_criteria.append(create_criterion(
                    "PP5", CriterionStatus.NOT_MET,
                    f"ClinVar reports {clinvar_sig} but review status ({clinvar_review}) not sufficient for PP5"
                ))
        else:
            all_criteria.append(create_criterion(
                "PP5", CriterionStatus.NOT_MET,
                f"ClinVar significance ({clinvar_sig}) is not pathogenic"
            ))
    else:
        all_criteria.append(create_criterion(
            "PP5", CriterionStatus.NOT_EVALUATED,
            "No ClinVar annotation available"
        ))

    # ========== BENIGN CRITERIA ==========

    # BA1: Very high population frequency (≥5%)
    if gnomad_af is not None and gnomad_af >= Thresholds.BA1_AF:
        all_criteria.append(create_criterion(
            "BA1", CriterionStatus.MET,
            f"Allele frequency {gnomad_af:.4%} is ≥5% in gnomAD (stand-alone benign)"
        ))
    elif gnomad_af is not None:
        all_criteria.append(create_criterion(
            "BA1", CriterionStatus.NOT_MET,
            f"Allele frequency {gnomad_af:.4%} is below 5% threshold"
        ))
    else:
        all_criteria.append(create_criterion(
            "BA1", CriterionStatus.NOT_MET,
            "Variant not found in gnomAD (cannot be common)"
        ))

    # BS1: Too common for disorder (>1%)
    if gnomad_af is not None and gnomad_af > Thresholds.BS1_AF:
        all_criteria.append(create_criterion(
            "BS1", CriterionStatus.MET,
            f"Allele frequency {gnomad_af:.4%} exceeds 1% (too common for most Mendelian disorders)"
        ))
    elif gnomad_af is not None:
        all_criteria.append(create_criterion(
            "BS1", CriterionStatus.NOT_MET,
            f"Allele frequency {gnomad_af:.4%} is below 1% threshold"
        ))
    else:
        all_criteria.append(create_criterion(
            "BS1", CriterionStatus.NOT_MET,
            "Variant not found in gnomAD"
        ))

    # BS2: Observed in healthy adult (for penetrant disorder)
    # Partially evaluable using gnomAD homozygote count (for recessive disorders)
    gnomad_hom = metadata.get("gnomad_hom")
    if gnomad_hom is not None and gnomad_hom > 0:
        all_criteria.append(create_criterion(
            "BS2", CriterionStatus.MET,
            f"Observed as homozygous in {int(gnomad_hom)} healthy individual(s) in gnomAD (supports benign for recessive disorders)"
        ))
    else:
        all_criteria.append(create_criterion(
            "BS2", CriterionStatus.NOT_EVALUATED,
            "No homozygotes observed in gnomAD. For dominant disorders, requires individual-level phenotype data."
        ))

    # BS3: Functional studies show no damaging effect
    all_criteria.append(create_criterion(
        "BS3", CriterionStatus.NOT_EVALUATED,
        "Requires review of published functional studies (not automated)"
    ))

    # BS4: Lack of segregation in affected family members
    all_criteria.append(create_criterion(
        "BS4", CriterionStatus.NOT_EVALUATED,
        "Requires family segregation data (not available from dbNSFP)"
    ))

    # BP1: Missense in gene where only truncating causes disease
    gnomad_pli = metadata.get("gnomad_pli")
    if "missense" in cons_lower:
        if gnomad_pli is not None:
            # High pLI (>0.9) + high missense o/e (>0.8) suggests truncating is primary mechanism
            # (LOF intolerant but missense tolerated)
            if gnomad_pli >= 0.9 and (gnomad_mis_oe is None or gnomad_mis_oe > 0.8):
                mis_oe_str = f"{gnomad_mis_oe:.2f}" if gnomad_mis_oe else "N/A"
                all_criteria.append(create_criterion(
                    "BP1", CriterionStatus.MET,
                    f"Missense in gene where truncating variants are primary disease mechanism (pLI = {gnomad_pli:.2f}, mis o/e = {mis_oe_str})"
                ))
            else:
                all_criteria.append(create_criterion(
                    "BP1", CriterionStatus.NOT_MET,
                    f"Gene constraint does not suggest truncating-only mechanism (pLI = {gnomad_pli:.2f})"
                ))
        else:
            all_criteria.append(create_criterion(
                "BP1", CriterionStatus.NOT_EVALUATED,
                "Gene constraint scores not available"
            ))
    elif consequence:
        all_criteria.append(create_criterion(
            "BP1", CriterionStatus.NOT_MET,
            f"Not a missense variant (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "BP1", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

    # BP2: Observed in cis/trans with pathogenic
    all_criteria.append(create_criterion(
        "BP2", CriterionStatus.NOT_EVALUATED,
        "Requires phasing data to determine cis/trans status (not available)"
    ))

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
                    all_criteria.append(create_criterion(
                        "BP3", CriterionStatus.MET,
                        f"In-frame indel in {repeat_class}/{repeat_family} repetitive region"
                    ))
                else:
                    all_criteria.append(create_criterion(
                        "BP3", CriterionStatus.NOT_MET,
                        "In-frame indel but not in a repetitive region"
                    ))
            else:
                all_criteria.append(create_criterion(
                    "BP3", CriterionStatus.NOT_EVALUATED,
                    "Position data not available for repeat region lookup"
                ))
        except Exception as e:
            all_criteria.append(create_criterion(
                "BP3", CriterionStatus.NOT_EVALUATED,
                f"Repeat region lookup failed: {e}"
            ))
    elif consequence:
        all_criteria.append(create_criterion(
            "BP3", CriterionStatus.NOT_MET,
            f"Not an in-frame variant (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "BP3", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

    # BP4: Computational evidence supports benign
    if predictor_details:
        if benign_count >= Thresholds.MIN_PREDICTORS_AGREE:
            all_criteria.append(create_criterion(
                "BP4", CriterionStatus.MET,
                f"{benign_count} of {len(predictor_details)} predictors support benign: {'; '.join(predictor_details)}"
            ))
        else:
            all_criteria.append(create_criterion(
                "BP4", CriterionStatus.NOT_MET,
                f"Only {benign_count} of {len(predictor_details)} predictors support benign (need ≥{Thresholds.MIN_PREDICTORS_AGREE}): {'; '.join(predictor_details)}"
            ))
    else:
        all_criteria.append(create_criterion(
            "BP4", CriterionStatus.NOT_EVALUATED,
            "No computational predictions available"
        ))

    # BP5: Variant found in case with alternate cause
    all_criteria.append(create_criterion(
        "BP5", CriterionStatus.NOT_EVALUATED,
        "Requires case-level data showing alternate molecular cause (not available)"
    ))

    # BP6: Reputable source reports benign
    if clinvar_sig:
        sig_lower = clinvar_sig.lower()
        has_high_review = clinvar_review and ("expert" in clinvar_review.lower() or "multiple" in clinvar_review.lower())
        if "benign" in sig_lower and "pathogenic" not in sig_lower:
            if has_high_review:
                all_criteria.append(create_criterion(
                    "BP6", CriterionStatus.MET,
                    f"ClinVar reports {clinvar_sig} with review status: {clinvar_review}"
                ))
            else:
                all_criteria.append(create_criterion(
                    "BP6", CriterionStatus.NOT_MET,
                    f"ClinVar reports {clinvar_sig} but review status ({clinvar_review}) not sufficient for BP6"
                ))
        else:
            all_criteria.append(create_criterion(
                "BP6", CriterionStatus.NOT_MET,
                f"ClinVar significance ({clinvar_sig}) is not benign"
            ))
    else:
        all_criteria.append(create_criterion(
            "BP6", CriterionStatus.NOT_EVALUATED,
            "No ClinVar annotation available"
        ))

    # BP7: Synonymous with no splice impact (enhanced with SpliceAI)
    if "synonymous" in cons_lower:
        # Get SpliceAI scores and compute max delta
        spliceai_ag = metadata.get("spliceai_ag")
        spliceai_al = metadata.get("spliceai_al")
        spliceai_dg = metadata.get("spliceai_dg")
        spliceai_dl = metadata.get("spliceai_dl")

        spliceai_scores = [s for s in [spliceai_ag, spliceai_al, spliceai_dg, spliceai_dl] if s is not None]
        spliceai_max = max(spliceai_scores) if spliceai_scores else None

        if spliceai_max is not None:
            if spliceai_max < 0.2:
                all_criteria.append(create_criterion(
                    "BP7", CriterionStatus.MET,
                    f"Synonymous variant ({consequence}) with no predicted splice impact (SpliceAI max={spliceai_max:.2f} < 0.2)"
                ))
            else:
                # Synonymous but potential splice impact - DO NOT apply BP7
                all_criteria.append(create_criterion(
                    "BP7", CriterionStatus.NOT_MET,
                    f"Synonymous variant but potential splice impact (SpliceAI max={spliceai_max:.2f} ≥ 0.2) - BP7 not applied"
                ))
        else:
            # No SpliceAI data available - apply with caution
            all_criteria.append(create_criterion(
                "BP7", CriterionStatus.MET,
                f"Synonymous variant ({consequence}) - no amino acid change (SpliceAI not available)"
            ))
    elif consequence:
        all_criteria.append(create_criterion(
            "BP7", CriterionStatus.NOT_MET,
            f"Not a synonymous variant (consequence: {consequence})"
        ))
    else:
        all_criteria.append(create_criterion(
            "BP7", CriterionStatus.NOT_EVALUATED,
            "Consequence type not available"
        ))

    return all_criteria

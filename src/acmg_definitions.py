"""Definitions and construction of the 28 ACMG evidence criteria."""

from src.acmg_types import ACMGCriterion, CriterionStatus, EvidenceStrength, EvidenceType

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
    evidence: str | None = None,
    met: bool | None = None,  # Deprecated, use status instead
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

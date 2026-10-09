"""Request and response schemas for the classification API."""

from pydantic import BaseModel, ConfigDict, Field

RUO_DISCLAIMER = (
    "FOR RESEARCH USE ONLY. Not for use in diagnostic procedures. "
    "Classifications require independent expert review before any clinical decision-making."
)


class VariantRequest(BaseModel):
    chr: str = Field(..., description="Chromosome (1-22, X, Y)")
    pos: int = Field(..., description="Position (1-based)")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    genome_build: str = Field(default="GRCh37", description="Genome build (GRCh37 or GRCh38)")

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "chr": "17",
                "pos": 41197801,
                "ref": "T",
                "alt": "A",
                "genome_build": "GRCh37",
            }
        }
    )


class ACMGCriteria(BaseModel):
    code: str = Field(..., description="ACMG criteria code (e.g., PVS1, PS1, PM2)")
    description: str = Field(..., description="Evidence description")


class ACMGCriteriaDetail(BaseModel):
    """Detailed ACMG criterion with full information from rule-based scoring."""

    code: str = Field(..., description="ACMG criteria code (e.g., PVS1, PS1, PM2)")
    name: str = Field(..., description="Criterion name")
    short_name: str = Field(..., description="Short display name")
    description: str = Field(..., description="Full description of the criterion")
    strength: str = Field(
        ..., description="Evidence strength (very_strong, strong, moderate, supporting)"
    )
    evidence_type: str = Field(..., description="Evidence type (pathogenic or benign)")
    category: str = Field(
        ..., description="Category (population, computation, intrinsic, clinical, literature)"
    )
    status: str = Field(
        default="met", description="Evaluation status: met, not_met, or not_evaluated"
    )
    evidence: str | None = Field(None, description="Specific evidence/explanation for this variant")


class ClinVarAnnotation(BaseModel):
    """ClinVar annotation data."""

    id: str | None = Field(None, description="ClinVar ID (e.g., RCV001234567)")
    significance: str | None = Field(None, description="Clinical significance")
    review_status: str | None = Field(
        None,
        description="Review status (e.g., 'criteria_provided,_multiple_submitters,_no_conflicts')",
    )
    trait: str | None = Field(None, description="Associated disease/phenotype")


class ClassificationResponse(BaseModel):
    variant_id: str = Field(..., description="Variant identifier (chr_pos_ref_alt)")
    gene: str | None = Field(None, description="Gene symbol")
    acmg_classification: str = Field(..., description="ACMG 5-tier classification")
    criteria: list[ACMGCriteria] = Field(
        default_factory=list, description="Applied ACMG criteria (legacy, from LLM)"
    )
    criteria_met: list[ACMGCriteriaDetail] = Field(
        default_factory=list, description="Met ACMG criteria with full details"
    )
    all_criteria: list[ACMGCriteriaDetail] = Field(
        default_factory=list, description="All 28 ACMG criteria with evaluation status"
    )
    rule_applied: str = Field(default="", description="ACMG combining rule that was applied")
    scoring_method: str = Field(
        default="rule_based", description="Classification method: 'rule_based' or 'llm_fallback'"
    )
    confidence: float = Field(..., description="Classification confidence (0-1)")
    interpretation: str = Field(..., description="Full interpretation text")
    scores: dict = Field(default_factory=dict, description="Pathogenicity scores")
    clinvar: ClinVarAnnotation | None = Field(None, description="ClinVar annotation if available")
    evidence_links: dict = Field(
        default_factory=dict, description="Links to external databases for evidence verification"
    )
    genome_build: str = Field(default="GRCh37", description="Genome build used for coordinates")
    disclaimer: str = Field(default=RUO_DISCLAIMER, description="Research use disclaimer")
    scoring_notes: list[str] = Field(
        default_factory=list, description="Additional notes from scoring"
    )


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
    database_loaded: bool
    variant_count: int
    genome_build: str = "GRCh37"
    dbnsfp_version: str = "5.3.1a"
    disclaimer: str = RUO_DISCLAIMER

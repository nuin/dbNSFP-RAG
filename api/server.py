#!/usr/bin/env python
"""
ACMG Variant Classification API

Exposes a REST API for variant classification using a fine-tuned LLM.

Usage:
    uvicorn api.server:app --host 0.0.0.0 --port 8000

    # Or with reload for development
    uvicorn api.server:app --reload

Endpoints:
    GET  /classify?chr=17&pos=41197801&ref=G&alt=A
    POST /classify  (JSON body: {"chr": "17", "pos": 41197801, "ref": "G", "alt": "A"})
    GET  /health
"""

import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

# Add project root to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import validation module
from src.validation import validate_variant, validate_reference_against_database

# Import ACMG scoring module
from src.acmg_scoring import score_variant_from_metadata, evaluate_all_criteria

# Import ClinVar API client (optional, for PS1/PM5 evaluation)
try:
    from src.clinvar_api import ClinVarClient
    CLINVAR_CLIENT = ClinVarClient(email="dbnsfp-api@localhost")
except ImportError:
    CLINVAR_CLIENT = None

# Import evidence links generator
from src.evidence_links import generate_evidence_links

# Research Use Only disclaimer
RUO_DISCLAIMER = (
    "FOR RESEARCH USE ONLY. Not for use in diagnostic procedures. "
    "Classifications require independent expert review before any clinical decision-making."
)

app = FastAPI(
    title="ACMG Variant Classifier",
    description=f"Classify genetic variants according to ACMG/AMP guidelines.\n\n**{RUO_DISCLAIMER}**",
    version="1.0.0",
)

# CORS for web access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve static files (web interface)
static_dir = Path(__file__).parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")


@app.get("/", include_in_schema=False)
async def root():
    """Serve the web interface."""
    index_file = static_dir / "index.html"
    if index_file.exists():
        return FileResponse(
            str(index_file),
            headers={"Cache-Control": "no-cache, no-store, must-revalidate"}
        )
    return {"message": "ACMG Variant Classifier API", "docs": "/docs"}


@app.get("/lookup.html", include_in_schema=False)
async def lookup_page():
    """Serve the dbNSFP lookup page."""
    lookup_file = static_dir / "lookup.html"
    if lookup_file.exists():
        return FileResponse(str(lookup_file))
    return {"message": "Lookup page not found"}


# Supported genome builds
SUPPORTED_BUILDS = {"GRCh37", "GRCh38"}


# Request/Response models
class VariantRequest(BaseModel):
    chr: str = Field(..., description="Chromosome (1-22, X, Y)")
    pos: int = Field(..., description="Position (1-based)")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    genome_build: str = Field(default="GRCh37", description="Genome build (GRCh37 or GRCh38)")

    class Config:
        json_schema_extra = {
            "example": {
                "chr": "17",
                "pos": 41197801,
                "ref": "T",
                "alt": "A",
                "genome_build": "GRCh37"
            }
        }


class ACMGCriteria(BaseModel):
    code: str = Field(..., description="ACMG criteria code (e.g., PVS1, PS1, PM2)")
    description: str = Field(..., description="Evidence description")


class ACMGCriteriaDetail(BaseModel):
    """Detailed ACMG criterion with full information from rule-based scoring."""
    code: str = Field(..., description="ACMG criteria code (e.g., PVS1, PS1, PM2)")
    name: str = Field(..., description="Criterion name")
    short_name: str = Field(..., description="Short display name")
    description: str = Field(..., description="Full description of the criterion")
    strength: str = Field(..., description="Evidence strength (very_strong, strong, moderate, supporting)")
    evidence_type: str = Field(..., description="Evidence type (pathogenic or benign)")
    category: str = Field(..., description="Category (population, computation, intrinsic, clinical, literature)")
    status: str = Field(default="met", description="Evaluation status: met, not_met, or not_evaluated")
    evidence: Optional[str] = Field(None, description="Specific evidence/explanation for this variant")


class ClinVarAnnotation(BaseModel):
    """ClinVar annotation data."""
    id: Optional[str] = Field(None, description="ClinVar ID (e.g., RCV001234567)")
    significance: Optional[str] = Field(None, description="Clinical significance")
    review_status: Optional[str] = Field(None, description="Review status (e.g., 'criteria_provided,_multiple_submitters,_no_conflicts')")
    trait: Optional[str] = Field(None, description="Associated disease/phenotype")


class ClassificationResponse(BaseModel):
    variant_id: str = Field(..., description="Variant identifier (chr_pos_ref_alt)")
    gene: Optional[str] = Field(None, description="Gene symbol")
    acmg_classification: str = Field(..., description="ACMG 5-tier classification")
    criteria: list[ACMGCriteria] = Field(default_factory=list, description="Applied ACMG criteria (legacy, from LLM)")
    criteria_met: list[ACMGCriteriaDetail] = Field(default_factory=list, description="Met ACMG criteria with full details")
    all_criteria: list[ACMGCriteriaDetail] = Field(default_factory=list, description="All 28 ACMG criteria with evaluation status")
    rule_applied: str = Field(default="", description="ACMG combining rule that was applied")
    scoring_method: str = Field(default="rule_based", description="Classification method: 'rule_based' or 'llm_fallback'")
    confidence: float = Field(..., description="Classification confidence (0-1)")
    interpretation: str = Field(..., description="Full interpretation text")
    scores: dict = Field(default_factory=dict, description="Pathogenicity scores")
    clinvar: Optional[ClinVarAnnotation] = Field(None, description="ClinVar annotation if available")
    evidence_links: dict = Field(default_factory=dict, description="Links to external databases for evidence verification")
    genome_build: str = Field(default="GRCh37", description="Genome build used for coordinates")
    disclaimer: str = Field(default=RUO_DISCLAIMER, description="Research use disclaimer")
    scoring_notes: list[str] = Field(default_factory=list, description="Additional notes from scoring")


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
    database_loaded: bool
    variant_count: int
    genome_build: str = "GRCh37"
    dbnsfp_version: str = "5.3.1a"
    disclaimer: str = RUO_DISCLAIMER


# Global model and database instances
_model = None
_databases = {}  # Cache databases by genome build


def get_model():
    """Load fine-tuned model (lazy loading).

    Tries backends in order:
    1. llama-cpp-python (GGUF model for CPU/GPU inference)
    2. mlx-lm (Apple Silicon)
    3. transformers (NVIDIA GPU or CPU)
    4. Ollama (external service fallback)
    """
    global _model
    if _model is None:
        model_path = os.environ.get("ACMG_MODEL_PATH", "models/acmg-classifier/model")
        gguf_model_path = os.environ.get("LLM_MODEL_PATH", "/app/models/model.gguf")
        use_bundled_llm = os.environ.get("USE_BUNDLED_LLM", "false").lower() == "true"

        # Try llama-cpp-python first (bundled GGUF model for Docker deployment)
        if use_bundled_llm or Path(gguf_model_path).exists():
            try:
                from llama_cpp import Llama
                n_ctx = int(os.environ.get("LLM_N_CTX", "4096"))
                n_threads = int(os.environ.get("LLM_N_THREADS", "4"))
                n_gpu_layers = int(os.environ.get("LLM_N_GPU_LAYERS", "0"))

                print(f"Loading GGUF model from {gguf_model_path}...")
                llm = Llama(
                    model_path=gguf_model_path,
                    n_ctx=n_ctx,
                    n_threads=n_threads,
                    n_gpu_layers=n_gpu_layers,
                    verbose=False,
                )
                _model = {"type": "llama_cpp", "model": llm}
                print(f"Loaded llama.cpp model (n_ctx={n_ctx}, threads={n_threads})")
                return _model
            except ImportError:
                print("llama-cpp-python not installed, trying other backends...")
            except Exception as e:
                print(f"llama.cpp model load failed: {e}")

        # Try mlx-lm (Apple Silicon)
        try:
            from mlx_lm import load, generate  # noqa: F401 - generate is used in generate_classification()
            model, tokenizer = load(model_path)
            _model = {"type": "mlx", "model": model, "tokenizer": tokenizer}
            print(f"Loaded MLX model from {model_path}")
        except ImportError:
            pass
        except Exception as e:
            print(f"MLX model load failed: {e}")

        # Fall back to transformers
        if _model is None:
            try:
                from transformers import AutoModelForCausalLM, AutoTokenizer
                import torch

                device = "cuda" if torch.cuda.is_available() else "cpu"
                tokenizer = AutoTokenizer.from_pretrained(model_path)
                model = AutoModelForCausalLM.from_pretrained(
                    model_path,
                    torch_dtype=torch.float16 if device == "cuda" else torch.float32,
                    device_map="auto" if device == "cuda" else None,
                )
                _model = {"type": "transformers", "model": model, "tokenizer": tokenizer, "device": device}
                print(f"Loaded Transformers model from {model_path}")
            except Exception as e:
                print(f"Could not load fine-tuned model: {e}")
                print("Falling back to Ollama...")
                _model = {"type": "ollama", "model": "llama3.2:3b"}

    return _model


def get_database(genome_build: str = "GRCh37"):
    """Load variant database for specified genome build (lazy loading).

    Uses SQLite backend (VariantDatabase) if the configured path is a .db file,
    otherwise falls back to FAISS (VariantVectorStore) for backwards compatibility.
    """
    global _databases

    if genome_build not in SUPPORTED_BUILDS:
        raise ValueError(f"Unsupported genome build: {genome_build}. Must be one of {SUPPORTED_BUILDS}")

    if genome_build not in _databases:
        from src.config import SQLITE_DB_PATH, VECTORDB_GRCH37_NGSGENES, VECTORDB_GRCH38_NGSGENES

        if genome_build == "GRCh37":
            db_path = os.environ.get("ACMG_DB_PATH", str(SQLITE_DB_PATH))
        else:
            db_path = os.environ.get("ACMG_DB_PATH_GRCH38", str(VECTORDB_GRCH38_NGSGENES))

        db_path = Path(db_path)

        if db_path.suffix == ".db" and db_path.exists():
            from src.variantdb import VariantDatabase
            _databases[genome_build] = VariantDatabase(db_path=db_path)
            print(f"Loaded {genome_build} SQLite database with {_databases[genome_build].count()} variants")
        else:
            # Fall back to FAISS for backwards compatibility
            from src.vectorstore import VariantVectorStore
            # If ACMG_DB_PATH pointed to a .db that doesn't exist, fall back to FAISS default
            if db_path.suffix == ".db":
                db_path = VECTORDB_GRCH37_NGSGENES if genome_build == "GRCh37" else VECTORDB_GRCH38_NGSGENES
            _databases[genome_build] = VariantVectorStore(db_path=db_path)
            print(f"Loaded {genome_build} FAISS database with {_databases[genome_build].count()} variants")

    return _databases[genome_build]


def generate_classification(variant_input: str, model_info: dict) -> str:
    """Generate ACMG classification using the model."""
    # Instruction matching training format
    instruction = "Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria."

    if model_info["type"] == "llama_cpp":
        # llama-cpp-python backend (GGUF model for Docker deployment)
        llm = model_info["model"]

        # System prompt for ACMG classification
        system_prompt = """You are an expert clinical geneticist performing ACMG/AMP variant classification.
Analyze the variant data provided and classify it according to the ACMG/AMP 2015 guidelines.
Provide the classification (Pathogenic, Likely_pathogenic, Uncertain_significance, Likely_benign, or Benign)
followed by the specific ACMG criteria that support your classification."""

        # Generate using chat format
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"{instruction}\n\n{variant_input}"}
        ]

        response = llm.create_chat_completion(
            messages=messages,
            max_tokens=512,
            temperature=0.1,
            top_p=0.9,
            repeat_penalty=1.1,
        )

        return response["choices"][0]["message"]["content"]

    elif model_info["type"] == "ollama":
        # Ollama uses raw prompt (not fine-tuned)
        prompt = f"{instruction}\n\n{variant_input}\n\nACMG Classification:"
        import requests
        ollama_host = os.environ.get("OLLAMA_HOST", "localhost:11434")
        if not ollama_host.startswith("http"):
            ollama_host = f"http://{ollama_host}"
        response = requests.post(
            f"{ollama_host}/api/generate",
            json={
                "model": model_info["model"],
                "prompt": prompt,
                "stream": False,
            },
            timeout=120,
        )
        return response.json().get("response", "")

    elif model_info["type"] == "mlx":
        from mlx_lm import generate
        from mlx_lm.sample_utils import make_sampler
        model = model_info["model"]
        tokenizer = model_info["tokenizer"]

        # Apply chat template to match training format
        messages = [
            {"role": "user", "content": f"{instruction}\n\n{variant_input}"}
        ]
        prompt = tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=True
        )
        # Use low temperature for consistent classification
        sampler = make_sampler(temp=0.1)
        return generate(model, tokenizer, prompt=prompt, max_tokens=500, sampler=sampler)

    elif model_info["type"] == "transformers":
        model = model_info["model"]
        tokenizer = model_info["tokenizer"]
        device = model_info["device"]

        # Apply chat template to match training format
        messages = [
            {"role": "user", "content": f"{instruction}\n\n{variant_input}"}
        ]
        prompt = tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=True
        )

        inputs = tokenizer(prompt, return_tensors="pt")
        if device == "cuda":
            inputs = {k: v.cuda() for k, v in inputs.items()}

        outputs = model.generate(
            **inputs,
            max_new_tokens=500,
            temperature=0.1,  # Low temperature for consistent output
            do_sample=True,
        )
        # Extract only the generated response (after the prompt)
        full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
        # Return just the assistant's response
        if "assistant" in full_output.lower():
            return full_output.split("assistant")[-1].strip()
        return full_output

    return "Classification unavailable"


def parse_classification_response(response: str) -> tuple[str, list[dict], float]:
    """Parse model response into structured data."""
    acmg_class = "Uncertain_significance"
    criteria = []
    confidence = 0.5

    response_lower = response.lower()

    # Extract classification
    if "pathogenic" in response_lower:
        if "likely pathogenic" in response_lower or "likely_pathogenic" in response_lower:
            acmg_class = "Likely_pathogenic"
            confidence = 0.75
        elif "benign" not in response_lower:
            acmg_class = "Pathogenic"
            confidence = 0.9

    if "benign" in response_lower:
        if "likely benign" in response_lower or "likely_benign" in response_lower:
            acmg_class = "Likely_benign"
            confidence = 0.75
        elif "pathogenic" not in response_lower:
            acmg_class = "Benign"
            confidence = 0.9

    # Extract criteria codes
    import re
    criteria_pattern = r'(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])[:.]?\s*([^.\n]+)'
    for match in re.finditer(criteria_pattern, response, re.IGNORECASE):
        criteria.append({
            "code": match.group(1).upper(),
            "description": match.group(2).strip()
        })

    return acmg_class, criteria, confidence


def generate_acmg_interpretation(
    acmg_score,
    gene: str,
    variant_id: str,
    scores: dict,
    clinvar_sig: str | None = None
) -> str:
    """Generate a clinical interpretation from ACMG rule-based scoring.

    Creates a clear, evidence-based interpretation that aligns with the
    ACMG/AMP guidelines and the criteria actually applied.
    """
    classification = acmg_score.classification.value.replace("_", " ")
    criteria_met = acmg_score.criteria_met
    rule = acmg_score.rule_applied

    # Build interpretation text
    lines = []

    # Classification summary
    lines.append(f"**ACMG Classification: {classification}**")
    lines.append("")

    if rule:
        lines.append(f"*Classification Rule: {rule}*")
        lines.append("")

    # Gene context
    if gene:
        lines.append(f"This variant in **{gene}** has been classified as **{classification}** based on the following evidence:")
    else:
        lines.append(f"This variant has been classified as **{classification}** based on the following evidence:")
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
        # Check concordance
        clinvar_lower = clinvar_sig.lower()
        classification_lower = classification.lower()
        if clinvar_lower in classification_lower or classification_lower in clinvar_lower:
            lines.append("✓ Rule-based classification is concordant with ClinVar.")
        elif "pathogenic" in clinvar_lower and "benign" in classification_lower:
            lines.append("⚠ Discordance: Rule-based classification differs from ClinVar. Expert review recommended.")
        elif "benign" in clinvar_lower and "pathogenic" in classification_lower:
            lines.append("⚠ Discordance: Rule-based classification differs from ClinVar. Expert review recommended.")
        lines.append("")

    # Key scores summary
    if scores:
        lines.append("**Key Scores:**")
        score_items = []
        if "gnomad_af" in scores and scores["gnomad_af"] is not None:
            af = scores["gnomad_af"]
            if af >= 0.01:
                score_items.append(f"gnomAD AF: {af*100:.2f}% (common)")
            elif af >= 0.001:
                score_items.append(f"gnomAD AF: {af*100:.3f}% (low frequency)")
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
    lines.append("*This classification is based on automated ACMG/AMP rule application. Expert review is required before clinical use.*")

    return "\n".join(lines)


@app.get("/health", response_model=HealthResponse)
async def health_check(genome_build: str = Query(default="GRCh37", description="Genome build")):
    """Check API health and model status."""
    model = get_model()
    try:
        store = get_database(genome_build)
        db_loaded = True
        count = store.count()
    except Exception:
        db_loaded = False
        count = 0

    return HealthResponse(
        status="healthy",
        model_loaded=model is not None,
        database_loaded=db_loaded,
        variant_count=count,
        genome_build=genome_build,
    )


@app.get("/classify", response_model=ClassificationResponse)
async def classify_variant_get(
    chr: str = Query(..., description="Chromosome"),
    pos: int = Query(..., description="Position"),
    ref: str = Query(..., description="Reference allele"),
    alt: str = Query(..., description="Alternate allele"),
    genome_build: str = Query(default="GRCh37", description="Genome build (GRCh37 or GRCh38)"),
):
    """Classify a variant using GET parameters."""
    return await classify_variant(VariantRequest(chr=chr, pos=pos, ref=ref, alt=alt, genome_build=genome_build))


@app.post("/classify", response_model=ClassificationResponse)
async def classify_variant(request: VariantRequest):
    """
    Classify a genetic variant according to ACMG/AMP guidelines.

    Returns the 5-tier classification (Pathogenic, Likely pathogenic,
    Uncertain significance, Likely benign, Benign) along with
    supporting evidence criteria.
    """
    # Validate genome build
    genome_build = request.genome_build
    if genome_build not in SUPPORTED_BUILDS:
        raise HTTPException(status_code=400, detail=f"Unsupported genome build: {genome_build}. Must be GRCh37 or GRCh38")

    # Validate input
    validation = validate_variant(request.chr, request.pos, request.ref, request.alt)
    if not validation.valid:
        raise HTTPException(
            status_code=400,
            detail={
                "error": validation.error,
                "suggestion": "Check chromosome (1-22, X, Y, M), position (positive integer), alleles (ACGT only)"
            }
        )

    # Use normalized variant ID from validation
    variant_id = validation.variant_id
    # Extract normalized chromosome from variant_id (format: chr_pos_ref_alt)
    chrom = variant_id.split("_")[0]

    # Try to find in database first
    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    scores = {}
    gene = None
    ref_normalized = request.ref.upper()
    alt_normalized = request.alt.upper()

    if variant_data:
        meta = variant_data["metadata"]

        # Validate reference allele against database
        stored_ref = meta.get("ref")
        if stored_ref:
            ref_check = validate_reference_against_database(variant_id, request.ref, stored_ref)
            if not ref_check.valid:
                raise HTTPException(
                    status_code=400,
                    detail={
                        "error": ref_check.error,
                        "expected_ref": stored_ref,
                        "provided_ref": request.ref,
                        "suggestion": "Verify your coordinates match the genome build"
                    }
                )

        gene = meta.get("gene")
        scores = {
            "cadd_phred": meta.get("cadd_phred"),
            "revel_score": meta.get("revel_score"),
            "gnomad_af": meta.get("gnomad_af"),
            "sift_pred": meta.get("sift_pred"),
            "polyphen_pred": meta.get("polyphen_pred"),
        }
        scores = {k: v for k, v in scores.items() if v is not None}

        # Extract ClinVar data if available
        clinvar_sig = meta.get("clinvar_sig")
        if clinvar_sig:
            clinvar_data = ClinVarAnnotation(
                id=meta.get("clinvar_id") or None,
                significance=clinvar_sig or None,
                review_status=meta.get("clinvar_review") or None,
                trait=meta.get("clinvar_trait") or None,
            )
        else:
            clinvar_data = None

        # Use rule-based ACMG scoring as PRIMARY classification
        acmg_score = score_variant_from_metadata(meta)

        # Convert criteria_met to response format
        criteria_met_details = [
            ACMGCriteriaDetail(
                code=c.code,
                name=c.name,
                short_name=c.short_name,
                description=c.description,
                strength=c.strength.value,
                evidence_type=c.evidence_type.value,
                category=c.category,
                status=c.status.value,
                evidence=c.evidence,
            )
            for c in acmg_score.criteria_met
        ]

        # Evaluate ALL 28 criteria with explanations
        # Pass ClinVar client for PS1/PM5 evaluation if available
        all_criteria_list = evaluate_all_criteria(meta, clinvar_client=CLINVAR_CLIENT)
        all_criteria_details = [
            ACMGCriteriaDetail(
                code=c.code,
                name=c.name,
                short_name=c.short_name,
                description=c.description,
                strength=c.strength.value,
                evidence_type=c.evidence_type.value,
                category=c.category,
                status=c.status.value,
                evidence=c.evidence,
            )
            for c in all_criteria_list
        ]

        # Generate interpretation from rule-based scoring
        clinvar_sig = meta.get("clinvar_sig")
        interpretation = generate_acmg_interpretation(
            acmg_score=acmg_score,
            gene=gene,
            variant_id=variant_id,
            scores=scores,
            clinvar_sig=clinvar_sig
        )

        # Generate external database links for evidence verification
        evidence_links = generate_evidence_links(meta, build=genome_build.lower())

        return ClassificationResponse(
            variant_id=variant_id,
            gene=gene,
            acmg_classification=acmg_score.classification.value,
            criteria=[],  # Rule-based scoring uses criteria_met instead
            criteria_met=criteria_met_details,
            all_criteria=all_criteria_details,
            rule_applied=acmg_score.rule_applied,
            scoring_method="rule_based",
            confidence=acmg_score.confidence,
            interpretation=interpretation,
            scores=scores,
            clinvar=clinvar_data,
            evidence_links=evidence_links,
            genome_build=genome_build,
            scoring_notes=acmg_score.notes,
        )
    else:
        # Variant not in database - use LLM fallback
        clinvar_data = None
        variant_input = f"""Variant: chr{chrom}:{request.pos} {ref_normalized}>{alt_normalized}
Note: This variant is not in the NGSgenes database. Limited evidence available."""

        # Get model and generate classification (LLM fallback)
        model = get_model()
        if model is None:
            raise HTTPException(status_code=503, detail="Model not available")

        interpretation = generate_classification(variant_input, model)
        acmg_class, criteria_list, confidence = parse_classification_response(interpretation)

        # Generate minimal evidence links for LLM fallback (basic variant info only)
        fallback_meta = {
            "chr": chrom,
            "pos": request.pos,
            "ref": ref_normalized,
            "alt": alt_normalized,
            "gene": gene,
        }
        evidence_links = generate_evidence_links(fallback_meta, build=genome_build.lower())

        return ClassificationResponse(
            variant_id=variant_id,
            gene=gene,
            acmg_classification=acmg_class,
            criteria=[ACMGCriteria(code=c["code"], description=c["description"]) for c in criteria_list],
            criteria_met=[],  # No rule-based criteria for fallback
            all_criteria=[],  # No criteria evaluation for LLM fallback
            rule_applied="",
            scoring_method="llm_fallback",
            confidence=confidence,
            interpretation=interpretation,
            scores=scores,
            clinvar=clinvar_data,
            evidence_links=evidence_links,
            genome_build=genome_build,
            scoring_notes=["Variant not found in database - using LLM classification as fallback"],
        )


@app.get("/variant/{variant_id}")
async def get_variant_info(
    variant_id: str,
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """Get raw variant information from database."""
    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(status_code=404, detail="Variant not found")

    return variant_data


def parse_hgvs_g(hgvs: str) -> tuple[str, int, str, str]:
    """Parse HGVS genomic notation (e.g., NC_000017.10:g.41197801T>A or 17:g.41197801T>A)."""
    import re

    # Pattern: optional NC_XXXXXX.X: or chrX: prefix, then g.POS REF>ALT
    pattern = r'(?:NC_0+(\d+)(?:\.\d+)?:)?(?:chr)?(\d+|X|Y)?:?g\.(\d+)([ACGT]+)>([ACGT]+)'
    match = re.match(pattern, hgvs.strip(), re.IGNORECASE)

    if not match:
        raise ValueError(f"Invalid HGVS format: {hgvs}")

    nc_chr, simple_chr, pos, ref, alt = match.groups()
    chrom = nc_chr or simple_chr
    if not chrom:
        raise ValueError(f"Could not parse chromosome from: {hgvs}")

    return chrom, int(pos), ref.upper(), alt.upper()


@app.get("/lookup")
async def lookup_variant(
    chr: str = Query(None, description="Chromosome"),
    pos: int = Query(None, description="Position"),
    ref: str = Query(None, description="Reference allele"),
    alt: str = Query(None, description="Alternate allele"),
    hgvs: str = Query(None, description="HGVS genomic notation (e.g., NC_000017.10:g.41197801T>A or 17:g.41197801T>A)"),
    genome_build: str = Query(default="GRCh37", description="Genome build (GRCh37 or GRCh38)"),
):
    """
    Get ALL stored dbNSFP data for a variant.

    Accepts either coordinate parameters (chr, pos, ref, alt) or HGVS genomic notation.
    Returns the full annotation document with all scores, predictions, and clinical data.
    """
    # Parse input - either HGVS or coordinates
    if hgvs:
        try:
            chrom, pos, ref, alt = parse_hgvs_g(hgvs)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
    elif chr and pos and ref and alt:
        chrom = chr
        # Values will be validated below
    else:
        raise HTTPException(status_code=400, detail="Provide either hgvs parameter or chr/pos/ref/alt parameters")

    # Validate input
    validation = validate_variant(chrom, pos, ref, alt)
    if not validation.valid:
        raise HTTPException(
            status_code=400,
            detail={
                "error": validation.error,
                "suggestion": "Check chromosome (1-22, X, Y, M), position (positive integer), alleles (ACGT only)"
            }
        )

    variant_id = validation.variant_id

    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(
            status_code=404,
            detail={
                "error": f"Variant {variant_id} not found in {genome_build} database",
                "suggestion": "Verify coordinates match the genome build, or try using HGVS notation"
            }
        )

    meta = variant_data["metadata"]

    # Validate reference allele against database
    stored_ref = meta.get("ref")
    user_ref = ref if not hgvs else ref  # ref from either coordinates or HGVS parsing
    if stored_ref and user_ref:
        ref_check = validate_reference_against_database(variant_id, user_ref, stored_ref)
        if not ref_check.valid:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": ref_check.error,
                    "expected_ref": stored_ref,
                    "provided_ref": user_ref,
                    "suggestion": "Verify your coordinates match the genome build"
                }
            )

    document = variant_data.get("document", "")

    # Parse normalized values from variant_id
    parts = variant_id.split("_")
    norm_chrom, norm_pos, norm_ref, norm_alt = parts[0], int(parts[1]), parts[2], parts[3]

    # Extract structured ClinVar data
    clinvar_sig = meta.get("clinvar_sig")
    clinvar_data = None
    if clinvar_sig:
        clinvar_data = {
            "id": meta.get("clinvar_id") or None,
            "significance": clinvar_sig or None,
            "review_status": meta.get("clinvar_review") or None,
            "trait": meta.get("clinvar_trait") or None,
        }

    # Generate external database links for evidence verification
    evidence_links = generate_evidence_links(meta, build=genome_build.lower())

    return {
        "variant_id": variant_id,
        "chromosome": norm_chrom,
        "position": norm_pos,
        "ref": norm_ref,
        "alt": norm_alt,
        "gene": meta.get("gene"),
        "transcript": meta.get("transcript"),
        "clinvar": clinvar_data,
        "metadata": meta,
        "full_annotation": document,
        "evidence_links": evidence_links,
        "genome_build": genome_build,
        "dbnsfp_version": "5.3.1a",
        "disclaimer": RUO_DISCLAIMER,
    }


@app.get("/gene/{gene_symbol}")
async def get_gene_variants(
    gene_symbol: str,
    limit: int = 100,
    genome_build: str = Query(default="GRCh37", description="Genome build"),
):
    """Get variants for a specific gene."""
    store = get_database(genome_build)
    results = store.search_by_gene(gene_symbol.upper(), k=limit)

    return {
        "gene": gene_symbol.upper(),
        "genome_build": genome_build,
        "count": len(results),
        "variants": [
            {
                "id": r["id"],
                "clinvar_sig": r["metadata"].get("clinvar_sig"),
                "cadd_phred": r["metadata"].get("cadd_phred"),
            }
            for r in results
        ]
    }


@app.get("/panels")
async def list_panels():
    """List available gene panels and variant counts."""
    from src.panels import list_panels as _list_panels

    panels_info = _list_panels()

    # Try to get variant counts from database if available
    panel_counts = {}
    try:
        db = get_database("GRCh37")
        if hasattr(db, "count_by_panel"):
            for name in panels_info:
                panel_counts[name] = db.count_by_panel(name)
    except Exception:
        pass

    return {
        "panels": [
            {
                "name": name,
                "gene_count": gene_count,
                "variant_count": panel_counts.get(name),
            }
            for name, gene_count in panels_info.items()
        ]
    }


@app.get("/panel/{panel_name}/genes")
async def get_panel_genes(panel_name: str):
    """Get the list of genes in a specific panel."""
    from src.panels import get_panel

    try:
        genes = get_panel(panel_name)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "panel": panel_name,
        "gene_count": len(genes),
        "genes": sorted(genes),
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

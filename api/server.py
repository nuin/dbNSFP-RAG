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
from pydantic import BaseModel, Field

# Add project root to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

app = FastAPI(
    title="ACMG Variant Classifier",
    description="Classify genetic variants according to ACMG/AMP guidelines",
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


# Request/Response models
class VariantRequest(BaseModel):
    chr: str = Field(..., description="Chromosome (1-22, X, Y)")
    pos: int = Field(..., description="Position (1-based)")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")

    class Config:
        json_schema_extra = {
            "example": {
                "chr": "17",
                "pos": 41197801,
                "ref": "T",
                "alt": "A"
            }
        }


class ACMGCriteria(BaseModel):
    code: str = Field(..., description="ACMG criteria code (e.g., PVS1, PS1, PM2)")
    description: str = Field(..., description="Evidence description")


class ClassificationResponse(BaseModel):
    variant_id: str = Field(..., description="Variant identifier (chr_pos_ref_alt)")
    gene: Optional[str] = Field(None, description="Gene symbol")
    acmg_classification: str = Field(..., description="ACMG 5-tier classification")
    criteria: list[ACMGCriteria] = Field(default_factory=list, description="Applied ACMG criteria")
    confidence: float = Field(..., description="Classification confidence (0-1)")
    interpretation: str = Field(..., description="Full interpretation text")
    scores: dict = Field(default_factory=dict, description="Pathogenicity scores")


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
    database_loaded: bool
    variant_count: int


# Global model and database instances
_model = None
_vectorstore = None


def get_model():
    """Load fine-tuned model (lazy loading)."""
    global _model
    if _model is None:
        model_path = os.environ.get("ACMG_MODEL_PATH", "models/acmg-classifier/model")

        # Try mlx-lm first (Apple Silicon)
        try:
            from mlx_lm import load, generate
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


def get_vectorstore():
    """Load vector database (lazy loading)."""
    global _vectorstore
    if _vectorstore is None:
        from src.vectorstore import VariantVectorStore
        from src.config import VECTORDB_GRCH37_NGSGENES

        db_path = os.environ.get("ACMG_DB_PATH", str(VECTORDB_GRCH37_NGSGENES))
        _vectorstore = VariantVectorStore(db_path=Path(db_path))
        print(f"Loaded database with {_vectorstore.count()} variants")

    return _vectorstore


def generate_classification(variant_input: str, model_info: dict) -> str:
    """Generate ACMG classification using the model."""
    # Use exact format from training data
    prompt = f"""Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria.

{variant_input}

ACMG Classification:"""

    if model_info["type"] == "ollama":
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
        model = model_info["model"]
        tokenizer = model_info["tokenizer"]
        return generate(model, tokenizer, prompt=prompt, max_tokens=300)

    elif model_info["type"] == "transformers":
        model = model_info["model"]
        tokenizer = model_info["tokenizer"]
        device = model_info["device"]

        inputs = tokenizer(prompt, return_tensors="pt")
        if device == "cuda":
            inputs = {k: v.cuda() for k, v in inputs.items()}

        outputs = model.generate(
            **inputs,
            max_new_tokens=500,
            temperature=0.7,
            do_sample=True,
        )
        return tokenizer.decode(outputs[0], skip_special_tokens=True)

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


@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Check API health and model status."""
    model = get_model()
    store = get_vectorstore()

    return HealthResponse(
        status="healthy",
        model_loaded=model is not None,
        database_loaded=store is not None,
        variant_count=store.count() if store else 0,
    )


@app.get("/classify", response_model=ClassificationResponse)
async def classify_variant_get(
    chr: str = Query(..., description="Chromosome"),
    pos: int = Query(..., description="Position"),
    ref: str = Query(..., description="Reference allele"),
    alt: str = Query(..., description="Alternate allele"),
):
    """Classify a variant using GET parameters."""
    return await classify_variant(VariantRequest(chr=chr, pos=pos, ref=ref, alt=alt))


@app.post("/classify", response_model=ClassificationResponse)
async def classify_variant(request: VariantRequest):
    """
    Classify a genetic variant according to ACMG/AMP guidelines.

    Returns the 5-tier classification (Pathogenic, Likely pathogenic,
    Uncertain significance, Likely benign, Benign) along with
    supporting evidence criteria.
    """
    # Normalize chromosome
    chrom = request.chr.replace("chr", "")
    variant_id = f"{chrom}_{request.pos}_{request.ref}_{request.alt}"

    # Try to find in database first
    store = get_vectorstore()
    variant_data = store.get_by_id(variant_id)

    scores = {}
    gene = None

    if variant_data:
        meta = variant_data["metadata"]
        gene = meta.get("gene")
        scores = {
            "cadd_phred": meta.get("cadd_phred"),
            "revel_score": meta.get("revel_score"),
            "gnomad_af": meta.get("gnomad_af"),
            "sift_pred": meta.get("sift_pred"),
            "polyphen_pred": meta.get("polyphen_pred"),
        }
        scores = {k: v for k, v in scores.items() if v is not None}

        # Format input for model
        variant_input = f"""Variant: chr{chrom}:{request.pos} {request.ref}>{request.alt}
Gene: {gene}"""
        for key, value in scores.items():
            variant_input += f"\n{key}: {value}"
    else:
        # Variant not in database - classify with minimal info
        variant_input = f"""Variant: chr{chrom}:{request.pos} {request.ref}>{request.alt}
Note: This variant is not in the NGSgenes database. Limited evidence available."""

    # Get model and generate classification
    model = get_model()
    if model is None:
        raise HTTPException(status_code=503, detail="Model not available")

    interpretation = generate_classification(variant_input, model)
    acmg_class, criteria_list, confidence = parse_classification_response(interpretation)

    return ClassificationResponse(
        variant_id=variant_id,
        gene=gene,
        acmg_classification=acmg_class,
        criteria=[ACMGCriteria(code=c["code"], description=c["description"]) for c in criteria_list],
        confidence=confidence,
        interpretation=interpretation,
        scores=scores,
    )


@app.get("/variant/{variant_id}")
async def get_variant_info(variant_id: str):
    """Get raw variant information from database."""
    store = get_vectorstore()
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(status_code=404, detail="Variant not found")

    return variant_data


@app.get("/lookup")
async def lookup_variant(
    chr: str = Query(..., description="Chromosome"),
    pos: int = Query(..., description="Position"),
    ref: str = Query(..., description="Reference allele"),
    alt: str = Query(..., description="Alternate allele"),
):
    """
    Get all stored data for a variant without classification.

    Returns all metadata stored in the database including scores,
    predictions, gene info, and ClinVar significance.
    """
    chrom = chr.replace("chr", "")
    variant_id = f"{chrom}_{pos}_{ref}_{alt}"

    store = get_vectorstore()
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(status_code=404, detail=f"Variant {variant_id} not found in database")

    meta = variant_data["metadata"]

    return {
        "variant_id": variant_id,
        "chromosome": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "gene": meta.get("gene"),
        "clinvar": {
            "significance": meta.get("clinvar_sig"),
        },
        "scores": {
            "cadd_phred": meta.get("cadd_phred"),
            "revel_score": meta.get("revel_score"),
            "gnomad_af": meta.get("gnomad_af"),
        },
        "predictions": {
            "sift": meta.get("sift_pred"),
            "polyphen": meta.get("polyphen_pred"),
        },
        "text": variant_data.get("text", ""),
    }


@app.get("/gene/{gene_symbol}")
async def get_gene_variants(gene_symbol: str, limit: int = 100):
    """Get variants for a specific gene."""
    store = get_vectorstore()
    results = store.search_by_gene(gene_symbol.upper(), k=limit)

    return {
        "gene": gene_symbol.upper(),
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


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

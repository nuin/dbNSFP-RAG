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

import logging
import os
import sqlite3

# Add project root to path
import sys
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

sys.path.insert(0, str(Path(__file__).parent.parent))

from api import browse, lookup
from api.browse import (
    get_gene_variants,
    get_panel_genes,
    get_variant_info,
    list_chromosomes,
    list_genes_on_chromosome,
    list_panels,
)
from api.classification import build_classification_response, missing_variant_response
from api.hgvs import parse_hgvs_g
from api.inference import generate_classification, load_model, parse_classification_response
from api.interpretation import generate_acmg_interpretation
from api.lookup import lookup_variant
from api.models import (
    RUO_DISCLAIMER,
    ACMGCriteria,
    ACMGCriteriaDetail,
    ClassificationResponse,
    ClinVarAnnotation,
    HealthResponse,
    VariantRequest,
)
from src.validation import validate_variant

logger = logging.getLogger(__name__)

# Import ClinVar API client (optional, for PS1/PM5 evaluation)
try:
    from src.clinvar_api import ClinVarClient

    CLINVAR_CLIENT = ClinVarClient(email="dbnsfp-api@localhost")
except ImportError:
    CLINVAR_CLIENT = None

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
            str(index_file), headers={"Cache-Control": "no-cache, no-store, must-revalidate"}
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


# Global model and database instances
_model = None
_databases = {}  # Cache databases by genome build


def get_model():
    """Load and cache the optional inference backend."""
    global _model
    if _model is None:
        _model = load_model()
    return _model


def get_database(genome_build: str = "GRCh37"):
    """Load variant database for specified genome build (lazy loading).

    Uses SQLite backend (VariantDatabase) if the configured path is a .db file,
    otherwise falls back to FAISS (VariantVectorStore) for backwards compatibility.
    """

    if genome_build not in SUPPORTED_BUILDS:
        raise ValueError(
            f"Unsupported genome build: {genome_build}. Must be one of {SUPPORTED_BUILDS}"
        )

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
            print(
                f"Loaded {genome_build} SQLite database with {_databases[genome_build].count()} variants"
            )
        else:
            # Fall back to FAISS for backwards compatibility
            from src.vectorstore import VariantVectorStore

            # If ACMG_DB_PATH pointed to a .db that doesn't exist, fall back to FAISS default
            if db_path.suffix == ".db":
                db_path = (
                    VECTORDB_GRCH37_NGSGENES
                    if genome_build == "GRCh37"
                    else VECTORDB_GRCH38_NGSGENES
                )
            _databases[genome_build] = VariantVectorStore(db_path=db_path)
            print(
                f"Loaded {genome_build} FAISS database with {_databases[genome_build].count()} variants"
            )

    return _databases[genome_build]


@app.get("/health", response_model=HealthResponse)
async def health_check(genome_build: str = Query(default="GRCh37", description="Genome build")):
    """Check API health and model status."""
    try:
        store = get_database(genome_build)
        db_loaded = True
        count = store.count()
    except (OSError, ValueError, RuntimeError, ImportError, sqlite3.Error) as exc:
        logger.warning("Database health check failed for %s: %s", genome_build, exc)
        db_loaded = False
        count = 0

    return HealthResponse(
        status="healthy",
        model_loaded=True,  # Rule-based scoring, no LLM needed
        database_loaded=db_loaded,
        variant_count=count,
        genome_build=genome_build,
    )


@app.get("/classify", response_model=ClassificationResponse)
async def classify_variant_get(
    chr: str = Query(None, description="Chromosome"),
    pos: int = Query(None, description="Position"),
    ref: str = Query(None, description="Reference allele"),
    alt: str = Query(None, description="Alternate allele"),
    hgvs: str = Query(None, description="HGVS genomic notation (e.g., 17:g.41197801T>A)"),
    genome_build: str = Query(default="GRCh37", description="Genome build (GRCh37 or GRCh38)"),
):
    """Classify a variant using GET parameters or HGVS notation."""
    if hgvs:
        try:
            chrom, position, reference, alternate = parse_hgvs_g(hgvs)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        return await classify_variant(
            VariantRequest(
                chr=chrom, pos=position, ref=reference, alt=alternate, genome_build=genome_build
            )
        )
    elif chr and pos and ref and alt:
        return await classify_variant(
            VariantRequest(chr=chr, pos=pos, ref=ref, alt=alt, genome_build=genome_build)
        )
    else:
        raise HTTPException(
            status_code=400, detail="Provide either hgvs parameter or chr/pos/ref/alt parameters"
        )


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
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported genome build: {genome_build}. Must be GRCh37 or GRCh38",
        )

    # Validate input
    validation = validate_variant(request.chr, request.pos, request.ref, request.alt)
    if not validation.valid:
        raise HTTPException(
            status_code=400,
            detail={
                "error": validation.error,
                "suggestion": "Check chromosome (1-22, X, Y, M), position (positive integer), alleles (ACGT only)",
            },
        )

    # Use normalized variant ID from validation
    variant_id = validation.variant_id

    # Try to find in database first
    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    if variant_data:
        return build_classification_response(
            variant_data["metadata"], variant_id, request, CLINVAR_CLIENT
        )
    return missing_variant_response(variant_id, genome_build)


app.include_router(browse.router)
app.include_router(lookup.router)

__all__ = [
    "RUO_DISCLAIMER",
    "ACMGCriteria",
    "ACMGCriteriaDetail",
    "ClassificationResponse",
    "ClinVarAnnotation",
    "HealthResponse",
    "VariantRequest",
    "app",
    "classify_variant",
    "classify_variant_get",
    "generate_acmg_interpretation",
    "generate_classification",
    "get_database",
    "get_gene_variants",
    "get_model",
    "get_panel_genes",
    "get_variant_info",
    "health_check",
    "list_chromosomes",
    "list_genes_on_chromosome",
    "list_panels",
    "lookup_page",
    "lookup_variant",
    "parse_classification_response",
    "parse_hgvs_g",
    "root",
]

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)

"""Lookup full dbNSFP records by coordinate or HGVS notation."""

from fastapi import APIRouter, HTTPException, Query

from api.dependencies import get_database
from api.hgvs import parse_hgvs_g
from api.models import RUO_DISCLAIMER
from src.evidence_links import generate_evidence_links
from src.validation import validate_reference_against_database, validate_variant

router = APIRouter()


@router.get("/lookup")
async def lookup_variant(
    chr: str = Query(None, description="Chromosome"),
    pos: int = Query(None, description="Position"),
    ref: str = Query(None, description="Reference allele"),
    alt: str = Query(None, description="Alternate allele"),
    hgvs: str = Query(
        None,
        description="HGVS genomic notation (e.g., NC_000017.10:g.41197801T>A or 17:g.41197801T>A)",
    ),
    genome_build: str = Query(default="GRCh37", description="Genome build (GRCh37 or GRCh38)"),
):
    """
    Get ALL stored dbNSFP data for a variant.

    Accepts either coordinate parameters (chr, pos, ref, alt) or HGVS genomic notation.
    Returns the full annotation document with all scores, predictions, and clinical data.
    """
    if hgvs:
        try:
            chrom, pos, ref, alt = parse_hgvs_g(hgvs)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
    elif chr and pos and ref and alt:
        chrom = chr
        # Values will be validated below
    else:
        raise HTTPException(
            status_code=400, detail="Provide either hgvs parameter or chr/pos/ref/alt parameters"
        )

    validation = validate_variant(chrom, pos, ref, alt)
    if not validation.valid:
        raise HTTPException(
            status_code=400,
            detail={
                "error": validation.error,
                "suggestion": "Check chromosome (1-22, X, Y, M), position (positive integer), alleles (ACGT only)",
            },
        )

    variant_id = validation.variant_id

    store = get_database(genome_build)
    variant_data = store.get_by_id(variant_id)

    if not variant_data:
        raise HTTPException(
            status_code=404,
            detail={
                "error": f"Variant {variant_id} not found in {genome_build} database",
                "suggestion": "Verify coordinates match the genome build, or try using HGVS notation",
            },
        )

    meta = variant_data["metadata"]

    stored_ref = meta.get("ref")
    user_ref = ref
    if stored_ref and user_ref:
        ref_check = validate_reference_against_database(variant_id, user_ref, stored_ref)
        if not ref_check.valid:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": ref_check.error,
                    "expected_ref": stored_ref,
                    "provided_ref": user_ref,
                    "suggestion": "Verify your coordinates match the genome build",
                },
            )

    document = variant_data.get("document", "")

    parts = variant_id.split("_")
    norm_chrom, norm_pos, norm_ref, norm_alt = parts[0], int(parts[1]), parts[2], parts[3]

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

"""Document chunking module - converts variant rows to structured text."""

import pandas as pd


def format_score(value, precision: int = 4) -> str:
    """Format a score value, handling NaN and multi-value fields."""
    if pd.isna(value):
        return "N/A"
    # Handle multi-value fields (e.g., ".;.;.;." or "0.5;0.6;0.7")
    if isinstance(value, str):
        if value == "." or value.replace(";", "").replace(".", "") == "":
            return "N/A"
        # Take first non-empty value
        parts = [p.strip() for p in value.split(";") if p.strip() and p.strip() != "."]
        if not parts:
            return "N/A"
        try:
            return f"{float(parts[0]):.{precision}f}"
        except ValueError:
            return parts[0]
    if isinstance(value, float):
        return f"{value:.{precision}f}"
    return str(value)


def format_prediction(pred, score=None) -> str:
    """Format prediction with optional score."""
    if pd.isna(pred):
        return "N/A"
    if score is not None and not pd.isna(score):
        return f"{pred} ({format_score(score)})"
    return str(pred)


def _parse_first_float(value) -> float | None:
    """Parse first valid float from a potentially multi-value field."""
    if pd.isna(value):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        if value == "." or value.replace(";", "").replace(".", "") == "":
            return None
        parts = [p.strip() for p in value.split(";") if p.strip() and p.strip() != "."]
        if not parts:
            return None
        try:
            return float(parts[0])
        except ValueError:
            return None
    return None


def interpret_cadd(phred_score) -> str:
    """Interpret CADD phred score."""
    score = _parse_first_float(phred_score)
    if score is None:
        return "N/A"
    if score >= 30:
        return f"{score:.1f} (top 0.1% deleterious)"
    elif score >= 20:
        return f"{score:.1f} (top 1% deleterious)"
    elif score >= 10:
        return f"{score:.1f} (top 10% deleterious)"
    else:
        return f"{score:.1f} (likely benign)"


def interpret_revel(score) -> str:
    """Interpret REVEL score."""
    s = _parse_first_float(score)
    if s is None:
        return "N/A"
    if s >= 0.75:
        return f"{s:.3f} (likely pathogenic)"
    elif s >= 0.5:
        return f"{s:.3f} (uncertain)"
    else:
        return f"{s:.3f} (likely benign)"


def interpret_gnomad_af(af) -> str:
    """Interpret gnomAD allele frequency."""
    freq = _parse_first_float(af)
    if freq is None:
        return "N/A (not observed)"
    if freq == 0:
        return "0 (not observed)"
    elif freq < 0.0001:
        return f"{freq:.2e} (ultra-rare)"
    elif freq < 0.001:
        return f"{freq:.2e} (rare)"
    elif freq < 0.01:
        return f"{freq:.4f} (low frequency)"
    else:
        return f"{freq:.4f} (common)"


def variant_to_text(row: pd.Series) -> str:
    """
    Convert a variant row to structured text for embedding.

    This format is optimized for:
    1. Semantic search (natural language descriptions)
    2. LLM context (structured but readable)
    3. Clinical interpretation
    """
    # Build variant ID
    chrom = row.get("#chr", "?")
    pos = row.get("pos(1-based)", "?")
    ref = row.get("ref", "?")
    alt = row.get("alt", "?")

    # Gene info
    gene = row.get("genename", "N/A")
    transcript = row.get("Ensembl_transcriptid", "")
    protein = row.get("Ensembl_proteinid", "")

    # Amino acid change
    aa_ref = row.get("aaref", "")
    aa_alt = row.get("aaalt", "")
    aa_pos = row.get("aapos", "")

    aa_change = "N/A"
    if aa_ref and aa_alt and aa_pos and not pd.isna(aa_ref):
        aa_change = f"p.{aa_ref}{aa_pos}{aa_alt}"

    # HGVS notation
    hgvsc = row.get("HGVSc_snpEff", "N/A")
    hgvsp = row.get("HGVSp_snpEff", "N/A")

    text = f"""Variant: chr{chrom}:{pos} {ref}>{alt}
Gene: {gene}
Transcript: {transcript if not pd.isna(transcript) else 'N/A'}
Protein: {protein if not pd.isna(protein) else 'N/A'}
HGVS coding: {hgvsc if not pd.isna(hgvsc) else 'N/A'}
HGVS protein: {hgvsp if not pd.isna(hgvsp) else 'N/A'}
Amino acid change: {aa_change}

=== Pathogenicity Predictions ===
SIFT: {format_prediction(row.get('SIFT_pred'), row.get('SIFT_score'))}
PolyPhen2 HDIV: {format_prediction(row.get('Polyphen2_HDIV_pred'), row.get('Polyphen2_HDIV_score'))}
PolyPhen2 HVAR: {format_prediction(row.get('Polyphen2_HVAR_pred'), row.get('Polyphen2_HVAR_score'))}
CADD: {interpret_cadd(row.get('CADD_phred'))}
REVEL: {interpret_revel(row.get('REVEL_score'))}
AlphaMissense: {format_prediction(row.get('AlphaMissense_pred'), row.get('AlphaMissense_score'))}
ClinPred: {format_prediction(row.get('ClinPred_pred'), row.get('ClinPred_score'))}
DANN: {format_score(row.get('DANN_score'))}
MetaSVM: {format_prediction(row.get('MetaSVM_pred'), row.get('MetaSVM_score'))}
MetaLR: {format_prediction(row.get('MetaLR_pred'), row.get('MetaLR_score'))}

=== Conservation ===
phyloP (100-way vertebrate): {format_score(row.get('phyloP100way_vertebrate'))}
phyloP (30-way mammalian): {format_score(row.get('phyloP30way_mammalian'))}
phastCons (100-way): {format_score(row.get('phastCons100way_vertebrate'))}
GERP++ RS: {format_score(row.get('GERP++_RS'))}

=== Population Frequency ===
gnomAD exomes: {interpret_gnomad_af(row.get('gnomAD_exomes_AF'))}
gnomAD genomes: {interpret_gnomad_af(row.get('gnomAD_genomes_AF'))}
1000 Genomes: {interpret_gnomad_af(row.get('1000Gp3_AF'))}

=== Clinical Annotation ===
ClinVar ID: {row.get('clinvar_id', 'N/A') if not pd.isna(row.get('clinvar_id')) else 'N/A'}
ClinVar significance: {row.get('clinvar_clnsig', 'N/A') if not pd.isna(row.get('clinvar_clnsig')) else 'N/A'}
ClinVar review status: {row.get('clinvar_review', 'N/A') if not pd.isna(row.get('clinvar_review')) else 'N/A'}
ClinVar trait: {row.get('clinvar_trait', 'N/A') if not pd.isna(row.get('clinvar_trait')) else 'N/A'}

=== Functional Annotation ===
InterPro domain: {row.get('Interpro_domain', 'N/A') if not pd.isna(row.get('Interpro_domain')) else 'N/A'}
GTEx gene: {row.get('GTEx_V8_gene', 'N/A') if not pd.isna(row.get('GTEx_V8_gene')) else 'N/A'}
GTEx tissue: {row.get('GTEx_V8_tissue', 'N/A') if not pd.isna(row.get('GTEx_V8_tissue')) else 'N/A'}
"""
    return text.strip()


def variant_to_id(row: pd.Series) -> str:
    """Generate unique variant ID."""
    chrom = str(row.get("#chr", "")).replace("chr", "")
    pos = row.get("pos(1-based)", "")
    ref = row.get("ref", "")
    alt = row.get("alt", "")
    return f"{chrom}_{pos}_{ref}_{alt}"


def variant_to_metadata(row: pd.Series) -> dict:
    """Extract metadata for filtering/retrieval."""

    def safe_str(val):
        return str(val) if not pd.isna(val) else ""

    def safe_float(val):
        if pd.isna(val):
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    return {
        "chr": safe_str(row.get("#chr", "")).replace("chr", ""),
        "pos": int(row.get("pos(1-based)", 0)) if not pd.isna(row.get("pos(1-based)")) else 0,
        "ref": safe_str(row.get("ref")),
        "alt": safe_str(row.get("alt")),
        "gene": safe_str(row.get("genename")),
        "transcript": safe_str(row.get("Ensembl_transcriptid")),
        "cadd_phred": safe_float(row.get("CADD_phred")),
        "revel_score": safe_float(row.get("REVEL_score")),
        "clinvar_sig": safe_str(row.get("clinvar_clnsig")),
        "gnomad_af": safe_float(row.get("gnomAD_exomes_AF")),
        "sift_pred": safe_str(row.get("SIFT_pred")),
        "polyphen_pred": safe_str(row.get("Polyphen2_HDIV_pred")),
        "alphamissense_pred": safe_str(row.get("AlphaMissense_pred")),
    }


def chunk_dataframe(df: pd.DataFrame) -> list[tuple[str, str, dict]]:
    """
    Convert DataFrame to list of (id, text, metadata) tuples.

    Returns:
        List of tuples: (variant_id, text_chunk, metadata_dict)
    """
    results = []
    for _, row in df.iterrows():
        var_id = variant_to_id(row)
        text = variant_to_text(row)
        metadata = variant_to_metadata(row)
        results.append((var_id, text, metadata))
    return results

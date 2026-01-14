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


def variant_to_text(row: pd.Series, use_grch37: bool = False) -> str:
    """
    Convert a variant row to structured text for embedding.

    This format is optimized for:
    1. Semantic search (natural language descriptions)
    2. LLM context (structured but readable)
    3. Clinical interpretation
    """
    # Build variant ID - use GRCh37 coordinates if specified
    if use_grch37:
        chrom = row.get("hg19_chr", row.get("#chr", "?"))
        pos = row.get("hg19_pos(1-based)", "?")
    else:
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
MutationTaster: {format_prediction(row.get('MutationTaster_pred'), row.get('MutationTaster_score'))}
BayesDel: {format_prediction(row.get('BayesDel_addAF_pred'), row.get('BayesDel_addAF_score'))}
PROVEAN: {format_prediction(row.get('PROVEAN_pred'), row.get('PROVEAN_score'))}

=== Conservation ===
phyloP (100-way vertebrate): {format_score(row.get('phyloP100way_vertebrate'))}
phyloP (470-way mammalian): {format_score(row.get('phyloP470way_mammalian'))}
phastCons (100-way): {format_score(row.get('phastCons100way_vertebrate'))}
GERP++ RS: {format_score(row.get('GERP++_RS'))}

=== Population Frequency ===
gnomAD v4 (joint): {interpret_gnomad_af(row.get('gnomAD4.1_joint_AF'))}
gnomAD v2 exomes: {interpret_gnomad_af(row.get('gnomAD2.1.1_exomes_controls_AF'))}
1000 Genomes: {interpret_gnomad_af(row.get('1000Gp3_AF'))}

=== Clinical Annotation ===
ClinVar ID: {row.get('clinvar_id', 'N/A') if not pd.isna(row.get('clinvar_id')) else 'N/A'}
ClinVar significance: {row.get('clinvar_clnsig', 'N/A') if not pd.isna(row.get('clinvar_clnsig')) else 'N/A'}
ClinVar review status: {row.get('clinvar_review', 'N/A') if not pd.isna(row.get('clinvar_review')) else 'N/A'}
ClinVar trait: {row.get('clinvar_trait', 'N/A') if not pd.isna(row.get('clinvar_trait')) else 'N/A'}

=== Functional Annotation ===
InterPro domain: {row.get('Interpro_domain', 'N/A') if not pd.isna(row.get('Interpro_domain')) else 'N/A'}
"""
    return text.strip()


def variant_to_id(row: pd.Series, use_grch37: bool = False) -> str:
    """Generate unique variant ID."""
    if use_grch37:
        chrom = str(row.get("hg19_chr", row.get("#chr", ""))).replace("chr", "")
        pos = row.get("hg19_pos(1-based)", "")
    else:
        chrom = str(row.get("#chr", "")).replace("chr", "")
        pos = row.get("pos(1-based)", "")
    ref = row.get("ref", "")
    alt = row.get("alt", "")
    return f"{chrom}_{pos}_{ref}_{alt}"


def variant_to_metadata(row: pd.Series, use_grch37: bool = False) -> dict:
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

    if use_grch37:
        chrom = safe_str(row.get("hg19_chr", row.get("#chr", ""))).replace("chr", "")
        pos_val = row.get("hg19_pos(1-based)", row.get("pos(1-based)", 0))
    else:
        chrom = safe_str(row.get("#chr", "")).replace("chr", "")
        pos_val = row.get("pos(1-based)", 0)

    # Prefer gnomAD v4, fall back to v2
    gnomad_af = safe_float(row.get("gnomAD4.1_joint_AF"))
    if gnomad_af is None:
        gnomad_af = safe_float(row.get("gnomAD2.1.1_exomes_controls_AF"))

    # Amino acid information for PS1/PM5 evaluation
    aa_pos_val = row.get("aapos")
    aa_pos_int = None
    if not pd.isna(aa_pos_val):
        try:
            aa_pos_int = int(str(aa_pos_val).split(";")[0])  # Take first if multi-value
        except (ValueError, TypeError):
            pass

    return {
        "chr": chrom,
        "pos": int(pos_val) if not pd.isna(pos_val) else 0,
        "ref": safe_str(row.get("ref")),
        "alt": safe_str(row.get("alt")),
        "gene": safe_str(row.get("genename")),
        "transcript": safe_str(row.get("Ensembl_transcriptid")),
        # Amino acid data (for PS1/PM5 ClinVar lookup)
        "aa_ref": safe_str(row.get("aaref")),
        "aa_alt": safe_str(row.get("aaalt")),
        "aa_pos": aa_pos_int,
        "cadd_phred": safe_float(row.get("CADD_phred")),
        "revel_score": safe_float(row.get("REVEL_score")),
        # ClinVar annotations
        "clinvar_id": safe_str(row.get("clinvar_id")),
        "clinvar_sig": safe_str(row.get("clinvar_clnsig")),
        "clinvar_review": safe_str(row.get("clinvar_review")),
        "clinvar_trait": safe_str(row.get("clinvar_trait")),
        # Population frequency
        "gnomad_af": gnomad_af,
        # Predictor results
        "sift_pred": safe_str(row.get("SIFT_pred")),
        "polyphen_pred": safe_str(row.get("Polyphen2_HDIV_pred")),
        "alphamissense_pred": safe_str(row.get("AlphaMissense_pred")),
        "mutationtaster_pred": safe_str(row.get("MutationTaster_pred")),
        "bayesdel_pred": safe_str(row.get("BayesDel_addAF_pred")),
        "provean_pred": safe_str(row.get("PROVEAN_pred")),
        # Consequence/effect type (for PVS1, PM4, BP7)
        "consequence": safe_str(row.get("Ensembl_consequence")),
        # Gene constraint scores (for PP2, BP1, PVS1)
        "gnomad_pli": safe_float(row.get("gnomAD_pLI")),
        "gnomad_mis_z": safe_float(row.get("gnomAD_mis_z")),
        "gnomad_lof_z": safe_float(row.get("gnomAD_lof_z")),
        "loeuf": safe_float(row.get("LOEUF")),  # LoF observed/expected upper (<0.6 = constrained)
        # Domain annotations (for PM1)
        "interpro_domain": safe_str(row.get("Interpro_domain")),
        # SpliceAI scores (for BP7 - synonymous splice impact)
        "spliceai_ag": safe_float(row.get("SpliceAI_pred_DS_AG")),  # Acceptor gain
        "spliceai_al": safe_float(row.get("SpliceAI_pred_DS_AL")),  # Acceptor loss
        "spliceai_dg": safe_float(row.get("SpliceAI_pred_DS_DG")),  # Donor gain
        "spliceai_dl": safe_float(row.get("SpliceAI_pred_DS_DL")),  # Donor loss
        # gnomAD homozygote count (for BS2 - healthy adult observation)
        "gnomad_hom": safe_float(row.get("gnomAD4.1_joint_AC_hom")),
    }


def chunk_dataframe(df: pd.DataFrame, use_grch37: bool = False) -> list[tuple[str, str, dict]]:
    """
    Convert DataFrame to list of (id, text, metadata) tuples.

    Args:
        df: DataFrame with variant data
        use_grch37: Use GRCh37/hg19 coordinates instead of GRCh38

    Returns:
        List of tuples: (variant_id, text_chunk, metadata_dict)
    """
    results = []
    for _, row in df.iterrows():
        # Skip rows with missing GRCh37 coordinates if using GRCh37
        if use_grch37:
            hg19_pos = row.get("hg19_pos(1-based)")
            if pd.isna(hg19_pos) or hg19_pos == ".":
                continue

        var_id = variant_to_id(row, use_grch37=use_grch37)
        text = variant_to_text(row, use_grch37=use_grch37)
        metadata = variant_to_metadata(row, use_grch37=use_grch37)
        results.append((var_id, text, metadata))
    return results

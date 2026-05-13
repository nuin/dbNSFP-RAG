"""Document chunking module - converts variant rows to structured text."""

from typing import Optional

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


def derive_consequence(row: pd.Series) -> str:
    """
    Derive variant consequence from available dbNSFP data.

    Since dbNSFP variant files don't have an explicit consequence column,
    we infer it from:
    1. HGVSp notation (stop_gained, frameshift, etc.)
    2. codon_degeneracy (2 = synonymous)
    3. aa_ref vs aa_alt comparison
    4. MutationTaster prediction (can indicate splice/nonsense)
    """
    # Check HGVSp for explicit consequence indicators
    hgvsp = str(row.get("HGVSp_snpEff", ""))
    hgvsc = str(row.get("HGVSc_snpEff", ""))

    # Frameshift detection
    if "fs" in hgvsp.lower() or "frameshift" in hgvsp.lower():
        return "frameshift_variant"

    # Get amino acid ref/alt early for multiple checks
    aa_ref = str(row.get("aaref", ""))
    aa_alt = str(row.get("aaalt", ""))

    # Stop gained (nonsense) - check aa_alt first as it's most reliable
    if aa_alt in ("*", "X", "Ter") and aa_ref not in ("*", "X", "Ter"):
        return "stop_gained"

    # Also check HGVSp for stop gained
    if "Ter" in hgvsp or "*" in hgvsp:
        if aa_ref not in ("*", "X", "Ter"):
            return "stop_gained"

    # Stop lost
    if aa_ref in ("*", "X", "Ter") and aa_alt not in ("*", "X", "Ter", ""):
        return "stop_lost"

    # Start lost
    if aa_ref == "M" and row.get("aapos") == 1:
        if aa_alt and aa_alt != "M":
            return "start_lost"

    # Synonymous detection via codon_degeneracy
    codon_deg = row.get("codon_degeneracy")
    if not pd.isna(codon_deg):
        try:
            # codon_degeneracy: 0=non-degenerate, 2=2-fold, 4=4-fold degenerate
            # If position is degenerate and aa doesn't change, it's synonymous
            deg = int(str(codon_deg).split(";")[0])
            if deg > 0 and aa_ref == aa_alt and aa_ref:
                return "synonymous_variant"
        except (ValueError, TypeError):
            pass

    # Synonymous - aa_ref equals aa_alt
    if aa_ref and aa_alt and aa_ref == aa_alt:
        return "synonymous_variant"

    # Missense - different amino acids
    if aa_ref and aa_alt and aa_ref != aa_alt:
        # Check for in-frame insertion/deletion
        ref_nt = str(row.get("ref", ""))
        alt_nt = str(row.get("alt", ""))
        if len(ref_nt) != len(alt_nt):
            len_diff = abs(len(ref_nt) - len(alt_nt))
            if len_diff % 3 == 0:
                if len(ref_nt) > len(alt_nt):
                    return "inframe_deletion"
                else:
                    return "inframe_insertion"
            else:
                return "frameshift_variant"
        return "missense_variant"

    # Splice site detection from HGVSc
    if hgvsc:
        # Check for splice site notation (e.g., c.1234+1G>A, c.1234-2T>C)
        import re
        splice_match = re.search(r'[+-][12][ACGT]>', hgvsc)
        if splice_match:
            if "+1" in hgvsc or "+2" in hgvsc:
                return "splice_donor_variant"
            if "-1" in hgvsc or "-2" in hgvsc:
                return "splice_acceptor_variant"

    # Check MutationTaster for splice predictions
    mt_pred = str(row.get("MutationTaster_pred", "")).lower()
    if "splice" in mt_pred:
        return "splice_region_variant"

    # Default: if we have ref/alt nucleotides but no aa change info
    if not aa_ref and not aa_alt:
        ref_nt = str(row.get("ref", ""))
        alt_nt = str(row.get("alt", ""))
        if ref_nt and alt_nt and len(ref_nt) == 1 and len(alt_nt) == 1:
            return "SNV"  # Generic single nucleotide variant

    return ""  # Unknown


def variant_to_metadata(
    row: pd.Series, use_grch37: bool = False, gene_data: Optional[dict] = None
) -> dict:
    """
    Extract metadata for filtering/retrieval.

    Args:
        row: Variant data row
        use_grch37: Use GRCh37/hg19 coordinates
        gene_data: Pre-loaded gene constraint data (from gene_data.py)
    """

    def safe_str(val):
        return str(val) if not pd.isna(val) else ""

    def safe_float(val):
        if pd.isna(val):
            return None
        try:
            # Handle multi-value fields
            if isinstance(val, str):
                parts = val.split(";")
                for p in parts:
                    p = p.strip()
                    if p and p != ".":
                        return float(p)
                return None
            return float(val)
        except (ValueError, TypeError):
            return None

    if use_grch37:
        chrom = safe_str(row.get("hg19_chr", row.get("#chr", ""))).replace("chr", "")
        pos_val = row.get("hg19_pos(1-based)", row.get("pos(1-based)", 0))
    else:
        chrom = safe_str(row.get("#chr", "")).replace("chr", "")
        pos_val = row.get("pos(1-based)", 0)

    # Always preserve hg38 coords so cross-build annotations (e.g. gnomAD v4.1,
    # SpliceAI) can be applied later without re-touching dbNSFP.
    hg38_chr_val = safe_str(row.get("#chr", "")).replace("chr", "") or None
    hg38_pos_raw = row.get("pos(1-based)")
    if pd.isna(hg38_pos_raw):
        hg38_pos_val = None
    else:
        try:
            hg38_pos_val = int(float(hg38_pos_raw))
        except (TypeError, ValueError):
            hg38_pos_val = None

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

    # Get gene constraint scores from gene_data if available
    gene_name = safe_str(row.get("genename"))
    gnomad_pli = None
    loeuf = None
    gnomad_mis_oe = None

    if gene_data and gene_name and gene_name in gene_data:
        gene_info = gene_data[gene_name]
        gnomad_pli = gene_info.get("gnomad_pli")
        loeuf = gene_info.get("loeuf")
        gnomad_mis_oe = gene_info.get("gnomad_mis_oe")

    # Derive consequence from available data
    consequence = derive_consequence(row)

    return {
        "chr": chrom,
        "pos": int(pos_val) if not pd.isna(pos_val) else 0,
        "ref": safe_str(row.get("ref")),
        "alt": safe_str(row.get("alt")),
        "gene": gene_name,
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
        # Consequence/effect type (derived from HGVSp/codon data)
        "consequence": consequence,
        # Gene constraint scores (from gene file, not variant file)
        "gnomad_pli": gnomad_pli,
        "gnomad_mis_oe": gnomad_mis_oe,  # Missense observed/expected
        "loeuf": loeuf,  # LoF observed/expected upper (<0.35 = highly constrained)
        # Domain annotations (for PM1)
        "interpro_domain": safe_str(row.get("Interpro_domain")),
        # GRCh38 coords always preserved (for cross-build annotation lookups)
        "hg38_chr": hg38_chr_val,
        "hg38_pos": hg38_pos_val,
        # gnomAD homozygote count (for BS2 - healthy adult observation)
        "gnomad_hom": safe_float(row.get("gnomAD4.1_joint_nhomalt")),
    }


def chunk_dataframe(
    df: pd.DataFrame, use_grch37: bool = False, gene_data: Optional[dict] = None
) -> list[tuple[str, str, dict]]:
    """
    Convert DataFrame to list of (id, text, metadata) tuples.

    Args:
        df: DataFrame with variant data
        use_grch37: Use GRCh37/hg19 coordinates instead of GRCh38
        gene_data: Pre-loaded gene constraint data (from gene_data.py)

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
        metadata = variant_to_metadata(row, use_grch37=use_grch37, gene_data=gene_data)
        results.append((var_id, text, metadata))
    return results

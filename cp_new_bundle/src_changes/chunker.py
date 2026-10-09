"""Document chunking module - converts variant rows to structured text."""

import pandas as pd

from .variant_text import (
    _parse_first_float as _parse_first_float,
)
from .variant_text import (
    format_prediction as format_prediction,
)
from .variant_text import (
    format_score as format_score,
)
from .variant_text import (
    interpret_cadd as interpret_cadd,
)
from .variant_text import (
    interpret_gnomad_af as interpret_gnomad_af,
)
from .variant_text import (
    interpret_revel as interpret_revel,
)
from .variant_text import (
    variant_to_text,
)


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
    hgvsp = str(row.get("HGVSp_snpEff", ""))
    hgvsc = str(row.get("HGVSc_snpEff", ""))

    # Frameshift detection
    if "fs" in hgvsp.lower() or "frameshift" in hgvsp.lower():
        return "frameshift_variant"

    aa_ref = str(row.get("aaref", ""))
    aa_alt = str(row.get("aaalt", ""))

    # Stop gained (nonsense) - check aa_alt first as it's most reliable
    if aa_alt in ("*", "X", "Ter") and aa_ref not in ("*", "X", "Ter"):
        return "stop_gained"

    # Also check HGVSp for stop gained
    if ("Ter" in hgvsp or "*" in hgvsp) and aa_ref not in ("*", "X", "Ter"):
        return "stop_gained"

    if aa_ref in ("*", "X", "Ter") and aa_alt not in ("*", "X", "Ter", ""):
        return "stop_lost"

    if aa_ref == "M" and row.get("aapos") == 1 and aa_alt and aa_alt != "M":
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

        splice_match = re.search(r"[+-][12][ACGT]>", hgvsc)
        if splice_match:
            if "+1" in hgvsc or "+2" in hgvsc:
                return "splice_donor_variant"
            if "-1" in hgvsc or "-2" in hgvsc:
                return "splice_acceptor_variant"

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
    row: pd.Series, use_grch37: bool = False, gene_data: dict | None = None
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
    df: pd.DataFrame, use_grch37: bool = False, gene_data: dict | None = None
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

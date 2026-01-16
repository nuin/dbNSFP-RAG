"""Generate clickable links to external databases for evidence verification.

This module creates URLs for major genomics databases so users can verify
the data sources used in ACMG variant classification.
"""

from typing import Optional
from urllib.parse import quote


def clinvar_url(clinvar_id: str) -> Optional[str]:
    """Generate link to ClinVar variant page.

    Args:
        clinvar_id: ClinVar variation ID (numeric, RCV, or VCV format)

    Returns:
        URL to ClinVar variant page or None if ID is invalid
    """
    if not clinvar_id or clinvar_id in ("", ".", "N/A"):
        return None
    # Handle RCV/VCV accession formats
    if clinvar_id.startswith(("RCV", "VCV")):
        return f"https://www.ncbi.nlm.nih.gov/clinvar/{clinvar_id}"
    # Numeric variation ID
    return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}"


def gnomad_variant_url(
    chrom: str, pos: int, ref: str, alt: str, build: str = "grch37"
) -> str:
    """Generate link to gnomAD variant page.

    Args:
        chrom: Chromosome (with or without 'chr' prefix)
        pos: Genomic position
        ref: Reference allele
        alt: Alternate allele
        build: Genome build ('grch37' or 'grch38')

    Returns:
        URL to gnomAD variant page
    """
    dataset = "gnomad_r2_1" if build == "grch37" else "gnomad_r4"
    chr_clean = str(chrom).replace("chr", "")
    return (
        f"https://gnomad.broadinstitute.org/variant/"
        f"{chr_clean}-{pos}-{ref}-{alt}?dataset={dataset}"
    )


def gnomad_gene_url(gene: str) -> str:
    """Generate link to gnomAD gene constraint page.

    Args:
        gene: Gene symbol (e.g., 'BRCA1')

    Returns:
        URL to gnomAD gene page with constraint metrics
    """
    return f"https://gnomad.broadinstitute.org/gene/{quote(gene)}?dataset=gnomad_r4"


def uniprot_url(gene: str) -> str:
    """Generate link to UniProt search for gene.

    Args:
        gene: Gene symbol

    Returns:
        URL to UniProt search results for human protein
    """
    return (
        f"https://www.uniprot.org/uniprotkb?"
        f"query={quote(gene)}+AND+organism_id:9606"
    )


def interpro_url(domain: str) -> Optional[str]:
    """Generate link to InterPro domain search.

    Args:
        domain: InterPro domain name or ID

    Returns:
        URL to InterPro search or None if domain is empty
    """
    if not domain or domain in ("", ".", "N/A"):
        return None
    return f"https://www.ebi.ac.uk/interpro/search/text/{quote(domain)}"


def ensembl_transcript_url(transcript: str) -> Optional[str]:
    """Generate link to Ensembl transcript page.

    Args:
        transcript: Ensembl transcript ID (e.g., 'ENST00000302118')

    Returns:
        URL to Ensembl transcript page or None if ID is empty
    """
    if not transcript or transcript in ("", ".", "N/A"):
        return None
    # Handle multi-value fields (take first transcript)
    if ";" in transcript:
        transcript = transcript.split(";")[0].strip()
    return f"https://ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript}"


def omim_url(gene: str) -> str:
    """Generate link to OMIM gene search.

    Args:
        gene: Gene symbol

    Returns:
        URL to OMIM search results
    """
    return f"https://omim.org/search?search={quote(gene)}"


def pubmed_functional_url(gene: str, variant_type: str = "pathogenic") -> str:
    """Generate link to PubMed search for functional studies.

    Args:
        gene: Gene symbol
        variant_type: Type of variant for search context

    Returns:
        URL to PubMed search results
    """
    query = f"{gene} {variant_type} functional study"
    return f"https://pubmed.ncbi.nlm.nih.gov/?term={quote(query)}"


def generate_evidence_links(metadata: dict, build: str = "grch37") -> dict:
    """Generate all available external links for a variant.

    Args:
        metadata: Variant metadata dictionary containing chr, pos, ref, alt,
                  gene, transcript, clinvar_id, interpro_domain, etc.
        build: Genome build for gnomAD link ('grch37' or 'grch38')

    Returns:
        Dictionary of link names to URLs (only non-None links included)
    """
    links = {}

    # ClinVar - direct variant lookup
    clinvar_id = metadata.get("clinvar_id")
    if clinvar_id:
        url = clinvar_url(clinvar_id)
        if url:
            links["clinvar"] = url

    # gnomAD variant page
    chrom = metadata.get("chr")
    pos = metadata.get("pos")
    ref = metadata.get("ref")
    alt = metadata.get("alt")
    if all([chrom, pos, ref, alt]):
        links["gnomad_variant"] = gnomad_variant_url(chrom, pos, ref, alt, build)

    # Gene-based links
    gene = metadata.get("gene")
    if gene and gene not in ("", ".", "N/A"):
        links["gnomad_gene"] = gnomad_gene_url(gene)
        links["uniprot"] = uniprot_url(gene)
        links["omim"] = omim_url(gene)
        links["pubmed"] = pubmed_functional_url(gene)

    # InterPro domain
    domain = metadata.get("interpro_domain")
    if domain:
        url = interpro_url(domain)
        if url:
            links["interpro"] = url

    # Ensembl transcript
    transcript = metadata.get("transcript")
    if transcript:
        url = ensembl_transcript_url(transcript)
        if url:
            links["ensembl"] = url

    return links


# Link display names for UI
LINK_DISPLAY_NAMES = {
    "clinvar": "ClinVar",
    "gnomad_variant": "gnomAD",
    "gnomad_gene": "Gene Constraint",
    "uniprot": "UniProt",
    "omim": "OMIM",
    "pubmed": "PubMed",
    "interpro": "InterPro",
    "ensembl": "Ensembl",
}

"""ClinVar API client for PS1/PM5 ACMG criteria evaluation.

This module provides a client for querying the NCBI E-utilities API to find
pathogenic variants at the same position (for PS1 - same AA change) or
different AA changes at the same position (for PM5).

API Documentation:
- https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- https://www.ncbi.nlm.nih.gov/books/NBK25500/ (E-utilities)
"""

import re
import time
from functools import lru_cache

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RATE_LIMIT_DELAY = 0.34  # 3 requests/sec without API key


class ClinVarResponseError(ValueError):
    """The upstream service returned data that does not match its response schema."""


def _response_object(value, name: str) -> dict:
    if not isinstance(value, dict):
        raise ClinVarResponseError(f"ClinVar {name} must be an object")
    return value


def _response_strings(value, name: str) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ClinVarResponseError(f"ClinVar {name} must be an array of strings")
    return value


def _response_text(value, name: str) -> str:
    if not isinstance(value, str):
        raise ClinVarResponseError(f"ClinVar {name} must be a string")
    return value


def _parse_variant_summary(uid: str, payload) -> dict:
    variant = _response_object(payload, f"variant {uid}")
    title = _response_text(variant.get("title", ""), "title")
    protein_change = _response_text(variant.get("protein_change") or "", "protein_change")
    if not protein_change:
        match = re.search(r"\(p\.([A-Za-z]+\d+[A-Za-z]+)\)", title)
        if match:
            protein_change = f"p.{match.group(1)}"
    significance = variant.get("clinical_significance")
    description = ""
    if significance is not None:
        significance = _response_object(significance, "clinical_significance")
        description = _response_text(significance.get("description", ""), "description")
    genes = variant.get("genes", [])
    if genes is None:
        genes = []
    if not isinstance(genes, list):
        raise ClinVarResponseError("ClinVar genes must be an array")
    gene_symbols = [
        _response_text(_response_object(gene, "gene").get("symbol", ""), "gene symbol")
        for gene in genes
    ]
    return {
        "uid": uid,
        "title": title,
        "clinical_significance": description,
        "protein_change": protein_change,
        "gene": gene_symbols[0] if gene_symbols else "",
    }


class ClinVarClient:
    """Client for querying ClinVar via NCBI E-utilities."""

    def __init__(self, email: str, api_key: str | None = None):
        """Initialize client with required email and optional API key.

        Args:
            email: Required by NCBI for tracking
            api_key: Optional API key for higher rate limits (10/sec vs 3/sec)
        """
        self.email = email
        self.api_key = api_key
        self._last_request = 0.0
        self.get_pathogenic_at_position = lru_cache(maxsize=10000)(self.get_pathogenic_at_position)
        self.get_pathogenic_in_gene = lru_cache(maxsize=1000)(self.get_pathogenic_in_gene)

    def _rate_limit(self) -> None:
        """Enforce rate limiting to comply with NCBI guidelines."""
        delay = RATE_LIMIT_DELAY if not self.api_key else 0.1
        elapsed = time.time() - self._last_request
        if elapsed < delay:
            time.sleep(delay - elapsed)
        self._last_request = time.time()

    def _make_request(self, endpoint: str, params: dict, timeout: int = 10) -> dict | None:
        """Make a rate-limited request to E-utilities."""
        self._rate_limit()

        params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            resp = requests.get(f"{EUTILS_BASE}/{endpoint}", params=params, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException:
            return None

    def get_pathogenic_at_position(self, chrom: str, pos: int, build: str = "37") -> list[dict]:
        """Find pathogenic variants at a genomic position.

        Used for PS1 (same AA change) and PM5 (same position, different AA).

        Args:
            chrom: Chromosome (e.g., "15" or "X")
            pos: Genomic position
            build: Genome build ("37" or "38")

        Returns:
            List of pathogenic variant records
        """
        # Normalize chromosome
        chrom = str(chrom).replace("chr", "")

        chrpos_field = f"chrpos{build}"
        query = f"{chrom}[chr] AND {pos}:{pos}[{chrpos_field}] AND pathogenic[clinsig]"

        params = {
            "db": "clinvar",
            "term": query,
            "retmode": "json",
            "retmax": 100,
        }

        ids = self._search_variant_ids(params)
        return self._get_variant_summaries(ids) if ids else []

    def get_pathogenic_in_gene(self, gene: str) -> list[dict]:
        """Get all pathogenic variants in a gene for PS1/PM5 lookup.

        Args:
            gene: Gene symbol (e.g., "FBN1")

        Returns:
            List of pathogenic variant records with protein changes
        """
        params = {
            "db": "clinvar",
            "term": f"{gene}[gene] AND pathogenic[clinsig]",
            "retmode": "json",
            "retmax": 1000,
        }

        ids = self._search_variant_ids(params)
        return self._get_variant_summaries(ids) if ids else []

    def _search_variant_ids(self, params: dict) -> list[str]:
        data = self._make_request("esearch.fcgi", params)
        if data is None:
            return []
        data = _response_object(data, "search response")
        if "esearchresult" not in data:
            return []
        search_result = _response_object(data["esearchresult"], "esearchresult")
        return _response_strings(search_result.get("idlist", []), "idlist")

    def _get_variant_summaries(self, ids: list[str]) -> list[dict]:
        """Fetch variant summaries for a list of ClinVar IDs.

        Args:
            ids: List of ClinVar variation IDs

        Returns:
            List of variant summary dictionaries
        """
        # Process in batches of 100
        all_variants = []
        for i in range(0, len(ids), 100):
            batch_ids = ids[i : i + 100]
            params = {
                "db": "clinvar",
                "id": ",".join(batch_ids),
                "retmode": "json",
            }

            data = self._make_request("esummary.fcgi", params)
            if data is None:
                continue
            data = _response_object(data, "summary response")
            results = _response_object(data.get("result", {}), "result")
            for uid in _response_strings(results.get("uids", []), "uids"):
                all_variants.append(_parse_variant_summary(uid, results.get(uid, {})))

        return all_variants


def parse_protein_change(prot_str: str) -> tuple[str, int, str] | None:
    """Parse a protein change string like 'p.Cys347Ter' or 'Cys347Ter'.

    Args:
        prot_str: Protein change notation

    Returns:
        Tuple of (ref_aa, position, alt_aa) or None if parsing fails
    """
    # Remove 'p.' prefix if present
    prot_str = prot_str.replace("p.", "")

    # Match pattern: letters + numbers + letters (e.g., Cys347Ter)
    match = re.match(r"([A-Za-z]+)(\d+)([A-Za-z]+)", prot_str)
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None


def evaluate_ps1_pm5(
    client: ClinVarClient,
    gene: str,
    aa_pos: int,
    aa_ref: str,
    aa_alt: str,
) -> tuple[bool, bool, str]:
    """Evaluate PS1 and PM5 criteria using ClinVar API.

    PS1: Same amino acid change as an established pathogenic variant
    PM5: Different amino acid change at the same position as established pathogenic

    Args:
        client: ClinVarClient instance
        gene: Gene symbol
        aa_pos: Amino acid position
        aa_ref: Reference amino acid (3-letter or 1-letter code)
        aa_alt: Alternate amino acid (3-letter or 1-letter code)

    Returns:
        Tuple of (ps1_met, pm5_met, evidence_string)
    """
    if not gene or not aa_pos or not aa_ref or not aa_alt:
        return False, False, "Insufficient data for PS1/PM5 evaluation"

    pathogenic_variants = client.get_pathogenic_in_gene(gene)

    ps1_met = False
    pm5_met = False
    evidence = []

    # Normalize input (handle both 3-letter and 1-letter codes)
    query_ref = aa_ref.upper()
    query_alt = aa_alt.upper()

    for var in pathogenic_variants:
        prot_change = var.get("protein_change", "")
        parsed = parse_protein_change(prot_change)
        if not parsed:
            continue

        var_ref, var_pos, var_alt = parsed

        # Check if same position
        if var_pos == aa_pos:
            # Normalize variant AA codes for comparison
            var_ref_upper = var_ref.upper()
            var_alt_upper = var_alt.upper()

            # Check for exact match (PS1)
            if var_ref_upper == query_ref and var_alt_upper == query_alt:
                ps1_met = True
                evidence.append(
                    f"PS1: Same AA change ({prot_change}) found pathogenic "
                    f"in ClinVar ({var['clinical_significance']})"
                )
            # Check for same position, different change (PM5)
            elif var_ref_upper == query_ref and var_alt_upper != query_alt:
                pm5_met = True
                evidence.append(
                    f"PM5: Different AA change at same position ({prot_change}) "
                    f"found pathogenic in ClinVar"
                )

    if not evidence:
        return False, False, "No pathogenic variants found at this position in ClinVar"

    return ps1_met, pm5_met, "; ".join(evidence)

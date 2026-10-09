"""Variant identities, local annotations, and cached public API access."""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "exports" / "cp_new"
PER_GENE_DIR = DATA
FINAL = DATA / "seqnext" / "cp_new_seqnext_FINAL.tsv"
VV_CACHE_PATH = DATA / ".vv_cache.json"
CLINVAR_BULK = Path.home() / "data" / "clinvar" / "variant_summary.txt.gz"
HG19_FA = Path.home() / "data" / "hg19" / "hg19.fa"
HG38_FA = Path.home() / "data" / "hg38" / "hg38.fa"


# =============================================================================
# Variant normalization
# =============================================================================

COORD_RE = re.compile(r"^(?:chr)?([\dXYM]+)[:-](\d+)[:-]([ACGTN-]+)[:-]([ACGTN-]+)$", re.IGNORECASE)
RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
HGVS_RE = re.compile(r"^(NM_\d+\.?\d*|ENST\d+\.?\d*)?:?(c\.[A-Za-z\d+\-*>\.]+)$")


_HTTP_CACHE_DIR = REPO / "data" / "exports" / "cp_new" / ".inspector_cache"
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json", "User-Agent": "cp_new-inspector/0.1"})


@dataclass
class Variant:
    """Canonicalized variant identity. Any field may be filled lazily."""

    gene: str | None = None
    transcript: str | None = None
    hgvs_c: str | None = None
    chr_hg19: str | None = None
    pos_hg19: int | None = None
    chr_hg38: str | None = None
    pos_hg38: int | None = None
    ref: str | None = None
    alt: str | None = None
    rsid: str | None = None
    build_input: str | None = None  # which build the user's input was in

    def key(self) -> str:
        if self.chr_hg19 and self.pos_hg19 and self.ref and self.alt:
            return f"{self.chr_hg19}-{self.pos_hg19}-{self.ref}-{self.alt}"
        if self.chr_hg38 and self.pos_hg38 and self.ref and self.alt:
            return f"hg38:{self.chr_hg38}-{self.pos_hg38}-{self.ref}-{self.alt}"
        return f"{self.gene or '?'}|{self.transcript or '?'}|{self.hgvs_c or '?'}"

    def coord_present(self) -> bool:
        return bool((self.chr_hg19 and self.pos_hg19) or (self.chr_hg38 and self.pos_hg38))


def parse_input(args: argparse.Namespace) -> Variant:
    """Normalize argparse args into a Variant. Multiple forms can be combined
    (e.g. --gene GOT2 --coord 16:X:G:A gives both hints)."""
    v = Variant()
    if args.gene:
        v.gene = args.gene.upper()
    if args.transcript:
        v.transcript = args.transcript
    if args.hgvs:
        v.hgvs_c = args.hgvs.strip()
    if args.coord:
        m = COORD_RE.match(args.coord.strip())
        if not m:
            sys.exit(f"--coord must look like chr:pos:ref:alt (got {args.coord!r})")
        chrom, pos, ref, alt = m.groups()
        if args.build == "hg38":
            v.chr_hg38, v.pos_hg38 = chrom, int(pos)
        else:
            v.chr_hg19, v.pos_hg19 = chrom, int(pos)
        v.ref, v.alt = ref.upper(), alt.upper()
        v.build_input = args.build
    if args.rsid:
        if not RSID_RE.match(args.rsid):
            sys.exit(f"--rsid must look like rsNNNN (got {args.rsid!r})")
        v.rsid = args.rsid.lower()
    if not any([v.gene, v.coord_present(), v.rsid, v.hgvs_c]):
        sys.exit("Need one of: --gene + --hgvs/--coord, --coord CHR:POS:REF:ALT, or --rsid rsX")
    return v


def load_vv_cache() -> dict:
    if not VV_CACHE_PATH.exists():
        return {}
    with open(VV_CACHE_PATH) as f:
        return json.load(f)


def find_in_per_gene(v: Variant) -> dict | None:
    """Look up variant in per-gene dbNSFP+gnomAD TSV."""
    if not v.gene:
        return None
    tsv = PER_GENE_DIR / f"{v.gene}.tsv"
    if not tsv.exists():
        return None
    import pandas as pd

    df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".", ""], low_memory=False)
    if v.chr_hg19 and v.pos_hg19:
        mask = (
            (df["hg19_chr"].astype(str) == str(v.chr_hg19))
            & (df["hg19_pos(1-based)"].astype(str) == str(v.pos_hg19))
            & (df["ref"] == v.ref)
            & (df["alt"] == v.alt)
        )
    elif v.chr_hg38 and v.pos_hg38:
        mask = (
            (df["#chr"].astype(str) == str(v.chr_hg38))
            & (df["pos(1-based)"].astype(str) == str(v.pos_hg38))
            & (df["ref"] == v.ref)
            & (df["alt"] == v.alt)
        )
    else:
        return None
    sub = df[mask]
    if sub.empty:
        return None
    return sub.iloc[0].to_dict()


def find_in_classifications(v: Variant) -> list[dict]:
    """Look up variant in cp_new_seqnext_FINAL.tsv. Tries successively looser
    matches: (gene, hgvs_c) -> (chr, pos, ref, alt) -> (gene, chr, pos).
    Returns the union of matches across strategies."""
    if not FINAL.exists():
        return []
    import pandas as pd

    df = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    masks = []
    if v.gene and v.hgvs_c:
        masks.append((df["gene"] == v.gene) & (df["hgvs_c"] == v.hgvs_c))
    if v.chr_hg19 and v.pos_hg19 and v.ref and v.alt:
        masks.append(
            (df["chr_grch37"] == str(v.chr_hg19))
            & (df["pos_grch37"] == str(v.pos_hg19))
            & (df["ref"] == v.ref)
            & (df["alt"] == v.alt)
        )
    if v.gene and v.chr_hg19 and v.pos_hg19:
        # Loose: same gene + position (catches ClinVar rows with empty ref/alt)
        masks.append(
            (df["gene"] == v.gene)
            & (df["chr_grch37"] == str(v.chr_hg19))
            & (df["pos_grch37"] == str(v.pos_hg19))
        )
    if not masks:
        return []
    combined = masks[0]
    for m in masks[1:]:
        combined = combined | m
    return (
        df[combined]
        .drop_duplicates(subset=["gene", "transcript", "hgvs_c", "source"])
        .to_dict("records")
    )


def _cache_path(key: str) -> Path:
    safe = re.sub(r"[^a-zA-Z0-9._-]", "_", key)[:200]
    return _HTTP_CACHE_DIR / f"{safe}.json"


def fetch_variantvalidator(v: Variant, transcript: str | None = None) -> dict | None:
    """Resolve HGVS on a transcript via VariantValidator REST."""
    if not (v.chr_hg19 and v.pos_hg19 and v.ref and v.alt):
        return None
    tx = transcript or v.transcript or "select"
    url = (
        f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/"
        f"GRCh37/{v.chr_hg19}-{v.pos_hg19}-{v.ref}-{v.alt}/{tx}"
    )
    return http_get_json(url, f"vv_{v.key()}_{tx}", sleep_after=0.3)


def fetch_myvariant(v: Variant) -> dict | None:
    """myvariant.info aggregator. Build hg19 query string from coords or HGVS."""
    if v.chr_hg19 and v.pos_hg19 and v.ref and v.alt:
        q = f"chr{v.chr_hg19}:g.{v.pos_hg19}{v.ref}>{v.alt}"
    elif v.rsid:
        q = v.rsid
    else:
        return None
    url = f"https://myvariant.info/v1/variant/{q}?assembly=hg19"
    return http_get_json(url, f"mv_{q}", sleep_after=0.1)


def fetch_ensembl_gene(gene_symbol: str) -> dict | None:
    """Ensembl REST lookup for gene metadata in GRCh37."""
    url = f"https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    return http_get_json(url, f"ens37_gene_{gene_symbol}", sleep_after=0.1)


def fetch_clinvar_eutils(v: Variant) -> dict | None:
    """ClinVar variation lookup via NCBI eutils -- gives submitter + assertion detail."""
    if not v.rsid:
        return None
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar"
        f"&term={v.rsid}&retmode=json"
    )
    return http_get_json(url, f"clinvar_eutils_{v.rsid}", sleep_after=0.5)


def http_get_json(url: str, cache_key: str, sleep_after: float = 0.0) -> dict | None:
    """GET a JSON object with on-disk caching; report cache and HTTP failures."""
    cp = _cache_path(cache_key)
    if cp.exists():
        try:
            cached = json.loads(cp.read_text())
            if isinstance(cached, dict):
                return cached
            logger.warning("Ignoring non-object API cache %s", cp)
        except (OSError, ValueError) as exc:
            logger.warning("Cannot read API cache %s: %s", cp, exc)
    try:
        response = SESSION.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        logger.warning("API request failed for %s: %s", url, exc)
        return None
    if not isinstance(data, dict):
        logger.warning("Expected JSON object from %s", url)
        return None
    try:
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps(data))
    except OSError as exc:
        logger.warning("Cannot write API cache %s: %s", cp, exc)
    if sleep_after:
        time.sleep(sleep_after)
    return data

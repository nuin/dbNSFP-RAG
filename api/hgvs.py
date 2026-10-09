"""Parse genomic HGVS coordinates accepted by the API."""


def parse_hgvs_g(hgvs: str) -> tuple[str, int, str, str]:
    """Parse HGVS genomic notation (e.g., NC_000017.10:g.41197801T>A or 17:g.41197801T>A)."""
    import re

    # Pattern: optional NC_XXXXXX.X: or chrX: prefix, then g.POS REF>ALT
    pattern = r"(?:NC_0+(\d+)(?:\.\d+)?:)?(?:chr)?(\d+|X|Y)?:?g\.(\d+)([ACGT]+)>([ACGT]+)"
    match = re.match(pattern, hgvs.strip(), re.IGNORECASE)

    if not match:
        raise ValueError(f"Invalid HGVS format: {hgvs}")

    nc_chr, simple_chr, pos, ref, alt = match.groups()
    chrom = nc_chr or simple_chr
    if not chrom:
        raise ValueError(f"Could not parse chromosome from: {hgvs}")

    return chrom, int(pos), ref.upper(), alt.upper()

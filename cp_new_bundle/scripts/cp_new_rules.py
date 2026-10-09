"""BED transcript selection and cp_new classification rules."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

INTRONIC_ROI_NEG = -15
INTRONIC_ROI_POS = 6
INTRONIC_OFFSET_RE = re.compile(r"c\.(?:\*?-?)?\d+([+\-])(\d+)")


def load_bed(
    bed_path: Path,
) -> tuple[
    dict[tuple[str, str], list[tuple[int, int, str]]],
    dict[str, str],
]:
    """Return (region_index, gene_default_transcript).

    region_index: {(chr_no_prefix, gene): [(start, end, transcript), ...]}
        BED 0-based half-open compared with dbNSFP 1-based pos as start < pos <= end.

    gene_default_transcript: {gene: most-common transcript across the gene's
        regions}. Used as a fallback when a variant's coordinates don't fall
        in any BED region but the gene IS in the BED.
    """
    region_idx: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    tx_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0].removeprefix("chr")
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            gene = parts[4]
            transcript = parts[5]
            region_idx[(chrom, gene)].append((start, end, transcript))
            tx_counts[gene][transcript] += 1
    for regions in region_idx.values():
        regions.sort()
    gene_default = {g: max(txs.items(), key=lambda kv: kv[1])[0] for g, txs in tx_counts.items()}
    return region_idx, gene_default


def lookup_bed_transcript(
    region_idx: dict[tuple[str, str], list[tuple[int, int, str]]],
    gene_default: dict[str, str],
    chrom: str,
    pos: int,
    gene: str,
) -> str | None:
    """Resolve the transcript for this variant:

    1. If the variant's hg19 position falls in a BED region for this gene,
       return that region's transcript (handles APC E01 / RAD51D alt-E3).
    2. Else if the gene is in the BED at all, return the gene's primary
       (most common) transcript -- so variants outside the assay's coverage
       windows still get a valid transcript anchor.
    3. Else (gene not in BED), return None.
    """
    regions = region_idx.get((chrom, gene))
    if regions:
        for start, end, tx in regions:
            if start < pos <= end:
                return tx
    return gene_default.get(gene)


def safe_float(v) -> float | None:
    if v is None or v == "" or v == "." or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def is_synonymous(row) -> bool:
    """Synonymous if HGVSp shows same AA (p.X=) or aaref == aaalt."""
    hgvsp = str(row.get("HGVSp_snpEff", ""))
    if "p.=" in hgvsp or hgvsp.endswith("="):
        return True
    aaref = str(row.get("aaref", "")).strip().upper()
    aaalt = str(row.get("aaalt", "")).strip().upper()
    if aaref and aaref == aaalt and aaref != "X":
        return True
    # codon_degeneracy: 2 = synonymous in dbNSFP convention (4-fold degenerate sites)
    cd = str(row.get("codon_degeneracy", "")).strip()
    return cd in ("2", "4")


def parse_intronic_offset(hgvsc: str) -> int | None:
    """Return the +/- offset from HGVSc, or None if it's a coding (non-intronic) c.

    'c.1207-1G>A'  -> -1
    'c.500+6A>T'   -> 6
    'c.500A>G'     -> None  (coding, no offset)
    'c.500-100G>T' -> -100  (deep intronic)
    """
    if not hgvsc:
        return None
    m = INTRONIC_OFFSET_RE.search(hgvsc)
    if not m:
        return None
    sign, mag = m.group(1), int(m.group(2))
    return -mag if sign == "-" else mag


def is_intronic(row, hgvsc_resolved: str | None = None) -> bool:
    """Intronic if HGVSc has a +/- offset, regardless of HGVSp.

    Prefers the VV-resolved hgvs_c when available (anchored to BED's NM_)
    over dbNSFP's snpEff value (which may be on a different transcript).
    """
    src = hgvsc_resolved or str(row.get("HGVSc_snpEff", ""))
    return parse_intronic_offset(src) is not None


def in_splice_roi(offset: int) -> bool:
    """True if the intronic offset is within the lab's reportable region
    (-15 to +6 around the exon boundary)."""
    return INTRONIC_ROI_NEG <= offset < 0 or 0 < offset <= INTRONIC_ROI_POS


def classify(row, hgvsc_resolved: str | None = None) -> tuple[str, str] | None:
    """Return (classification, rule) or None if no workflow matches.

    Args:
        row: variant row from per-gene TSV.
        hgvsc_resolved: VariantValidator-resolved HGVSc on the BED's NM_,
            if known. Used for accurate intronic detection. If None, falls
            back to dbNSFP's snpEff HGVSc (may be on a different transcript).

    SpliceAI handling: a missing spliceai_ds_max_masked is treated as 0
    (no detected splice signal). SpliceAI v1.3 silently skips variants
    outside its GENCODE V24 canonical gene model.

    Intronic ROI: variants whose intronic offset falls outside [-15, +6]
    are not eligible for W2 -- deep intronic positions need separate
    splice/regulatory evaluation, not benign-synonymous auto-call.
    """
    faf = safe_float(row.get("gnomad_v41_faf95_grpmax"))
    revel = safe_float(row.get("REVEL_score"))
    spliceai_raw = safe_float(row.get("spliceai_ds_max_masked"))
    spliceai = spliceai_raw if spliceai_raw is not None else 0.0
    phastcons = safe_float(row.get("phastCons100way_vertebrate"))
    phylop = safe_float(row.get("phyloP100way_vertebrate"))

    # Workflow 1 -- doesn't need SpliceAI
    if faf is not None and faf > 0.05:
        return "Benign", "FAF >5%"

    # Workflow 2 -- synonymous / silent benign
    if is_synonymous(row) and spliceai <= 0.1 and phastcons is not None and phastcons < 1.0:
        # Use resolved c. for intronic detection if available
        offset = parse_intronic_offset(hgvsc_resolved or row.get("HGVSc_snpEff", ""))
        if offset is not None:
            # Intronic: enforce ROI + PhyloP gate (BOTH PhastCons AND PhyloP required)
            if not in_splice_roi(offset):
                return None  # deep intronic -- not eligible for W2 auto-classify
            if phylop is None or phylop >= 0.1:
                return None  # intronic but PhyloP fails or unknown
            return (
                "Benign",
                f"synonymous, SpliceAI<=0.1, PhastCons<1.0, PhyloP<0.1, intronic_ROI({offset:+d})",
            )
        else:
            # Coding-synonymous (no offset)
            return "Benign", "synonymous, SpliceAI<=0.1, PhastCons<1.0"

    # Workflow 3 -- rare LB
    if faf is not None and faf > 0.001 and revel is not None and revel < 0.290 and spliceai <= 0.1:
        return "Likely_benign", "FAF >0.1%, REVEL <0.29, SpliceAI<=0.1"

    return None

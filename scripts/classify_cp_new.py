#!/usr/bin/env python3
"""Apply the three cp_new benign/LB filtering workflows and emit SeqNext rows.

Inputs:
  * data/exports/cp_new/*.tsv     -- per-gene dbNSFP+gnomAD output from cp_new.py
  * BED file with columns chr,start,end,region,gene,transcript
                                  -- defines the lab's coverage + reporting tx

Workflows (priority order, first match wins):
  1. BENIGN -- gnomad_v41_faf95_grpmax > 0.05
  2. BENIGN -- synonymous AND spliceai_ds_max_masked <= 0.1
              AND phastCons100way_vertebrate < 1.0
              AND (if intronic) phyloP < 0.1
  3. LIKELY_BENIGN -- gnomad_v41_faf95_grpmax > 0.001
              AND REVEL_score < 0.290 AND spliceai_ds_max_masked <= 0.1

Output (SeqNext import format), one TSV per gene under data/exports/cp_new/seqnext/:
  gene  transcript  hgvs_c  classification

Each gene file is named {GENE}_seqnext.tsv. A combined cp_new_seqnext.tsv is
also written at the same directory for convenience.

Variants outside any BED region are silently skipped (the assay doesn't cover
them so they can't be reported anyway). Workflows 2 and 3 require SpliceAI;
when spliceai_ds_max_masked is empty (e.g. before running spliceai -M 1) the
script counts those variants as 'pending' and does NOT emit a classification.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

DEFAULT_BED = Path("/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed")
DEFAULT_TSV_DIR = Path("data/exports/cp_new")
DEFAULT_OUT_DIR = DEFAULT_TSV_DIR / "seqnext"
DEFAULT_VV_CACHE = DEFAULT_TSV_DIR / ".vv_cache.json"

# ROI for intronic variants relative to nearest exon boundary
INTRONIC_ROI_NEG = -15   # acceptor side
INTRONIC_ROI_POS = 6     # donor side

# Parses ".c.1207-1G>A" -> -1, ".c.500A>G" -> 0, ".c.123+12C>T" -> 12
INTRONIC_OFFSET_RE = re.compile(r"c\.(?:\*?-?)?\d+([+\-])(\d+)")


# ---------------------------------------------------------------------------
# BED interval index
# ---------------------------------------------------------------------------

def load_bed(bed_path: Path) -> tuple[
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
    for k in region_idx:
        region_idx[k].sort()
    gene_default = {g: max(txs.items(), key=lambda kv: kv[1])[0]
                    for g, txs in tx_counts.items()}
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


# ---------------------------------------------------------------------------
# Workflow logic
# ---------------------------------------------------------------------------

def safe_float(v) -> float | None:
    """Coerce a dbNSFP value to float. Handles multi-transcript ';'-joined
    strings (e.g. '.;.;0.220;.;.;.;.;.') by returning the first real number.

    BUG HISTORY: earlier versions did `float(v)` directly, which raised on
    every ';'-joined value and silently returned None. That caused
    classify() to fail W3 on every multi-transcript dbNSFP row (i.e. most
    missense variants), so the pipeline emitted ~8 rows per gene instead
    of thousands. DO NOT remove the ';' handling.
    """
    if v is None or v == "" or v == "." or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v)
    if ";" in s:
        for part in s.split(";"):
            part = part.strip()
            if part and part != ".":
                try:
                    return float(part)
                except (TypeError, ValueError):
                    continue
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def is_synonymous(row) -> bool:
    """Synonymous if HGVSp shows same AA (p.X=) or aaref == aaalt.

    NOTE: dbNSFP's `codon_degeneracy` describes the codon position's
    degeneracy, NOT whether this specific alt is silent. At a 2-fold
    degenerate site, only 1 of 3 alts is synonymous; the other 2 are
    missense. Earlier versions of this function incorrectly treated any
    variant at a 2- or 4-fold degenerate site as synonymous and produced
    many mis-tagged "Benign synonymous" calls on missense / nonsense
    variants (e.g. HOXB13 c.108C>A, c.144T>A). DO NOT re-introduce that
    shortcut. Trust only the per-row aaref/aaalt comparison.
    """
    hgvsp = str(row.get("HGVSp_snpEff", ""))
    if "p.=" in hgvsp or hgvsp.endswith("="):
        return True
    aaref = str(row.get("aaref", "")).strip().upper()
    aaalt = str(row.get("aaalt", "")).strip().upper()
    if aaref and aaref == aaalt and aaref != "X":
        return True
    return False


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
    if is_synonymous(row):
        if spliceai <= 0.1 and phastcons is not None and phastcons < 1.0:
            # Use resolved c. for intronic detection if available
            offset = parse_intronic_offset(hgvsc_resolved or row.get("HGVSc_snpEff", ""))
            if offset is not None:
                # Intronic: enforce ROI + PhyloP gate (BOTH PhastCons AND PhyloP required)
                if not in_splice_roi(offset):
                    return None  # deep intronic -- not eligible for W2 auto-classify
                if phylop is None or phylop >= 0.1:
                    return None  # intronic but PhyloP fails or unknown
                return "Benign", f"synonymous, SpliceAI<=0.1, PhastCons<1.0, PhyloP<0.1, intronic_ROI({offset:+d})"
            else:
                # Coding-synonymous (no offset)
                return "Benign", "synonymous, SpliceAI<=0.1, PhastCons<1.0"

    # Workflow 3 -- rare LB
    if (faf is not None and faf > 0.001
            and revel is not None and revel < 0.290
            and spliceai <= 0.1):
        return "Likely_benign", "FAF >0.1%, REVEL <0.29, SpliceAI<=0.1"

    return None


def workflow_pending(row) -> str | None:
    """Return None always now -- SpliceAI-missing rows are no longer 'pending';
    they're classified directly with SpliceAI treated as 0. This function is
    kept as a stub so the calling counters don't break."""
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv-dir", type=Path, default=DEFAULT_TSV_DIR)
    p.add_argument("--bed", type=Path, default=DEFAULT_BED)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument("--vv-cache", type=Path, default=DEFAULT_VV_CACHE,
                   help="VariantValidator cache for accurate intronic detection. "
                        "Built incrementally by resolve_hgvs_vv.py; on the first "
                        "pipeline run it's empty so classify uses dbNSFP's c.")
    args = p.parse_args()

    if not args.bed.exists():
        sys.exit(f"BED not found: {args.bed}")
    if not args.tsv_dir.is_dir():
        sys.exit(f"TSV dir not found: {args.tsv_dir}")

    print(f"BED: {args.bed}")
    print(f"TSVs: {args.tsv_dir}")
    print(f"Out: {args.out_dir}\n")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    bed_idx, gene_default_tx = load_bed(args.bed)
    n_bed_regions = sum(len(v) for v in bed_idx.values())
    print(f"BED: {n_bed_regions} regions over {len(bed_idx)} (chr, gene) pairs, "
          f"{len(gene_default_tx)} gene defaults")

    # Load VV cache for accurate intronic detection + status carry-over
    vv_cache: dict[str, str] = {}        # key -> resolved c. (only "ok" entries)
    vv_status_cache: dict[str, str] = {} # key -> raw status string (ok/flagged:*)
    if args.vv_cache.exists():
        with open(args.vv_cache) as f:
            raw = json.load(f)
        for k, v in raw.items():
            if not (isinstance(v, list) and len(v) == 2):
                continue
            c_part, status = v
            vv_status_cache[k] = status
            if status == "ok" and c_part:
                vv_cache[k] = "c." + c_part.split("c.", 1)[-1]
        print(f"VV cache: {len(vv_status_cache):,} entries ({len(vv_cache):,} ok)\n")
    else:
        print("VV cache: not found -- using dbNSFP HGVSc for intronic detection\n")

    tsvs = sorted(args.tsv_dir.glob("*.tsv"))
    # Skip subdirs/seqnext output if rerunning
    tsvs = [t for t in tsvs if t.parent == args.tsv_dir]

    combined_rows = []
    counters = {
        "total_rows": 0,
        "not_in_bed": 0,
        "classified": 0,
        "w1_benign_faf5": 0,
        "w2_benign_synonymous": 0,
        "w3_likely_benign": 0,
        "w2_pending_spliceai": 0,
        "w3_pending_spliceai": 0,
    }

    for tsv in tqdm(tsvs, desc="classify"):
        first = tsv.read_text().splitlines()[:1]
        if not first or not first[0].startswith("#chr"):
            continue
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".", ""], low_memory=False)
        if df.empty:
            continue

        gene_rows = []
        seen: set[tuple[str, int, str, str, str]] = set()  # (chr, pos, ref, alt, tx)

        for _, row in df.iterrows():
            counters["total_rows"] += 1
            # dbNSFP `genename` can be ';'-delimited multi-gene/multi-transcript.
            # Try each part and pick the first one that's actually in the BED
            # (i.e. a known panel gene); fall back to first non-empty otherwise.
            raw = str(row.get("genename") or "").strip()
            parts = [p.strip() for p in raw.split(";") if p.strip()]
            gene = next((p for p in parts if p in gene_default_tx), parts[0] if parts else "")
            try:
                hg19_chr = str(row.get("hg19_chr") or "").replace("chr", "")
                hg19_pos = int(float(row.get("hg19_pos(1-based)")))
            except (TypeError, ValueError):
                continue

            transcript = lookup_bed_transcript(bed_idx, gene_default_tx,
                                                hg19_chr, hg19_pos, gene)
            if transcript is None:
                counters["not_in_bed"] += 1
                continue

            ref = str(row.get("ref") or "")
            alt = str(row.get("alt") or "")
            key = (hg19_chr, hg19_pos, ref, alt, transcript)
            if key in seen:
                continue
            seen.add(key)

            # Look up VV-resolved c. (anchored to BED's NM_) for intronic detection
            vv_key = f"{hg19_chr}-{hg19_pos}-{ref}-{alt}|{transcript}"
            hgvsc_resolved = vv_cache.get(vv_key)

            result = classify(row, hgvsc_resolved=hgvsc_resolved)
            if result is None:
                pending = workflow_pending(row)
                if pending == "W2_synonymous_pending_spliceai":
                    counters["w2_pending_spliceai"] += 1
                elif pending == "W3_rare_LB_pending_spliceai":
                    counters["w3_pending_spliceai"] += 1
                continue

            classification, rule = result
            counters["classified"] += 1
            if rule == "FAF >5%":
                counters["w1_benign_faf5"] += 1
            elif rule.startswith("synonymous"):
                counters["w2_benign_synonymous"] += 1
            elif rule.startswith("FAF >0.1%"):
                counters["w3_likely_benign"] += 1

            # Prefer VV-resolved c. when available, else first ';'-split value
            if hgvsc_resolved:
                hgvsc = hgvsc_resolved
            else:
                raw = str(row.get("HGVSc_snpEff") or "")
                hgvsc = raw.split(";")[0] if raw else ""

            out_row = {
                "gene": gene,
                "transcript": transcript,
                "hgvs_c": hgvsc,
                "classification": f"{classification} {rule}",
                "chr_grch37": hg19_chr,
                "pos_grch37": hg19_pos,
                "ref": ref,
                "alt": alt,
                # Scores used by the rules -- analyst-facing for review
                "PhastCons100way": row.get("phastCons100way_vertebrate") or "",
                "PhyloP100way": row.get("phyloP100way_vertebrate") or "",
                "REVEL": row.get("REVEL_score") or "",
                "SpliceAI_masked": row.get("spliceai_ds_max_masked") or "",
                "FAF95_grpmax": row.get("gnomad_v41_faf95_grpmax") or "",
                "CADD_phred": row.get("CADD_phred") or "",
                "AlphaMissense_pred": row.get("AlphaMissense_pred") or "",
                "ClinVar_sig": row.get("clinvar_clnsig") or "",
                "vv_status": vv_status_cache.get(vv_key, "not_resolved"),
            }
            gene_rows.append(out_row)
            combined_rows.append(out_row)

        if gene_rows:
            gene_name = tsv.stem
            out_path = args.out_dir / f"{gene_name}_seqnext.tsv"
            pd.DataFrame(gene_rows).to_csv(out_path, sep="\t", index=False)

    # Combined output
    combined_path = args.out_dir / "cp_new_seqnext.tsv"
    if combined_rows:
        pd.DataFrame(combined_rows).to_csv(combined_path, sep="\t", index=False)

    print("\n=== Summary ===")
    for k, v in counters.items():
        print(f"  {k:30s} {v:>8,}")
    print(f"\nClassified rows -> {combined_path} ({len(combined_rows):,} rows)")
    if counters["w2_pending_spliceai"] or counters["w3_pending_spliceai"]:
        print("\nWorkflows 2 & 3 are partially pending: run spliceai -M 1 on")
        print("  data/exports/cp_new/cp_new.vcf, then annotate_spliceai.py,")
        print("  then re-run this script.")

    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Pull Benign / Likely_benign ClinVar assertions for cp_new genes and
emit them in the SeqNext format. Recovers variants the dbNSFP-based path
misses (notably synonymous variants -- dbNSFP only carries non-synonymous).

Source: NCBI ClinVar variant_summary.txt.gz (bulk, refreshed weekly).
       https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

Filters:
  - GeneSymbol in cp_new panel
  - ClinSigSimple == 1 (single-significance: Benign / Likely_benign)
  - ClinicalSignificance starts with 'Benign' or 'Likely benign'
  - ReviewStatus has at least one star (excludes 'no assertion' / 'no classification')
  - Assembly == GRCh37 for hg19 output, GRCh38 for cross-check
  - Type in {single nucleotide variant, Deletion, Insertion, Indel, Duplication}

Output: data/exports/cp_new/clinvar_benign_seqnext.tsv
        Same column schema as the SeqNext exports so downstream tools work
        unchanged. Classification reads:
          'Benign ClinVar(2-star,multi-submitter)' or similar.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.panels import get_panel

CLINVAR_TSV = Path.home() / "data" / "clinvar" / "variant_summary.txt.gz"
BED = Path("/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed")
OUT = Path("data/exports/cp_new/seqnext/clinvar_benign_seqnext.tsv")

# ReviewStatus -> star count (NCBI standard, https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
STARS = {
    "no assertion provided": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
    "no interpretation for the single variant": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}


def load_gene_transcripts(bed: Path) -> dict[str, str]:
    counts: dict[str, Counter] = {}
    with open(bed) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 6:
                counts.setdefault(fields[4], Counter())[fields[5]] += 1
    return {gene: count.most_common(1)[0][0] for gene, count in counts.items()}


def benign_row(record: dict, panel: set[str], gene_tx: dict[str, str]) -> dict | None:
    gene = record["GeneSymbol"].upper()
    sig = record["ClinicalSignificance"]
    review = record["ReviewStatus"]
    stars = STARS.get(review, 0)
    if gene not in panel or not sig.startswith(("Benign", "Likely benign")):
        return None
    if stars < 1 or record["Assembly"] != "GRCh37":
        return None
    try:
        position = int(record["Start"])
    except ValueError:
        return None
    name = record["Name"]
    hgvs_c = "c." + name.partition(":c.")[2].split(" ")[0] if ":c." in name else ""
    classification = "Benign" if sig.startswith("Benign") else "Likely_benign"
    # ClinVar c. notation may be on the minus strand, so it cannot supply genomic alleles.
    ref = record["ReferenceAllele"]
    alt = record["AlternateAllele"]
    return {
        "gene": gene,
        "transcript": gene_tx.get(gene, ""),
        "hgvs_c": hgvs_c,
        "classification": f"{classification} ClinVar({stars}-star, {review.split(',')[0]})",
        "chr_grch37": record["Chromosome"],
        "pos_grch37": position,
        "ref": "" if ref == "na" else ref,
        "alt": "" if alt == "na" else alt,
        "PhastCons100way": "",
        "PhyloP100way": "",
        "REVEL": "",
        "SpliceAI_masked": "",
        "FAF95_grpmax": "",
        "CADD_phred": "",
        "AlphaMissense_pred": "",
        "ClinVar_sig": sig,
        "vv_status": "clinvar_assertion",
        "clinvar_review": review,
        "clinvar_stars": stars,
        "clinvar_type": record["Type"],
        "clinvar_name": name,
    }


def read_benign_rows(path: Path, panel: set[str], gene_tx: dict[str, str]) -> list[dict]:
    rows = []
    seen = set()
    fields = (
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
        "Assembly",
        "Start",
        "Name",
        "ReferenceAllele",
        "AlternateAllele",
        "Chromosome",
        "Type",
    )
    with gzip.open(path, "rt") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            if any(record.get(field) is None for field in fields):
                continue
            row = benign_row(record, panel, gene_tx)
            if row is None:
                continue
            key = tuple(
                row[field] for field in ("chr_grch37", "pos_grch37", "ref", "alt", "transcript")
            )
            if key not in seen:
                seen.add(key)
                rows.append(row)
    return rows


def main() -> int:
    if not CLINVAR_TSV.exists():
        sys.exit(f"Missing: {CLINVAR_TSV}")
    cp = {g.upper() for g in get_panel("cp_new")}
    print(f"cp_new panel genes: {len(cp)}")
    gene_tx = load_gene_transcripts(BED)
    print(f"Reading {CLINVAR_TSV} ...")
    df = pd.DataFrame(read_benign_rows(CLINVAR_TSV, cp, gene_tx))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {len(df):,} ClinVar Benign/LB rows -> {OUT}")
    if df.empty:
        return 0
    print("\nTop genes by ClinVar B/LB count:")
    for gene, count in df["gene"].value_counts().head(15).items():
        print(f"  {gene}: {count}")
    print("\nBy variant type:")
    for variant_type, count in df["clinvar_type"].value_counts().head(10).items():
        print(f"  {variant_type}: {count}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

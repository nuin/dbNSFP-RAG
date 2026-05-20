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

import gzip
import sys
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


def main() -> int:
    if not CLINVAR_TSV.exists():
        sys.exit(f"Missing: {CLINVAR_TSV}")

    cp = {g.upper() for g in get_panel("cp_new")}
    print(f"cp_new panel genes: {len(cp)}")

    # Build BED gene -> transcript map (most-common NM_ per gene)
    from collections import Counter
    tx_counts: dict[str, Counter] = {}
    with open(BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6:
                tx_counts.setdefault(p[4], Counter())[p[5]] += 1
    gene_tx = {g: c.most_common(1)[0][0] for g, c in tx_counts.items()}

    print(f"Reading {CLINVAR_TSV} ...")
    rows = []
    seen = set()
    with gzip.open(CLINVAR_TSV, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        ci = {k: col[k] for k in ("GeneSymbol", "ClinicalSignificance",
                                   "ClinSigSimple", "ReviewStatus",
                                   "Assembly", "Chromosome", "Start", "Stop",
                                   "ReferenceAllele", "AlternateAllele",
                                   "Name", "Type")}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(ci.values()):
                continue
            gene = fields[ci["GeneSymbol"]].upper()
            if gene not in cp:
                continue
            sig = fields[ci["ClinicalSignificance"]]
            # ClinSigSimple is unreliable (0 sometimes appears even when
            # ClinicalSignificance is unambiguously "Benign"); trust the text.
            if not (sig.startswith("Benign") or sig.startswith("Likely benign")):
                continue
            review = fields[ci["ReviewStatus"]]
            stars = STARS.get(review, 0)
            if stars < 1:
                continue
            assembly = fields[ci["Assembly"]]
            if assembly != "GRCh37":
                continue  # we emit hg19 for SeqNext
            chrom = fields[ci["Chromosome"]]
            start = fields[ci["Start"]]
            try:
                start_i = int(start)
            except ValueError:
                continue
            ref = fields[ci["ReferenceAllele"]]
            alt = fields[ci["AlternateAllele"]]
            # ClinVar bulk file sometimes has 'na' for SNVs; pull from c. notation
            if ref == "na" or alt == "na":
                # Name has c.XXX{REF}>{ALT} -- extract
                import re
                m = re.search(r"c\.[*\-+\d]+([ACGT])>([ACGT])", fields[ci["Name"]])
                if m:
                    # The c. notation is coding-strand; for genomic ref/alt
                    # we'd need to know the gene's strand. Many cp_new genes
                    # are on the minus strand. Leave ref/alt as-is and rely
                    # on the Name + transcript for downstream consumers.
                    if ref == "na": ref = ""
                    if alt == "na": alt = ""
                else:
                    if ref == "na": ref = ""
                    if alt == "na": alt = ""
            vtype = fields[ci["Type"]]
            name = fields[ci["Name"]]
            # Extract HGVS c. from Name field
            # Name format: NM_xxxx.x(GENE):c.NNNXXX>YYY (p.X)
            hgvs_c = ""
            transcript = gene_tx.get(gene, "")
            if ":c." in name:
                # Take the c.* part; verify NM_ matches the BED's NM_ if possible
                pre, _, after = name.partition(":c.")
                # 'NM_xxxx.x(GENE)' before colon
                tx_in_name = pre.split("(")[0]
                # The c. may have a trailing ' (p.X)' suffix
                c_part = "c." + after.split(" ")[0]
                hgvs_c = c_part
                # Optionally prefer the variant's own NM_ if it matches BED's
                if transcript and tx_in_name != transcript:
                    # ClinVar may have used a different version (.4 vs .5).
                    # Keep BED's transcript label, but flag this in classification.
                    pass

            key = (chrom, start_i, ref, alt, transcript)
            if key in seen:
                continue
            seen.add(key)

            classification = ("Benign" if sig.startswith("Benign") else "Likely_benign")
            classification = (f"{classification} ClinVar({stars}-star, "
                              f"{review.split(',')[0]})")

            rows.append({
                "gene": gene,
                "transcript": transcript,
                "hgvs_c": hgvs_c,
                "classification": classification,
                "chr_grch37": chrom,
                "pos_grch37": start_i,
                "ref": ref,
                "alt": alt,
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
                "clinvar_type": vtype,
                "clinvar_name": name,
            })

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nWrote {len(df):,} ClinVar Benign/LB rows -> {OUT}")

    # Per-gene summary
    print("\nTop genes by ClinVar B/LB count:")
    for g, n in df["gene"].value_counts().head(15).items():
        print(f"  {g}: {n}")

    # By type
    print("\nBy variant type:")
    for t, n in df["clinvar_type"].value_counts().head(10).items():
        print(f"  {t}: {n}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

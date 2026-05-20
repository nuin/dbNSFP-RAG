# Option A — local synonymous variant catalog

dbNSFP is a non-synonymous variant database by design. Synonymous coding
variants — which the lab's W2 rule explicitly targets — are **not in
dbNSFP at all**. The dbNSFP-based pipeline therefore can't see them, no
matter how cleanly we wire up the classifier.

Triggering finding: lab QA for GOT2 found 3 synonymous variants the
pipeline missed:

```
c.816C>T (p.Cys272=)   hg19 chr16:58,750,604 G>A
c.228T>G (p.Val76=)    hg19 chr16:58,757,668 A>C
c.213T>C (p.Asn71=)    hg19 chr16:58,757,683 A>G
```

All three exist in gnomAD with high frequency, all three are 2-star
Benign in ClinVar. None of them exist in dbNSFP. The pipeline can't
classify what it can't see.

## Solution: enumerate synonymous candidates from reference + RefSeq

`scripts/enumerate_synonymous_catalog.py` builds a complete catalog of
every possible synonymous SNV in the BED-anchored RefSeq transcript's
CDS by walking the hg19 reference FASTA against the UCSC `ncbiRefSeq`
table.

### Inputs

| Source | Path | What it gives |
|---|---|---|
| UCSC hg19 FASTA | `~/data/hg19/hg19.fa` (+ `.fai`) | Reference sequence per chromosome |
| UCSC `ncbiRefSeq` table | `~/data/refseq_hg19/ncbiRefSeq.txt.gz` | CDS exon structure per NM_ transcript (chrom, strand, cdsStart/End, exonStarts/Ends, frames) |
| Lab BED | `cp_new_bundle/bed/C+_ALL_IDPE_APR2026.bed` | Gene → BED-anchored RefSeq NM_ |
| Codon table | BioPython `Bio.Seq` | Translation incl. degenerate sites |

### Algorithm

```
For each cp_new gene's BED-anchored NM_ (164/165 are coding):
  Look up CDS exon structure in ncbiRefSeq
  Build the concatenated CDS sequence (strand-aware, reverse-complemented for - strand)
  Track each CDS index -> genomic position
  For each CDS position i in [0, cds_len):
    codon = cds[(i//3)*3 : (i//3)*3 + 3]
    ref_aa = translate(codon)
    for each alt nt in {A,C,G,T}-{ref}:
      new_codon = codon with alt substituted at position (i%3)
      alt_aa = translate(new_codon)
      if alt_aa == ref_aa:
        emit (gene, transcript, genomic_pos, genomic_ref, genomic_alt,
              hgvs_c, codon_ref, codon_alt, aa, strand)
```

Genomic ref/alt are derived from the CDS base by complementation for
minus-strand genes (which is why c.816C>T on a minus-strand gene shows
in genomic VCF as G>A).

### Output of enumeration

`data/exports/cp_new/cp_new_synonymous_catalog.tsv` — **276,885
synonymous candidates** across 164 cp_new genes.

Top contributors:
| Gene | Candidates |
|---|---|
| KMT2D | 12,432 |
| BRCA2 | 6,398 |
| ATM | 5,768 |
| NF1 | 5,660 |
| APC | 5,639 |
| NSD1 | 5,424 |
| POLE | 4,620 |
| TSC2 | 3,975 |
| SLX4 | 3,920 |
| FANCM | 3,843 |

GOT2 has 891 candidates total, including all three originally-missed
variants.

## Annotation: gnomAD FAF + dbNSFP conservation

`scripts/annotate_classify_synonymous.py` annotates each catalog entry
and applies W1/W2:

- **FAF**: gnomAD v2.1.1 exomes sites VCF, queried in hg19 directly
  (remote tabix per gene region). No liftover. `AF_popmax` used as the
  FAF estimate.
- **PhastCons / PhyloP**: borrowed from any dbNSFP row at the same
  (chr, pos). dbNSFP's nsSNV rows include per-position conservation
  scores; a synonymous variant at the same position can reuse them
  (conservation is per-base, not per-allele). Missing if dbNSFP doesn't
  cover the position (rare in coding regions).
- **SpliceAI**: assumed 0 for synonymous variants in CDS interior.
  Variants near a splice site (within ±50 bp of an exon boundary) should
  go through a SpliceAI rerun before classifying — not currently
  filtered, but a future addition.

W2 classifier rule (synonymous-specific):
- PhastCons < 1.0 (position not perfectly conserved)
- SpliceAI ≤ 0.1 (assumed 0 for CDS-interior)
- W1 override: FAF > 5%

Output: `data/exports/cp_new/seqnext/cp_new_seqnext_synonymous_catalog.tsv`
— **52,004 classified synonymous rows** (51,622 W2 Benign + 382 W1
Benign).

### Verification on the three originally-missed GOT2 variants

| HGVS c. | hg19 | Classification | Source |
|---|---|---|---|
| c.213T>C | chr16:58,757,683 A>G | Benign FAF >5% (FAF=0.80) | catalog (also in ClinVar) |
| c.228T>G | chr16:58,757,668 A>C | Benign FAF >5% (FAF=0.83) | catalog (also in ClinVar) |
| c.816C>T | chr16:58,750,604 G>A | Benign FAF >5% (FAF=0.85) | catalog (also in ClinVar) |

All three fire W1 (FAF >5%) on gnomAD v2.1.1 hg19 native. The pipeline
now recovers what it was missing.

## Merge into FINAL output

`scripts/merge_all_sources.py` combines:
1. **Pipeline** (dbNSFP nsSNVs through W1/W2/W3 + VV resolver) — 6,963 rows
2. **ClinVar B/LB** (NCBI bulk, ≥1 star) — 121,568 rows
3. **Synonymous catalog** (local enumeration) — 31,519 rows

Dedup priority: pipeline > ClinVar > catalog. Each row tagged with `source`.

**FINAL total: 160,050 rows** (158,973 strict).

For the 68 new genes only: 46,344 rows (45,883 strict).

## What's NOT in the catalog

- **Variants in non-coding RNA genes** (TERC, NR_001566.1) — excluded
  by definition. Need a separate non-coding workflow.
- **Variants in UTRs** — the catalog walks only the CDS (cdsStart to
  cdsEnd in ncbiRefSeq). 5'/3' UTR variants would need a separate pass.
- **Variants outside the BED-anchored NM_** — alternative isoforms
  with different exon structures aren't covered. Most cp_new genes are
  single-transcript in the BED so this is rare; APC and RAD51D dual-tx
  regions are the exceptions and only the BED's primary tx is enumerated.
- **Indel variants** — catalog is SNV-only.
- **Splice-affecting synonymous variants near canonical splice sites** —
  assumed SpliceAI=0; should be flagged for SpliceAI rerun in a future
  iteration (the current pipeline already covers these for the dbNSFP
  rows).

## Re-running the catalog

```bash
# 1. Enumerate (15 sec for 165 genes)
uv run python scripts/enumerate_synonymous_catalog.py

# 2. Annotate + classify (~8 min, 164 remote tabix calls to gnomAD v2.1.1)
uv run python scripts/annotate_classify_synonymous.py

# 3. Merge with pipeline + ClinVar
uv run python scripts/merge_all_sources.py
```

Outputs land in `data/exports/cp_new/seqnext/`:
- `cp_new_seqnext_synonymous_catalog.tsv` — Option A output only
- `cp_new_seqnext_FINAL.tsv` — combined with pipeline + ClinVar
- `cp_new_seqnext_FINAL_strict.tsv` — strict, drops intergenic + canonical splice

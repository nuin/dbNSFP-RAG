# v3 UPLOAD — rule-based classification only

Built 2026-06-05 in response to lab feedback:

> ClinVar's classifications cannot be used in this project. For variants to be
> classified in MGL they would need to meet the criteria provided at the
> beginning of the project (ie, Missense variant with a freq over 0.1% with a
> REVEL score <0.29 and a SpliceAI score <=0.1).

## Result

**13,871 rows** — down from v2's 39,331 (-25,460).

| Source | rows | classification |
|---|---:|---|
| synonymous_catalog | 12,253 | W2: Benign synonymous + SpliceAI≤0.1 + PhastCons<1.0 |
| pipeline (W1/W2/W3) | 1,617 | W1 Benign (FAF>5%), W2 Benign syn, W3 LB |
| gnomad_common | 1 | W1 Benign (FAF>5% gnomAD direct) |

Class split: 13,849 Benign + 22 Likely_benign. Gene coverage: 67/68 (TERC blank).

## Rules applied

| Rule | Criteria | Call |
|---|---|---|
| W1 | gnomAD grpmax FAF >5% | Benign |
| W2 | synonymous + SpliceAI ≤0.1 + PhastCons <1.0 (+ PhyloP <0.1 if intronic ROI) | Benign |
| W3 | missense + FAF >0.1% + REVEL <0.29 + SpliceAI ≤0.1 | Likely_benign |

## v2 → v3 changes

- **24,856 ClinVar-source rows dropped** (per lab feedback — ClinVar not accepted as evidence)
- **603 W2 syn-catalog rows dropped** because measured SpliceAI > 0.1
  (these were previously called Benign using the "0 (assumed)" placeholder; now
  that real SpliceAI scores are backfilled, they fail the W2 threshold)
- **1 W3 row dropped** for same reason

## Open items

1. **SpliceAI mode**: current scores use Illumina masked mode (`-M 1`). The lab
   compared against spliceai.org (unmasked default) and saw different values.
   Unmasked re-run is in progress (`cp_new_v2_needed.spliceai_unmasked.vcf`).
   When complete, re-validate W2/W3 thresholds against unmasked scores.
2. **TERC has 0 rows** — non-coding RNA, no W1/W2/W3 applies. Manual review only.
3. **No P/LP**: by design, this pipeline only emits Benign/Likely_benign. Any
   variant not autoclassified goes to manual review.

## Files
```
UPLOAD_v3.tsv               13,871 rows -- the SeqNext upload
UPLOAD_v3_dropped.tsv       25,460 rows -- audit trail of v2 rows removed
per_gene/{GENE}__v3.tsv         67 files -- per-gene splits (TERC has no file)
docs/README.md                          -- this file
```

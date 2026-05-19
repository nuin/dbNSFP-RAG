# cp_new — net-new genes only (vs. old CP panel)

Filtered slice of the full cp_new pipeline output, restricted to the **68
genes added in CP_new (APR2026)** that weren't in the previous CP panel
(01JUN2021). Same classifications, same scores, same review piles — just
the new-panel subset.

## Source diff

```
CP (01JUN2021):  104 genes
CP_new (APR2026): 165 genes
NEW in CP_new:    68 genes
DROPPED from CP:   7 genes
```

Gene list at `_new_genes.txt`.

## Counts (this subset)

| File | Rows |
|---|---|
| `cp_new_seqnext_strict.tsv` (recommended upload) | **1,728** |
| `cp_new_seqnext_minimal.tsv` (everything) | 2,192 |
| `cp_new_seqnext_combined.tsv` (audit with chr/pos/ref/alt/vv_status) | 2,192 |
| `_review_canonical_splice.tsv` | 93 |
| `_review_intergenic.tsv` | 442 |
| `per_gene/{GENE}_seqnext.tsv` files | 66 |

2 of the 68 new genes had zero Benign/LB classifications (silent — heavy
purifying selection / conservation prevents the rules from firing).

## Schema

Every row carries the four core SeqNext fields plus the scores that
fired the rule:

```
gene  transcript  hgvs_c  classification
PhastCons100way  PhyloP100way  REVEL  SpliceAI_masked  FAF95_grpmax
CADD_phred  AlphaMissense_pred  ClinVar_sig
```

Intronic calls also tag the offset in the classification string
(`..., intronic_ROI(-5)`) so the analyst can see at a glance the variant
sits in the -15 to +6 splice region with PhyloP also confirmed.

## Three classification rules (same as full panel)

| Rule | Logic |
|---|---|
| **W1 Benign** | gnomAD FAF >5% |
| **W2 Benign** | synonymous AND SpliceAI≤0.1 AND PhastCons<1.0; intronic also requires PhyloP<0.1 AND offset in [-15, +6] |
| **W3 Likely_benign** | FAF >0.1% AND REVEL<0.29 AND SpliceAI≤0.1 |

## Reproduce this subset

```bash
uv run python scripts/filter_new_genes.py
```

Reads `/Users/nuin/Projects/ahs/new_bed/CP/C+_ALL_IDPE_01JUN2021.bed` and
`/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed`, diffs
the gene lists, filters every SeqNext output to the new-only set.

For the full 165-gene pipeline output, see the parent `outputs/` directory
and `cp_new_bundle/SUMMARY.md`.

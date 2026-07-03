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

## Pipeline diagram

```
 ┌──────────────────────────┐     ┌──────────────────────────────┐
 │  dbNSFP 5.3.1a           │     │  gnomAD v4.1 joint sites VCF │
 │  per-chromosome chr*.gz  │     │  (remote tabix over HTTPS,   │
 │  ~1B nsSNV rows total    │     │   per-gene region pull)      │
 └────────────┬─────────────┘     └────────────────┬─────────────┘
              │                                    │
              └────────────┐    ┌──────────────────┘
                           ▼    ▼
              ┌─────────────────────────────┐
              │  scripts/cp_new.py          │
              │  - filter to 165 panel genes│
              │  - join gnomAD FAF + AF     │
              │  - emit per-gene TSVs       │
              │  - emit consolidated VCF    │
              └────────────┬────────────────┘
                           ▼
              165 per-gene TSVs (1,057,357 rows)
              data/exports/cp_new/cp_new.vcf (1,045,306 unique variants)
                           │
                           ▼
              ┌─────────────────────────────┐
              │  SpliceAI v1.3 (-M 1)       │
              │  local install, GENCODE V24 │
              │  4-way parallel on M1 Ultra │
              │  ~4 hr wall time            │
              └────────────┬────────────────┘
                           ▼
              cp_new.spliceai.vcf
                           │
                           ▼
              ┌─────────────────────────────┐
              │  annotate_spliceai.py       │
              │  join masked DS_MAX into    │
              │  the per-gene TSVs in-place │
              └────────────┬────────────────┘
                           ▼
              ┌─────────────────────────────┐                ┌──────────────────────────┐
              │  classify_cp_new.py         │ ◄──────────────┤  BED (APR2026)           │
              │  apply W1/W2/W3 rules:      │  transcript    │  C+_ALL_IDPE_APR2026.bed │
              │  - missing SpliceAI = 0     │  lookup        │  165 genes / 2,379 lines │
              │  - intronic ROI: -15 to +6  │                │  RefSeq NM_ per region   │
              │  - intronic needs PhyloP    │                └──────────────────────────┘
              │  - reads VV cache for       │
              │    BED-anchored c.          │
              └────────────┬────────────────┘
                           ▼
              cp_new_seqnext.tsv (6,971 classified)
                           │
                           ▼
              ┌─────────────────────────────┐
              │  resolve_hgvs_vv.py         │
              │  VariantValidator REST per  │
              │  row -> authoritative c. on │
              │  BED's NM_                  │
              │  (1.2s/call, ~5 hr cold)    │
              └────────────┬────────────────┘
                           ▼
              cp_new_seqnext.tsv with vv_status + corrected hgvs_c
              + cache at .vv_cache.json (reused by future classify runs)
                           │
                           ▼
              ┌─────────────────────────────┐
              │  build_cp_new_exports.py    │
              │  --strict drops:            │
              │   - vv_status != ok         │
              │   - canonical splice -1/-2  │
              │                  /+1/+2     │
              └────────────┬────────────────┘
                           ▼
              ┌─ cp_new_seqnext_strict.tsv   (5,899 cp_new-wide)
              ├─ cp_new_seqnext_minimal.tsv  (6,971)
              ├─ cp_new_seqnext_combined.tsv (6,971 + audit cols)
              ├─ cp_new_all_annotations.tsv.gz (1,057,357 x 67 cols)
              ├─ cp_new_viz.db (visualization SQLite)
              └─ per_gene/{GENE}_seqnext.tsv (163 files)
                           │
                           ▼
              ┌─────────────────────────────┐                ┌──────────────────────────┐
              │  scripts/filter_new_genes.py│ ◄──────────────┤  CP (01JUN2021) BED      │
              │  diff CP_new ∖ CP gene sets │  gene diff     │  104 genes               │
              │  slice every output to the  │                └──────────────────────────┘
              │  68 net-new genes only      │
              └────────────┬────────────────┘
                           ▼
                  new_genes_only/  (this folder)
                  ┌─ cp_new_seqnext_strict.tsv   1,728 rows   ◄── upload to SeqNext
                  ├─ cp_new_seqnext_minimal.tsv  2,192 rows
                  ├─ cp_new_seqnext_combined.tsv 2,192 rows
                  ├─ _review_canonical_splice.tsv   93 rows   ── manual review
                  ├─ _review_intergenic.tsv        442 rows   ── drop or re-resolve
                  └─ per_gene/*.tsv  66 files (2 silent genes)
```

## Counts (this subset of 68 new genes)

### By rule

| Rule | Logic | Hits |
|---|---|---|
| **W1 Benign** | gnomAD FAF >5% | **126** |
| **W2 Benign** | synonymous + SpliceAI≤0.1 + PhastCons<1.0 (+ PhyloP & ROI for intronic) | **2,040** |
| **W3 Likely_benign** | FAF >0.1% + REVEL<0.29 + SpliceAI≤0.1 | **26** |
| **Total** | | **2,192** |

### By output file

| File | Rows |
|---|---|
| `cp_new_seqnext_strict.tsv` (recommended upload) | **1,728** |
| `cp_new_seqnext_minimal.tsv` (everything) | 2,192 |
| `cp_new_seqnext_combined.tsv` (audit with chr/pos/ref/alt/vv_status) | 2,192 |
| `_review_canonical_splice.tsv` | 93 |
| `_review_intergenic.tsv` | 442 |
| `per_gene/{GENE}_seqnext.tsv` files | 66 |

### Top contributors

| Gene | Calls |
|---|---|
| TRPV6 | 451 |
| GCM2 | 397 |
| HOXB13 | 180 |
| GNA11 | 171 |
| CBL | 86 |
| AP2S1 | 79 |
| SMARCA4 | 71 |
| SUCLG2 | 68 |
| MITF | 65 |
| KIF1B | 55 |

### Silent genes (0 calls)

5 of the 68 new genes had zero Benign/LB classifications — heavy purifying
selection / conservation prevents the rules from firing, or the gene is
non-coding (TERC is a non-coding RNA so it's not in dbNSFP at all):

`CDKN2C, ETV6, PTH, SPRED1, TERC`

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

# cp_new_bundle/outputs/ — layout

Five top-level dirs. **The current SeqNext upload is `new_genes_68_v2/UPLOAD_v2_scored.tsv`.**

```
outputs/
├── _raw/                   2.7 GB   source intermediates (regenerable, kept for audit)
├── _shared/                376 KB   cross-version inputs (gene whitelist, ranges, review piles)
├── all_panel_165/           89 MB   full 165-gene cp_new panel scope
├── new_genes_68_v1/         39 MB   May/early-Jun delivery (68 net-new genes)
└── new_genes_68_v2/         20 MB   current delivery (68 net-new genes, gnomAD-inclusive)
```

## new_genes_68_v2/ — current delivery
```
UPLOAD_v2_scored.tsv             39,331 rows -- the SeqNext upload (column D = classification)
UPLOAD_v2.tsv                    39,331 rows -- same, pre-backfill scores
gnomad_common_standalone.tsv        172 rows -- gnomAD FAF>5% pull with HGVSc resolved
gnomad_common_needs_hgvs.tsv     17,243 rows -- noncoding gnomAD-commons pending VV
per_gene/{GENE}__v2.tsv             68 files -- per-gene splits
docs/
  README.md                       overview + v1 vs v2 comparison
  EMAIL_UPDATE.md                 lab-facing summary
  V2_ISSUES.md                    QA findings (7 items)
```

## new_genes_68_v1/ — historical reference
```
UPLOAD.tsv                       5,796 rows -- v1 production
UPLOAD_conservative.tsv          1,335 rows -- v1 undercall (A+B+C+D filter)
UPLOAD_conservative_dropped.tsv  4,461 rows -- audit trail
FINAL.tsv                       46,344 rows -- pre-filter merge of all sources
FINAL_strict.tsv                            -- strict variant
gnomad_common.tsv                            -- initial gnomAD pull (pre-reannotation)
gnomad_common.raw.tsv                        -- backup of original pull
intermediate/                    8 files     -- combined/minimal/strict steps + review piles
per_gene/                      272 files     -- 4 tag variants x 68 genes
                                                ({GENE}__full|pipeline|with_clinvar|conservative.tsv)
docs/
  SUMMARY.md, UPLOAD_ISSUES.md, EMAIL_UPDATE.md, EMAIL_REPLY.md
```

## all_panel_165/ — full panel scope (165 genes)
```
FINAL.tsv                       160,050 rows
FINAL_strict.tsv
intermediate/                   6 progression TSVs
per_gene/                       163 files ({GENE}_seqnext.tsv)
```

## _shared/
```
new_genes.txt                   68-gene whitelist for new_genes_68_*
gene_hg38_ranges.json           hg38 coords per gene (used by gnomAD pulls)
review_canonical_splice.tsv     manual-review pile
review_intergenic.tsv           manual-review pile
```

## _raw/ — source dumps (regenerable, kept for audit)
```
clinvar_benign.tsv              29 MB    ClinVar bulk B/LB pull
synonymous_catalog.tsv          21 MB    Option A local enumeration
all_annotations.tsv.gz         155 MB    dbNSFP + gnomAD + SpliceAI per variant
viz.db                         2.7 GB    SQLite for viz tool
gnomad_v2_summary.txt                    v2 vs v4 FAF comparison report
```

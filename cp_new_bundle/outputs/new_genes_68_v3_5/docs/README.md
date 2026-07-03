# v3.5 — consolidated release (all genes)

v3.5 = v3.4 (c.-based gnomAD anchor, liftover-free) with every fix applied,
promoted to a single clean release. This is the definitive deliverable.

## What's in it (all fixes through 2026-07)

- gnomAD-anchored on `(gene, NM_, c.)` — liftover-free (position fallback only
  for DPYD, GPR161, MITF, SMARCA4 where project NM_ ≠ MANE)
- ClinVar variants kept but reclassified purely by W1/W2/W3 (no ClinVar evidence)
- Intronic ROI [-20, +10]
- 5'UTR (c.-N) variants removed
- PhastCons + PhyloP filled for 100% of rows (UCSC hg19 100-way tracks)
- Synonymous with PhastCons ≥ 1.0 kept but labeled "conserved (review)"
- Missense SpliceAI measured for all (masked + unmasked); those >0.1 dropped
- g_hg19 / gnomad_id / g_hg38 / p_hgvs columns for gnomAD v4.1 lookup

## Numbers

Total: **30,201 rows**, 67/68 genes (TERC = non-coding RNA, empty per-gene file).

| classification | rows |
|---|---:|
| Benign W2 (synonymous, PhastCons<1.0) | 21,430 |
| Synonymous conserved (review, not auto-benign) | 7,886 |
| Benign W1 (FAF>5%) | 469 |
| Likely_benign W3 (missense) | 416 |

## Columns

gene, transcript, hgvs_c, classification, PhastCons100way, PhyloP100way, REVEL,
SpliceAI_masked, SpliceAI_unmasked, FAF95_grpmax, CADD_phred, AlphaMissense_pred,
ClinVar_sig, chr_grch37, pos_grch37, ref, alt, g_hg19, gnomad_id, g_hg38,
p_hgvs, vv_status, source

- **classification** (col D) is what SeqNext imports
- **gnomad_id** (e.g. 5-80768028-G-A) pastes into gnomAD v4.1 search
- **g_hg19** native; **g_hg38** via liftover (hg19 shop)

## Files
```
UPLOAD_v3_5.tsv            30,201 rows -- the SeqNext upload
per_gene/{GENE}__v3_5.tsv  68 files (TERC empty)
docs/README.md             this file
```

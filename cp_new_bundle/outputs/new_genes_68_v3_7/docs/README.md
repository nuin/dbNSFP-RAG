# v3.7 — consolidated release (real SpliceAI everywhere)

v3.7 = v3.4 (c.-based gnomAD anchor, liftover-free) with the SpliceAI honesty
fixes from the 2026-07-15 lab review. This is the current definitive deliverable.

## What changed from v3.5

The lab flagged that many rows showed SpliceAI "0 (assumed)" — a placeholder,
not a measured score (HOXB13 c.600A>G was the example). Fixed:

1. **Ran real SpliceAI (masked + unmasked) on all 20,425 "assumed" positions.**
   No row is scored by assumption anymore.
2. **Fixed stale labels** — 5,722 rows had a *measured* score but the label text
   still said "(assumed 0)". Now reads "(measured)".
3. **Excluded genuinely-unscored rows** (SpliceAI couldn't score them — ref
   mismatch / contig edge): 149 synonymous + ~7 missense. Per the lab: "if these
   scores are actually just assumed and not based on SpliceAI data, we would want
   these rows excluded." Saved to `*_assumed_excluded.tsv`.
4. **1,555 synonymous with real SpliceAI > 0.1** (were wrongly "assumed 0 →
   benign") relabeled "splice signal - review, not benign" with the real score.
5. **All 414 missense now carry both masked + unmasked SpliceAI.**

## Numbers

Total: **30,050 rows**, 67/68 genes (TERC = non-coding RNA, empty per-gene file).

| classification | rows |
|---|---:|
| Benign W2 (synonymous, SpliceAI≤0.1 measured, PhastCons<1.0) | 20,554 |
| Synonymous conserved (review, PhastCons≥1.0) | 7,061 |
| Synonymous splice-signal (review, SpliceAI>0.1) | 1,555 |
| Benign W1 (FAF>5%) | 466 |
| Likely_benign W3 (missense) | 414 |

Integrity: every W2/W3 row has real measured SpliceAI (masked + unmasked);
0 "assumed" anywhere; 100% PhastCons + gnomad_id; no 5'UTR; ROI [-20,+10].
(Benign W1 common variants have no SpliceAI by design — splice is irrelevant to
a >5% population frequency call.)

## Columns

gene, transcript, hgvs_c, classification, PhastCons100way, PhyloP100way, REVEL,
SpliceAI_masked, SpliceAI_unmasked, FAF95_grpmax, CADD_phred, AlphaMissense_pred,
ClinVar_sig, chr_grch37, pos_grch37, ref, alt, g_hg19, gnomad_id, g_hg38,
p_hgvs, vv_status, source

- **classification** (col D) = the SeqNext import value
- **SpliceAI_masked** (`-M 1`) / **SpliceAI_unmasked** (`-M 0`, = spliceai.org default); unmasked ≥ masked always — see `docs/SPLICEAI_MASKED_VS_UNMASKED.md`
- **gnomad_id** (e.g. 5-80768028-G-A) pastes into gnomAD v4.1
- **g_hg19** native; **g_hg38** via liftover (hg19 shop)

## Open items (do not block this release)

- phyloP threshold source for the conservation gate (needs lab confirmation of
  which track their 7.2 is calibrated to — phyloP is source-dependent, unlike
  phastCons)
- Benign vs Likely_benign tier for BP4+BP7 synonymous
- keep-and-flag vs drop for the conserved / splice-signal review piles

## Files
```
UPLOAD_v3_7.tsv            30,050 rows -- the SeqNext upload
per_gene/{GENE}__v3_7.tsv  68 files (TERC empty)
docs/README.md             this file
```

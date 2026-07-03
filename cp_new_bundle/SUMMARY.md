# cp_new — one-page overview

**Goal**: Classify variants in 165 hereditary-cancer genes as Benign / Likely_benign using three lab-specific rules, output SeqNext-ready files anchored to the lab's BED transcripts.

## Pipeline (7 stages)

```
dbNSFP per-chr → cp_new.py → per-gene TSVs (+ consolidated VCF)
                                       ↓
gnomAD v4.1 (remote tabix) — joined in same step (FAF95_grpmax)
                                       ↓
SpliceAI -M 1, 4-way parallel local run, ~4 hr
                                       ↓
annotate_spliceai.py → joins masked DS_MAX into per-gene TSVs
                                       ↓
classify_cp_new.py → applies W1/W2/W3, ROI filter, VV-aware intronic detection
                                       ↓
resolve_hgvs_vv.py → VariantValidator HGVS resolution on BED's RefSeq NM_
                                       ↓
build_cp_new_exports.py → minimal + strict + all-annotations + viz SQLite
```

## Three rules

| Rule | Logic | Hits |
|---|---|---|
| **W1 Benign** | gnomAD FAF >5% | 328 |
| **W2 Benign** | synonymous AND SpliceAI≤0.1 AND PhastCons<1.0; intronic also requires PhyloP<0.1 + offset in [-15,+6] | 6,524 |
| **W3 Likely_benign** | FAF >0.1% AND REVEL<0.29 AND SpliceAI≤0.1 | 119 |

## Final numbers

- Raw dbNSFP rows across the 165 panel genes: **1,057,357**
- Classified Benign / Likely_benign total: **6,971**
- Upload-clean (strict, drops intergenic + canonical-splice review piles): **5,899**
- Manual-review files: 234 canonical-splice + 981 intergenic = **1,215**

## Bugs found and fixed during build

1. **dbNSFP multi-transcript `genename`** (e.g. `POLE;POLE;POLE`) was being exact-matched in extraction → 90% of variants silently dropped. Fixed: match any `;`-split part.
2. **Classifier was dropping variants outside the BED interval** for their own gene — liftover variance lost legitimate calls. Fixed: use the gene's primary BED transcript as fallback for HGVS reporting.
3. **Intronic detection on dbNSFP's HGVSc (often non-MANE)** missed variants that are intronic on the BED's NM_ — PhyloP gate never fired. Fixed: classifier loads the VariantValidator cache and uses the BED-anchored c. for intronic detection.
4. **Deep-intronic variants were being auto-classified as W2 Benign**. Fixed: ROI [-15, +6] required; positions outside are ineligible for W2 and go to ACMG review.
5. **SpliceAI v1.3 silently skips alt-transcript regions** (APC E01 on NM_001127511.3, RAD51D alt-E3 on NM_001142571.2). Fixed: missing SpliceAI is treated as 0 (no detected splice signal).

## Deliverables (in `cp_new_bundle/`, 2.7 GB)

| File | Purpose | Rows | Cols |
|---|---|---|---|
| `outputs/cp_new_seqnext_strict.tsv` | **Recommended for SeqNext upload.** Excludes intergenic + canonical-splice review piles. | 5,899 | 12 |
| `outputs/cp_new_seqnext_minimal.tsv` | Everything classified, no QA filter. | 6,971 | 12 |
| `outputs/cp_new_seqnext_combined.tsv` | Audit version with chr/pos/ref/alt/vv_status columns appended. | 6,971 | 16 |
| `outputs/cp_new_all_annotations.tsv.gz` | Complete per-variant annotation dump. | 1,057,357 | 67 |
| `outputs/cp_new_viz.db` | Standalone visualization SQLite (Datasette / DB Browser / Metabase). | 1,057,357 | — |
| `outputs/per_gene/{GENE}_seqnext.tsv` | Same data sliced by gene (163 files). | — | 12 |
| `outputs/review/_review_canonical_splice.tsv` | Manual review: canonical splice positions. | 234 | — |
| `outputs/review/_review_intergenic.tsv` | Manual review: non-coding on BED's NM_. | 981 | — |

**Every SeqNext row carries the scores that fired the rule**:
`PhastCons100way`, `PhyloP100way`, `REVEL`, `SpliceAI_masked`, `FAF95_grpmax`, `CADD_phred`, `AlphaMissense_pred`, `ClinVar_sig`.

Rule strings show the intronic offset for intronic calls, e.g. `Benign synonymous, SpliceAI<=0.1, PhastCons<1.0, PhyloP<0.1, intronic_ROI(-5)` — analyst can see at a glance that all gates fired.

## Open items

- **234 canonical-splice** calls (-1/-2/+1/+2 positions): SpliceAI masked is near-zero at canonical sites by design, so these pass the gate trivially. Need clinical review before upload.
- **981 intergenic flags**: VariantValidator reports the variant as non-coding on the BED's NM_ (it's on an alternative isoform). Drop from upload or re-resolve.
- **~1.05M unclassified variants** in panel genes don't match any of the three rules. They need traditional ACMG/AMP evaluation via `src/acmg_scoring.py` — out of scope for this pipeline.

## More

- Full reference doc: `docs/cp_new_pipeline.md`
- Process / decision narrative: `docs/cp_new_process_journal.md`
- Bundle layout & per-file index: `README.md`

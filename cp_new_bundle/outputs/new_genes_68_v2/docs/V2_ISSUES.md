# v2 UPLOAD — data quality issues

Input: `cp_new_seqnext_UPLOAD_v2.tsv` (39,331 rows, 68/68 genes).

## Sanity checks (passed)

- No residual duplicate (gene, transcript, hgvs_c) keys after the two-pass dedup
- 0 rows missing transcript
- All 68 genes have ≥1 row (TERC has the minimum at 1 — expected for non-coding RNA)

## Issues to surface (priority-ordered)

### #1 — 12,753 W2 "SpliceAI assumed 0" calls (32% of upload)
These are synonymous_catalog rows where SpliceAI wasn't actually run; the
pipeline assumed DS_MAX = 0 because no prediction was available. At 46 rows (v1)
this was a defensible shortcut. At 12,753 rows it's a sizable assumption — any
of them could have a true SpliceAI hit we didn't measure.

**Fix options**: (a) run SpliceAI in batch on the synonymous catalog (Illumina
masked, ~few hours); (b) leave as-is and accept the assumption; (c) flag these
in SeqNext as `Benign provisional — SpliceAI not measured`.

### #2 — 15,693 ClinVar 1-star calls (40% of upload)
Same single-submitter concern as v1, now at larger scale. 9,163 are 2+star
(safer); 14,475 are pipeline-derived no-star (W2/W3 — reliable rules, fine).

**Fix option**: same as v1 — optionally split into `_2star_plus.tsv` vs
`_1star.tsv` so the lab can choose.

### #3 — 7 rows with ClinVar "Conflicting interpretations"
Up from 4 in v1 (the syn-coding unfilter brought 3 more through). Easy fix:
exclude in the ClinVar bulk pull.

### #4 — 428 rows where predictors lean pathogenic
AlphaMissense P/LP, REVEL >0.5, or CADD phred >25 on rows we now ship as Benign.
Concentrated in 4 genes (HOXB13: 112, GNA11: 97, TRPV6: 97, GCM2: 97 — together
~95% of the 428). Worth a per-gene look; common variants in these genes may all
be legitimate population polymorphisms that predictors flag as "damaging"
because they don't have AF features.

### #5 — 30 rows missing HGVSc
Small, but check before upload — these may break SeqNext import.

### #6 — KMT2D is 10% of the upload (4,054 rows)
Largest single gene. Expected given KMT2D size + ClinVar density, but the lab
may want a stricter ClinVar floor (e.g., 2-star only) just for KMT2D.

### #7 — TERC = 1 row only
TERC is a non-coding RNA, not in dbNSFP. Only one classification row exists.
Worth confirming the lab is OK with this (or pulling TERC coords from another
source if more coverage is needed).

## What's NOT an issue anymore

- ✅ Gene coverage: 68/68 (was 67/68 in v1; CDKN2C, ELOC, H3-3A now populated)
- ✅ GOT2 lab-flagged variants (c.816C>T, c.213T>C, c.228T>G) all present
- ✅ Synonymous coding included as Benign (was filtered out in v1)

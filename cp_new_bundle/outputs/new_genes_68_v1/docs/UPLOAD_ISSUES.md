# cp_new UPLOAD — data quality & possible issues

Input: `cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv`  (5,796 rows, 67 of 68 new genes)

## 1. Composition

| Source | Rows | % | What it is |
|---|---:|---:|---|
| ClinVar bulk (B/LB ≥1-star) | 5,622 | 97.0% | Trust the assertion — no scores attached |
| pipeline rules W1/W2/W3 | 174 | 3.0% | Scored: FAF / synonymous+conservation / FAF+REVEL+SpliceAI |
| synonymous_catalog (Option A) | 0 | 0% | **Eliminated by design** (see Issue #1) |

Class split: Benign 1,837 (31.7%) / Likely_benign 3,959 (68.3%). No P/LP/VUS — by design (benign-only autoclassify).

ClinVar review-star mix:

| Stars | Rows | % | Meaning |
|---|---:|---:|---|
| 1-star | 4,448 | 76.7% | criteria provided, **single submitter** |
| 2-star | 1,172 | 20.2% | multiple submitters, no conflicts |
| 3-star | 2 | 0.0% | expert panel |
| no-star | 174 | 3.0% | pipeline rows (not ClinVar) |

---

## 2. Possible issues (priority-ordered)

### #1 — Synonymous catalog produces zero upload rows
Option A enumerated 276,885 synonymous SNVs across cp_new; 12,761 made it into FINAL for the 68 new genes; **all 12,761 were filtered out** of UPLOAD by the "drop pure synonymous coding" rule in `build_new_genes_upload.py` step 5. The catalog *only* survives where the variant is in intronic ROI, but by definition catalog variants are coding-only.

**Impact**: any synonymous variant outside ClinVar will NOT be in upload. The lab's manual review must catch them. The original GOT2 misses (c.816C>T, c.228T>G, c.213T>C) — verify ClinVar covers them now; if not, they'll be missed again.

**Fix options**: (a) leave as-is (matches your "synonyms are clinically useless" feedback) — accept downstream review burden; (b) ship Option A as a *separate* supplementary file so analysts have it as a lookup table without polluting the autoclassify upload.

---

### #2 — ClinVar 1-star calls dominate (76.7%)
4,448 of 5,796 rows come from a single submitter with no expert review. ClinVar 1-star B/LB can later get reclassified or contradicted.

**Fix option**: tag 1-star rows for periodic re-pull, or split the upload into `_seqnext_UPLOAD_2star_plus.tsv` (1,172 rows, high confidence) and `_seqnext_UPLOAD_1star.tsv` (4,448 rows, accept-with-flag).

---

### #3 — ELOC has zero coverage in upload
36 ELOC rows in FINAL, 0 in UPLOAD:
- 25 dropped: `vv_status: flagged:intergenic` (BED-vs-transcript coord mismatch)
- 10 dropped: `vv_status: local_enumeration` (synonymous catalog — see #1)
- 1 dropped: `c.148+219G>C` — deep intronic, outside [-15,+6]

**Fix**: investigate ELOC `NM_005648.4` transcript boundaries vs the BED; the `flagged:intergenic` calls look like the variants ARE in the gene transcript but the BED interval excludes them. Either re-check coords or whitelist ELOC and re-run.

---

### #4 — Predictor conflicts on W1 (FAF>5%) Benign calls
W1 calls anything with grpmax FAF > 5% Benign — but several rows have strong pathogenic-leaning predictors:

| Gene | HGVSc | REVEL | CADD | AlphaMissense (per-tx) |
|---|---|---:|---:|---|
| POT1 | c.548T>G | 0.748 | 25.5 | …LP;P;LP;LP |
| SAMD9L | c.866T>C | 0.441 | 24.7 | LP;.;P;… |
| SUCLG2 | c.1186A>G | 0.757 | 25.1 | …P… |
| SRP72 | c.19G>T | 0.508 | 28.8 | …P;A |
| DPYD | c.496A>G | 0.523 | 26.3 | …B… |

Total: 5 rows where AlphaMissense calls P/LP, 7 rows with CADD>25, on rows we ship as Benign by W1.

**Interpretation**: high population frequency overrides predictor signal in ACMG (BA1 is stand-alone), so these *can* be benign — common variants often score "deleterious" because predictors don't have AF features. But these specific genes are cancer-predisposition; worth manual review before the lab commits to "Benign autoclassify" on `POT1 c.548T>G` etc.

**Fix option**: add a `predictor_conflict` flag column when AlphaMissense=P/LP or REVEL>0.5 on W1 calls.

---

### #5 — ClinVar source rows carry NO scores
5,622 of 5,796 rows (97%) have empty `REVEL`, `CADD_phred`, `SpliceAI_masked`, `PhastCons100way`, `PhyloP100way`, `FAF95_grpmax`, `AlphaMissense_pred`. We're trusting the ClinVar assertion blind, without showing analysts the supporting evidence.

**Fix**: left-join ClinVar rows back to per-gene dbNSFP TSVs in `data/exports/cp_new/{GENE}.tsv` by (gene, chr, pos, ref, alt) and populate the score columns. Already works for the inspector tool; just bake it into the merge.

---

### #6 — "Conflicting interpretations" in ClinVar slipped through
4 rows have `Conflicting interpretations` in `ClinVar_sig`. These should not be autoclassified Benign without review.

**Fix**: in `pull_clinvar_benign.py`, exclude `clnsig` containing "Conflict".

---

### #7 — AlphaMissense column shows per-transcript multi-values
107 rows display values like `.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;LP;P;LP;LP` — the dbNSFP per-transcript pipe across all annotated transcripts. Display-ugly in SeqNext and easy for an analyst to misread (which transcript's call applies?).

**Fix**: collapse to the canonical transcript's prediction (use `Ensembl_transcriptid` index to pick the matching entry).

---

### #8 — KMT2D over-represents the upload
1,237 of 5,796 rows (21.3%) are KMT2D. ClinVar has it heavily annotated and our 1-star-or-better cut lets a lot through. Worth confirming the lab wants the full set, vs. trimming KMT2D to 2-star or applying a stricter FAF floor for it.

---

### #9 — Pipeline rule W3 ran only 23 times
W3 (FAF>0.1% + REVEL<0.290 + SpliceAI≤0.1 → Likely_benign) produced 23 rows. REVEL is missense-only, so W3 misses synonymous/intronic candidates. Probably correct, but if the lab expected W3 to cover more, this is the explanation.

---

### #10 — 40,548 rows dropped between FINAL and UPLOAD
Breakdown by reason (within new-gene scope):
- 19,234 ClinVar rows: pure synonymous coding (no ROI)
- 12,761 synonymous_catalog: see #1
- 6,538 ClinVar deep intronic (outside [-15, +6])
- 2,015 pipeline rows (mostly intergenic-flagged or deep intronic)
- 1 ClinVar HTTP 429 (retry failed)
- 439 flagged:intergenic

**Note**: 7× drop ratio is expected given the aggressive synonymous/intronic filtering, but the lab should know they're shipping ~14% of what FINAL contained.

---

## 3. Sanity checks (passed)

- No duplicate (gene, chr, pos, ref, alt) keys
- No coords with >1 transcript (no fan-out)
- 99.2% have ClinVar_sig populated; the 174 missing are pipeline-only rows (expected)
- All transcripts have version suffix (`NM_XXXX.N`)
- 67 of 68 new genes have ≥1 row (only ELOC missing — see #3)

## 4. Suggested next actions

1. Investigate ELOC `flagged:intergenic` (#3).
2. Decide on 1-star handling (#2) — split file or accept.
3. Add `predictor_conflict` flag for W1 Benigns with P/LP predictors (#4).
4. Exclude ClinVar "Conflicting interpretations" (#6) and re-merge.
5. Backfill scores onto ClinVar rows from per-gene dbNSFP cache (#5).
6. Collapse AlphaMissense to canonical-transcript only (#7).

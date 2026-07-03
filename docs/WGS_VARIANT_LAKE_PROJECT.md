# Project: Whole-Genome Variant Annotation Lake (+ optional LoRA)

A step-by-step buildout plan you can execute manually. Goal: turn the ad-hoc
per-panel scripts into one queryable, tested, genome-wide annotation lake, so
classifying any gene panel (or a whole exome/genome) becomes a query — not a
week of coordinate/transcript/liftover whack-a-mole.

Starting point already in this repo:
- `lake_poc/` — working dbt-duckdb POC on 5 chr5 genes (98.8% match to hand-built
  v3.5, 15/15 tests). This is your template; every phase below extends it.
- `scripts/` — the reference logic (rules, anchoring, conservation fill) to port.

Assumptions (adjust to taste):
- Primary build = hg38; hg19 carried via liftover (you're an hg19 shop, so both
  coordinate columns stay in every table).
- Start local (DuckDB), graduate to BigQuery only when data outgrows one box.
- Clinical rule: the lake is 100% PUBLIC data (no PHI) → cloud is fine; patient
  VCFs are annotated by joining against it and never enter the lake.

---

## Phase 0 — Foundations (½ day)

**Do**
1. Adopt the POC layout as the repo standard: `raw/` (parquet sources),
   `models/{staging,intermediate,marts}`, `tests/`.
2. Pin tooling: Python 3.12 venv, `dbt-duckdb`, `duckdb`, `pyarrow`, `pyBigWig`,
   `pyliftover`. (dbt-core does NOT support 3.14 — keep a separate 3.12 venv.)
3. Decide the canonical variant keys used everywhere:
   - genomic: `(chrom, pos_hg38, ref, alt)` and `(chrom, pos_hg19, ref, alt)`
   - transcript: `(gene, nm_base, cdot)`  ← the liftover-free join key
4. Version the raw sources (record source URL + release + download date in a
   `raw/MANIFEST.md`). Reproducibility for a clinical pipeline.

**Done when** `dbt build` on the POC is green on your machine.

---

## Phase 1 — Raw layer, genome-wide (2–4 days, mostly downloads)

Load each PUBLIC source once as partitioned parquet (partition by chromosome).
None contain PHI.

| Source | What | Where | Size |
|---|---|---|---|
| dbNSFP 5.x | all coding non-syn SNVs + scores (REVEL, CADD, AlphaMissense, conservation, spliceai) | dbnsfp download | ~few hundred GB unz. |
| gnomAD v4.1 | observed variants + FAF + VEP (HGVSc, MANE_SELECT, HGVSp) | GCS `gcp-public-data--gnomad`; also a BigQuery public dataset | ~TB (VCF) |
| SpliceAI | **precomputed** all-SNV masked+unmasked scores | Illumina/Basespace precomputed VCFs | ~28 GB SNV + indels |
| UCSC conservation | phastCons100way + phyloP100way bigWig (hg38 + hg19) | hgdownload.soe.ucsc.edu | ~10–20 GB each |
| ClinVar | variant_summary.txt.gz (weekly) | NCBI | ~few hundred MB |
| RefSeq/MANE | NM_ ↔ ENST ↔ MANE mapping | NCBI MANE | small |

**Do**
1. For each: download, normalize columns to the canonical keys, write
   `raw/<source>/chrom=<N>/*.parquet`. Reuse the POC's `raw_*` extraction logic
   (`scripts/` has working parsers for dbNSFP, gnomAD tabix, SpliceAI VCF, UCSC
   bigWig).
2. Synonymous catalog: port `scripts/enumerate_synonymous_catalog.py` to run
   genome-wide over RefSeq CDS (BioPython codon table). Store as a raw source.
3. Do NOT re-run SpliceAI — load the precomputed scores. (Re-running genome-wide
   is compute-prohibitive; you only ran it locally because a few hundred panel
   variants were missing from the precomputed set.)

**Gotchas (learned the hard way this project)**
- gnomAD's displayed `c.` is MANE Select — not always your project NM_. Extract
  HGVSc for the specific NM_ via the VEP `MANE_SELECT` field.
- dbNSFP EXCLUDES synonymous → conservation must come from the UCSC track, not
  dbNSFP, for synonymous positions.
- Large tabix/bigWig range reads: read per-gene (or per-bin), never per-whole-
  chromosome-span (a 190 Mb `.values()` call will OOM).
- Every big pull needs a non-empty-result assertion — 4 gene pulls silently
  returned empty this project and dropped real variants.

**Done when** `dbt source freshness` passes and row counts match published
source totals (±expected).

---

## Phase 2 — Staging + intermediate models (2–3 days)

Port the POC models to genome scale.

**Do**
1. `stg_*` views: one per raw source, typed and cleaned (see POC).
2. `int_variants_annotated`: the join layer.
   - union coding (dbNSFP) + synonymous (catalog)
   - gnomAD **c.-anchor** on `(gene, nm_base, cdot)` → gnomad_id, g_hg38, p_hgvs,
     `in_gnomad` flag; position-fallback (liftover) only where project NM_ ≠ MANE
   - conservation backfill: `coalesce(dbnsfp.phastcons, ucsc.phastcons)` — the
     single most important join (makes synonymous classifiable)
   - SpliceAI join (masked + unmasked)
   - liftover columns: carry both hg19 and hg38 on every row
3. Materialize `int_*` as partitioned parquet/table.

**Done when** a spot-check gene reproduces the hand-built v3.5 ≥99%.

---

## Phase 3 — Marts / classification (1–2 days)

**Do**
1. `mart_upload`: the W1/W2/W3 rules as SQL `case` (see POC mart). Keep the rule
   thresholds in one dbt `var` block so the lab can tune them without editing SQL.
2. Filter models (make each a discrete, testable step):
   - intronic ROI [-20,+10]
   - 5'UTR removal (c.-N)
   - SpliceAI gate (unmasked ≤ 0.1; drop measured > 0.1)
   - nonsense exclusion
3. `mart_panel_<name>`: parameterized by a gene list → per-panel upload. cp_new
   becomes `where gene in (…)`; any future panel is just a new gene list.

**Done when** `mart_panel_cp_new` matches v3.5 within the known, documented diffs.

---

## Phase 4 — Tests & CI (1 day, highest ROI)

Every bug this project was a data-integrity failure. Encode them as tests so the
build goes red instead of the lab finding it weeks later.

**Do** — schema tests on the marts:
- `not_null(gnomad_id)` — anchored
- `not_null(phastcons)` + `accepted_range(0,1)` — conservation present/valid
- `not_5utr(cdot)` — custom generic test (in POC)
- `accepted_values(classification)` — sanctioned labels only
- `accepted_values(consequence in missense,synonymous)` — nonsense excluded
- `unique(gene||cdot)` — no dupes
- relationship test: every mart row exists in `raw_gnomad_observed`
- freshness: fail if ClinVar/gnomAD releases are stale

Wire `dbt build` into CI (GitHub Actions) on every change to models or raw
manifest.

**Done when** CI runs `dbt build` green and a deliberately-broken model fails it.

---

## Phase 5 — Scale-out to BigQuery (optional, 2–3 days)

Only when local DuckDB stops being comfortable (true genome-wide multi-source).

**Do**
1. Swap the profile to `dbt-bigquery`. Models mostly unchanged (watch DuckDB-only
   SQL: `exclude`, `any_value`, `regexp_matches` → BQ equivalents).
2. gnomAD v4.1 is already a BQ public dataset — reference it directly, no upload.
3. Load dbNSFP / SpliceAI / conservation / ClinVar into BQ once.
4. Keep patient VCFs OUT of BQ — annotate on-prem by exporting the relevant lake
   slice, or via a controlled reverse-ETL. Document the PHI boundary.

**Done when** the same `dbt build` produces identical marts on BQ.

---

## Phase 6 — Patient annotation service (2–3 days)

**Do**
1. Given a patient VCF: normalize → join against the lake on
   `(chrom,pos_hg38,ref,alt)` → return annotations + rule classification.
2. Package as a CLI (extend `scripts/inspect_variant.py`) and/or the Streamlit
   inspector already in the repo.
3. This is where the lake pays off daily: any variant, any gene, instant.

**Done when** an analyst can paste a variant and get the full annotation + call.

---

## Phase 7 — LoRA (optional, later — needs labeled data first)

Do NOT start until you have a few hundred+ lab-verified verdicts accumulated.

**Do**
1. Capture every lab review verdict (Ainslay's yellow/grey/purple decisions) as
   a raw source in the lake: `(variant, evidence snapshot, verdict, rationale)`.
   This is the gold training set — not dbNSFP.
2. Build training data from that table (the lake gives you the evidence context
   per variant for free).
3. Fine-tune (the repo already has the LoRA harness: `training/finetune_acmg.py`,
   Llama-3.2-3B). Scope the model to the AMBIGUOUS middle (VUS, conserved
   synonymous, borderline missense) — keep deterministic rules for clean
   auto-benign calls. A model there only adds non-determinism.
4. Optional RAG assist: local FAISS/DuckDB-VSS over variant embeddings +
   literature for BP4/BP7 evidence summaries. Self-hosted, on-prem (no Pinecone
   for clinical data).

**Done when** the model measurably helps on the review pile (agreement with lab
verdicts on a held-out set), used as decision-support, not authority.

---

## Effort summary

| Phase | Effort | Blocking? |
|---|---|---|
| 0 Foundations | ½ day | — |
| 1 Raw layer | 2–4 days | downloads |
| 2 Staging/int | 2–3 days | Phase 1 |
| 3 Marts | 1–2 days | Phase 2 |
| 4 Tests/CI | 1 day | Phase 3 |
| 5 BigQuery | 2–3 days | optional |
| 6 Patient service | 2–3 days | Phase 3 |
| 7 LoRA | weeks | needs labels first |

Minimum viable lake (Phases 0–4, local DuckDB): **~1.5–2 weeks**. That already
replaces the script pile and kills the rework loop. Everything after is upside.

## Guiding principles (from this project's scars)

1. **It's a join problem, not a modeling problem.** Every failure was a coord/
   transcript/staleness join. Fix joins with tested models; don't reach for ML.
2. **The `c.` is the liftover-free key.** Anchor on `(gene, NM_, c.)`; use
   position/liftover only as fallback.
3. **Conservation for synonymous comes from UCSC, not dbNSFP.**
4. **Public lake, private patients.** Clean governance boundary.
5. **Tests are the deliverable.** A green build the lab can trust beats a
   correct-today spreadsheet nobody can reproduce.

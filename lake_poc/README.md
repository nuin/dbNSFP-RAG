# cp_new variant lake — dbt-duckdb proof of concept

Single-"chromosome" POC (chr5 cp_new genes: CTNNA1, DDX41, MSH3, TERT, TRIP13)
demonstrating the data-lake + dbt architecture for the whole-genome vision.
Proves the pattern before scaling to 24 chromosomes / BigQuery.

## What it proves

The pile of `build_upload_v3*.py` + backfill/filter scripts becomes a tested,
lineage-tracked DAG. Everything that broke this week (stale copies, missing
conservation, 5'UTR leakage, unanchored variants, transcript joins) is now
either a model join or a schema test that fails the build instead of reaching
the lab.

**Result:** reproduces the hand-built v3.5 at **98.8% on MSH3 (611/618)**,
builds in <1s, **15/15 tests pass**. Remaining MSH3 diffs are the filters not
yet ported (ROI [-20,+10], SpliceAI gate, gnomad_common indel source).

## Architecture (mirrors the WGS target)

```
raw/  (precomputed parquet — the "lake", all PUBLIC data, no PHI)
  raw_dbnsfp.parquet          non-syn coding variants + scores
  raw_syn_catalog.parquet     enumerated synonymous SNVs
  raw_gnomad_observed.parquet position-level observed + FAF
  raw_gnomad_hgvsc.parquet    (gene, NM_, c.) -> gnomad_id/g./p.  (the anchor)
  raw_conservation.parquet    UCSC phastCons/phyloP per base

models/
  staging/     stg_*  — typed, cleaned views over each raw source
  intermediate/int_variants_annotated — THE join layer:
                 union coding+synonymous, gnomAD c.-anchor, conservation backfill
  marts/       mart_cp_new_upload — W1/W2/W3 rules -> the SeqNext upload
```

## Tests (the guardrails)

- `not_null(gnomad_id)` — every row must be gnomAD-anchored
- `not_null(phastcons)` + `accepted_range(0,1)` — conservation present & valid
- `not_5utr(cdot)` — custom: no c.-N 5'UTR variants
- `accepted_values(classification)` — only the sanctioned W1/W2/W3 labels
- `accepted_values(consequence in missense,synonymous)` — nonsense excluded

Each maps to a specific bug from this week's rework.

## Run

```bash
python3.12 -m venv .venv-dbt && .venv-dbt/bin/pip install dbt-duckdb
cd lake_poc && ../.venv-dbt/bin/dbt deps && ../.venv-dbt/bin/dbt build
```
Output: `lake.duckdb` (query `mart_cp_new_upload`).

## Scaling to whole genome

- **Engine**: swap the duckdb profile for **dbt-bigquery**. gnomAD v4.1 is
  already a BigQuery public dataset — join server-side, no TB download.
- **Raw layer**: load dbNSFP (~84M coding SNVs), ClinVar, SpliceAI *precomputed*
  all-SNV scores (~28GB, do NOT re-run the model), conservation bigWig once.
  Partition parquet by chromosome.
- **Public/private split**: the entire lake is public (no PHI) so it can live in
  the cloud freely; patient VCFs are annotated by joining against it, on-prem.
- **Same models, same tests** — only the sources and the profile change.

## Not yet ported (for parity with v3.5)

- ROI [-20,+10] filter and 5'UTR removal as models (currently only the mart)
- SpliceAI masked/unmasked as a source + the W3 splice gate
- gnomad_common indel source (high-FAF dels/dups not in dbNSFP/syn_catalog)
- the lab's verified verdicts as a source (future LoRA training table)

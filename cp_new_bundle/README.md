# cp_new pipeline bundle

Self-contained snapshot of the cp_new pipeline work — scripts, library
changes, documentation, and produced artifacts. Built 2026-05-13.

## What's where

| Path | Contents |
|---|---|
| `docs/cp_new_pipeline.md` | **Reference doc** — what the pipeline is, scope, workflows, variant funnel, rule rationale, deliverables, manual-review lists, known issues. 700 lines. |
| `docs/cp_new_process_journal.md` | **Process narrative** — how we got here. Chronological story of discoveries (gnomAD v4 lost VEP CSQ, dbNSFP transcript bug, SpliceAI install patches, 4-way parallel, 255→13 splice count), decisions made and why, what I'd do differently. 334 lines. |
| `scripts/` | Seven scripts that make up the pipeline (see below) |
| `src_changes/` | Modified library files (`src/chunker.py`, `src/panels.py`, `src/variantdb.py`) — reference copies; the real changes live in the parent repo |
| `bed/` | The lab's BED file — `C+_ALL_IDPE_OCT2025.bed` (166 genes, 2380 regions) |
| `outputs/` | All produced artifacts (see below) |

## Outputs at a glance

**Main upload**: `outputs/cp_new_seqnext_minimal.tsv` — 4-column SeqNext-ready
TSV. 4,567 rows. Columns: `gene, transcript, hgvs_c, classification`.

**For audit / review**:
- `outputs/cp_new_seqnext_combined.tsv` — same 4,567 rows plus QA columns
  (chr, pos, ref, alt, vv_status). Use this to investigate any specific call.
- `outputs/review/_review_canonical_splice.tsv` — **13 rows** at canonical
  splice positions (`-1`/`-2`/`+1`/`+2`). Filter artifact of SpliceAI masked
  mode; needs clinical review before upload. 6 of these are real CASR cases;
  7 also flagged as intergenic (drop under either filter).
- `outputs/review/_review_intergenic.tsv` — **70 rows** that VariantValidator
  reports as non-coding on the BED's RefSeq NM_. Mostly FANCD2 (62) +
  SMARCA4 (8). Drop or re-resolve against an alternative isoform.

**For analysis**:
- `outputs/cp_new_all_annotations.tsv.gz` — full per-variant annotations,
  **100,034 rows × 69 columns**. dbNSFP scores + gnomAD v4.1 FAF + SpliceAI
  masked + classification + VV status. Gzipped (6.2 MB).
- `outputs/cp_new_viz.db` — standalone **SQLite** for visualization tools
  (Datasette, DB Browser, Metabase). Two indexed tables: `cp_new_variants`
  (100,034 rows, wide layout) + `cp_new_classifications` (4,567 rows from
  SeqNext). 38 MB.

**Per-gene splits**:
- `outputs/per_gene_raw/{GENE}.tsv` — **raw dbNSFP + gnomAD + SpliceAI
  annotations per gene, one file per gene (165 files, all panel genes).**
  Every variant for that gene with all 59+ columns. Use for QA / debugging
  / re-classification.
- `outputs/per_gene/{GENE}_seqnext.tsv` — post-classification, 4-column
  SeqNext-ready format, one file per gene **(100 files only — genes that
  had at least one W1/W2/W3 hit).** The other 65 genes had no Benign/LB
  calls, so they only appear in `per_gene_raw/`, not here.

## Scripts (run order)

| # | Script | What it does | Runtime |
|---|---|---|---|
| 1 | `cp_new.py` | Pulls dbNSFP + gnomAD v4.1 FAF via remote tabix; emits per-gene TSVs + consolidated VCF for SpliceAI | ~30-60 min |
| 2 | `annotate_spliceai.py` | Joins externally-computed SpliceAI scores back into per-gene TSVs | seconds |
| 3 | `classify_cp_new.py` | Applies the three workflows (W1/W2/W3); emits SeqNext per-gene + combined | seconds |
| 4 | `resolve_hgvs_vv.py` | Re-anchors HGVS c. to the BED's RefSeq NM_ via VariantValidator REST | ~50-80 min |
| 5 | `extract_review_lists.py` | Extracts canonical-splice + intergenic review lists | seconds |
| 6 | `build_cp_new_exports.py` | Builds the minimal SeqNext TSV, complete annotations gzip, and viz SQLite | ~10 sec |
| — | `load_cp_new_sqlite.py` | (optional) Loads the cp_new panel into the production API SQLite database | ~15 sec |

SpliceAI itself is run separately (not a script in this bundle); see the
"SpliceAI setup" section of `docs/cp_new_pipeline.md` for instructions.

## Counts summary

```
100,034   raw dbNSFP rows across 165 cp_new genes
 99,436   unique variants
 53,984   in BED (assay-covered)
 53,574   with SpliceAI masked DS_MAX
  4,567   classified Benign / Likely_benign by the 3 rules
            20  W1   FAF >5%
         4,460  W2   synonymous + low splice + low conservation
            87  W3   rare + low REVEL + low splice
  4,497   c. resolved cleanly on BED transcript via VariantValidator
     70   flagged intergenic (drop / re-resolve)
     13   at canonical splice positions (clinical review)
```

## Workflows (full detail in `docs/cp_new_pipeline.md`)

1. **Benign** — `gnomad_v41_faf95_grpmax > 0.05`
2. **Benign** — synonymous AND SpliceAI masked ≤ 0.1 AND PhastCons < 1.0
   AND (if intronic) PhyloP < 0.1
3. **Likely_benign** — FAF > 0.001 AND REVEL < 0.290 AND SpliceAI masked ≤ 0.1

## Known issues / manual-review items

- **13 canonical-splice W2 calls** — SpliceAI masked is zero by design at
  canonical sites. These pass W2 spuriously. 6 real CASR cases need review.
- **70 intergenic flags** — variants non-coding on the BED's NM_ per
  VariantValidator. Mostly FANCD2 NM_033084.5 (62). Drop or re-resolve.
- **TSC1 c.1334-2A>G** in W1 — splice-acceptor at FAF >5%. Pop freq
  overrides splice prediction at >5% but worth confirming the FAF.

See `docs/cp_new_pipeline.md` for the full known-issues list and recovery
options.

## What this pipeline does NOT do

The pipeline pre-classifies **only** the easy-benign / easy-likely-benign
subset of variants in cp_new genes. **~49,000 in-BED variants don't match
any of the three rules** and still need traditional ACMG/AMP classification
through the project's `src/acmg_scoring.py` pipeline. cp_new short-circuits
the low-hanging cases — it doesn't replace ACMG.

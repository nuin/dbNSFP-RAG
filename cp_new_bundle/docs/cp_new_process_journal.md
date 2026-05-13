# cp_new — process journal

Chronological narrative of how the cp_new pipeline was built and the
decisions that shaped it. Companion to `cp_new_pipeline.md` (the reference
doc covering what the pipeline IS). This one covers HOW WE GOT HERE — the
discoveries, the dead-ends, what we changed mid-stream, and why.

Single session, ~6 hours of build + ~5 hours of SpliceAI compute time.

## 0. Starting position

User had:
- An existing dbNSFP RAG pipeline (FAISS vector store + SQLite backend) for
  the 314-gene `NGSgenes` panel, ~270k variants loaded.
- A new 165-gene hereditary cancer panel (`cp_new`), 67 genes net-new vs
  `NGSgenes`.
- Three lab-specific benign/LB filtering rules to apply to the panel.
- An output requirement: per-variant rows of `gene, transcript, hgvs_c,
  classification` for SeqNext upload.
- A BED file with per-region RefSeq transcript anchoring.
- Production constraint: liftover hg19↔hg38 is unreliable and not acceptable.

The three rules (paraphrased):
1. Benign: gnomAD grpmax FAF >5%.
2. Benign: synonymous variants with SpliceAI ≤0.1 (masked), PhastCons <1.0,
   and for intronic, PhyloP <0.1.
3. Likely_benign: FAF >0.1% AND REVEL <0.29 AND SpliceAI ≤0.1.

## 1. Initial data-source survey

First question: does dbNSFP alone have everything we need?

**Finding**: dbNSFP 5.3.1a carries dual hg19+hg38 coordinates per row, REVEL,
CADD, PhastCons, PhyloP, ClinVar — most of what we need. But:
- **No SpliceAI variant-level scores** (dbNSFP supports SpliceAI as an
  *attached* DB but doesn't ship the scores).
- **No grpmax FAF** — only POPMAX_AF (the older raw max, not the 95% CI
  lower bound that gnomAD v4 uses as the recommended filter).

**Decision**: dbNSFP as backbone; pull FAF + SpliceAI from elsewhere.

## 2. Panel registration

Added `cp_new` (165 genes) to `src/panels.py`. Cross-check against existing
data: 96 of 165 genes had variants already in the SQLite DB via `NGSgenes`
overlap (71,270 rows). 69 genes net-new.

## 3. First-pass extraction script (`scripts/cp_new.py`)

Initial implementation: per-gene loop, each gene scans every dbNSFP chr file.
Smoke test on TP53: **38 minutes per gene**. Projected: 100+ hours for the
full panel. Untenable.

**Fix**: chromosome-grouped pass. Read each `chr*.gz` *once*; filter all
genes on that chromosome in a single streaming scan. Built a
`gene→chromosome` map from `dbNSFP5.3_gene.gz` (40,606 genes) at startup.

Result: full 165-gene extraction in ~30-40 minutes wall-clock. 100,034
variant rows; 99,436 unique after dedup.

## 4. gnomAD v4.1 integration — the surprise

Plan was to remote-tabix the gnomAD v4.1 joint sites VCF over HTTPS (no
local download), grab FAF + SpliceAI from the VEP CSQ field per variant
region.

**Discovery #1**: gnomAD v4.1 joint sites VCF has **no VEP CSQ field at all**.
664 INFO fields, all population frequencies and QC. The v4 release split
functional annotations into separate Hail tables. So:
- FAF: still recoverable (`fafmax_faf95_max_joint`).
- SpliceAI: **gone**. Not in v4 sites VCFs. Was relying on this as the
  unmasked fallback for SpliceAI.

**Discovery #2**: gnomAD v4 INFO field naming uses `_joint` suffix
consistently (`AF_joint`, `fafmax_faf95_max_joint`, `nhomalt_joint`). My
initial code tried the unsuffixed names first — worked due to fallback
chain, but reordered for clarity.

Throughput: ~2.5 s per gene's region query. 165 genes ≈ 5-15 min added to
the dbNSFP pass. Acceptable.

## 5. SpliceAI sourcing — the decision tree

Without gnomAD's SpliceAI, three paths:
- **(a)** Use gnomAD's unmasked → not in v4, dead.
- **(b)** Download Illumina's precomputed masked VCFs (~80 GB) from
  BaseSpace.
- **(c)** Run SpliceAI locally on just the cp_new variants we extracted.

User picked (b) initially → I pointed out the BaseSpace download is
auth-walled and ~80 GB. User then picked (c).

**Rationale for (c)**: we have ~99k variants, not the entire genome. Local
SpliceAI on a small variant set is faster than downloading precomputed
files for every possible variant. `-M 1` flag enables masking on local
runs, matching the spec.

## 6. SpliceAI install pain

Three install gotchas, all undocumented in SpliceAI's README:

1. **`pkg_resources` import**: SpliceAI 1.3's `__init__.py` imports
   `pkg_resources` which setuptools 81+ removed/relocated. Pin
   `setuptools<80`.
2. **`np.fromstring` binary mode**: SpliceAI 1.3's one-hot encoder uses
   `np.fromstring(seq, np.int8)` which numpy removed years ago. One-line
   patch in `<venv>/lib/python3.11/site-packages/spliceai/utils.py`:
   ```python
   # was: return map[np.fromstring(seq, np.int8) % 5]
   return map[np.frombuffer(seq.encode("ascii") if isinstance(seq, str) else seq, np.int8) % 5]
   ```
3. **TF on M1 Ultra**: `tensorflow>=2.13` works natively, but no Metal
   acceleration without separately installing `tensorflow-metal`. Stock
   install runs on CPU.

Isolated SpliceAI into its own venv at `~/data/spliceai/venv-spliceai` so
its heavy TF dep doesn't disturb the project's `.venv`. Reference FASTA
at `~/data/hg38/hg38.fa` (1 GB compressed → 3 GB uncompressed from UCSC).

## 7. SpliceAI runtime — single-thread vs 4-way parallel

First run, single-threaded: ~110ms per inference, but each variant takes
multiple inferences. Observed throughput: **~1.4 variants/sec**. For
99,436 variants that's ~19.7 hours. Unacceptable.

Two mitigations:

**(a) BED-pre-filter the SpliceAI input.** The classifier drops variants
outside the BED anyway, so scoring them wastes time. Filter:
- BED is hg19 coords; cp_new.vcf is hg38 (dbNSFP primary). Initial filter
  attempt compared hg38 VCF against hg19 BED → 5% match, wrong axis.
- Correct filter: use per-gene TSVs' `hg19_chr` / `hg19_pos(1-based)` for
  the BED lookup, emit hg38 coords back into the VCF for SpliceAI input.
- Result: 99,436 → **53,984** unique variants (matches the classifier's
  in-BED count).

**(b) 4-way parallel SpliceAI on chunks.** Split the BED-filtered VCF into
4 equal chunks (~13,496 each), run 4 spliceai processes in parallel with
`TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2`. M1 Ultra schedules 2
chunks onto P-cores and 2 onto E-cores; first 2 finish in ~3 hr, latter 2
in ~4 hr. Total wall time ~4 hours vs ~10 hours single-thread.

Concatenate chunk outputs (header from chunk_00, data from all four).
Result: **53,574 of 53,984 (99.2%) scored**. 410 missed on chunk-split
boundaries; lived with the loss.

## 8. Classification — first run

Wrote `scripts/classify_cp_new.py`. Applies the three workflows in priority
order, per-variant per-BED-region. Edge cases:

- **Synonymous detection** via `HGVSp_snpEff`, `aaref == aaalt`, or
  `codon_degeneracy in {2, 4}`. Quirk: intronic variants have empty aaref/
  aaalt, which pandas reads as NaN, which `str(NaN).upper() == "NAN"`, so
  the `aaref == aaalt` check returns True. This is *intentional* in spec
  (W2 explicitly handles intronic via the PhyloP branch), but the side
  effect is canonical-splice variants pass the synonymous gate too.
- **Intronic detection** via `+/-` offset in HGVSc with no protein change.

First-run numbers:
- 100,034 total rows; 45,803 out-of-BED; **4,567 classified** (20 W1,
  4,460 W2, 87 W3); 58 still pending SpliceAI.

I had previously surfaced "255 W2 calls at canonical splice positions" as
a SpliceAI-masked artifact. That count came from the OLD (pre-VV) c.
notation; see §10 for what happened to it.

## 9. The transcript bug

User asked: "didn't you get the transcripts we are using?"

**Discovery**: `classify_cp_new.py` was emitting the BED's RefSeq NM_ in
the `transcript` column, paired with `HGVSc_snpEff` from dbNSFP in the
`hgvs_c` column. But dbNSFP records each variant only on whichever
Ensembl transcript snpEff chose during annotation — **often an
alternative isoform, not MANE Select**.

Concrete case: BRCA1 hg19:41197801 T>A appears in dbNSFP with
`Ensembl_transcriptid = ENST00000468300`, `HGVSc_snpEff = c.2100A>T`. But
MANE Select / `NM_007294.4` is `ENST00000357654` — not present in the
dbNSFP transcript list for this position at all (the variant might be
intronic on canonical and exonic only on the alternative).

So my SeqNext output paired `BRCA1 NM_007294.4` with `c.2100A>T` —
mixing the BED's transcript ID with a c. notation from a different
transcript. Wrong c. coordinates for an unknown number of variants.

Why this happens at extraction: `cp_new.py` dedupes by `variant_id`
(chr_pos_ref_alt), keeping whichever transcript row pandas saw first.
For genes where dbNSFP carries multiple transcripts per variant, the
surviving row is essentially arbitrary.

## 10. VariantValidator resolver

Built `scripts/resolve_hgvs_vv.py`. For every SeqNext row, query
`rest.variantvalidator.org` with chr-pos-ref-alt + the BED's RefSeq NM_.
VV returns the authoritative c. on that exact transcript. Overwrite
`hgvs_c` with the resolved value; tag rows where the variant is
non-coding on the requested transcript as `vv_status=flagged:intergenic`.

Cache results to `data/exports/cp_new/.vv_cache.json` so re-runs are
incremental.

**Rate-limiting incident**: First run, 50ms sleep per call. ~20% of calls
came back HTTP 429 (rate-limited) — VV started throttling. 4,567 calls
took 52 min, but 900 of them failed. Increased sleep to 1.2 s/call,
cleared the 429 entries from cache, retried just those 900 — completed
in 23 min, all OK on retry.

Final state: **4,497 OK**, **70 flagged:intergenic**, 0 errors. Per-gene
files updated in 3 sec (all cache hits).

**The 255 → 13 reduction**: with the corrected c. notation on the BED
transcripts, re-counting canonical-splice positions (`-1`/`-2`/`+1`/`+2`)
dropped from 255 to 13. The original 255 were artifacts of looking at c.
notation on dbNSFP's alternative transcripts whose exon structures put
genomic positions at false splice offsets. Of the 13 survivors:
- 6 real CASR canonical-splice variants (`c.1609-1G>X`, `c.1609-2A>X`)
- 7 dual-flagged as intergenic (drop under either filter)

So the genuine clinical-review subset is **6 CASR + 70 intergenic**.

## 11. Schema additions

Added 4 columns to the SQLite `variants` table:
- `gnomad_v41_faf95_grpmax` REAL
- `spliceai_ds_max_masked` REAL
- `hg38_chr` TEXT, `hg38_pos` INTEGER

The `hg38_chr/pos` columns are the lesson-learned column: the existing
DB stored only GRCh37 coords, which meant any future cross-build
annotation (like gnomAD v4.1 hg38 lookups) would have to re-extract from
dbNSFP to get hg38. Now those coords are preserved.

Idempotent ALTER TABLE migration in `VariantDatabase._init_schema()`
applied to the existing 272,875-row DB without rebuilding.

Indexes added: `idx_faf95` on `gnomad_v41_faf95_grpmax`, `idx_spliceai`
on `spliceai_ds_max_masked`.

## 12. Final SQLite load

`scripts/load_cp_new_sqlite.py` runs the per-gene TSVs through the
existing chunker (same code path as `build-sqlite`) and inserts via
`add_variants`. Reuses chunker to ensure consistent metadata extraction.

After load: DB grew from 272,875 → **303,628 rows**. `cp_new` panel
registered with 99,433 variants in `variant_panels`.

## 13. Exports for hand-off

Three deliverables built by `scripts/build_cp_new_exports.py`:

1. **`cp_new_seqnext_minimal.tsv`** — the SeqNext upload file. 4 columns,
   4,567 rows.
2. **`cp_new_all_annotations.tsv.gz`** — complete dump, 100,034 rows × 69
   columns. For analysis / QA / re-classification.
3. **`cp_new_viz.db`** — standalone SQLite with `cp_new_variants` (wide
   layout, 100,034 rows) + `cp_new_classifications` (4,567 SeqNext rows).
   Indexed for Datasette / DB Browser / Metabase.

Plus per-gene splits (100 files) and the two review files (canonical
splice, intergenic).

## Decisions made and rationale

| Decision | Rationale |
|---|---|
| dbNSFP as backbone, not just gnomAD | dbNSFP has dual hg19+hg38 coords per row, plus REVEL/CADD/conservation. gnomAD v4.1 has FAF but no functional annotation. |
| Remote tabix on gnomAD (not local download) | Joint sites VCF is ~hundreds of GB total. Remote range queries are seconds each; per-gene batching makes this trivial. |
| Local SpliceAI run, not Illumina precomputed | 80 GB precomputed download vs running on the 54k filtered variants ourselves — local is faster end-to-end for this size. |
| Isolated SpliceAI venv | TensorFlow is heavy and conflicts with project deps; isolation keeps `.venv` clean. |
| 4-way parallel SpliceAI | M1 Ultra has plenty of cores; serial run was projected at 20 hr, parallel ~5 hr. |
| BED-filter before SpliceAI | Variants outside BED can't be reported anyway; scoring them wastes compute. |
| Keep all transcripts (no canonical filter at extract) | "You don't know which transcript the sample is using" — preserve multi-transcript dbNSFP data, anchor HGVS on BED transcript at output time only. |
| VariantValidator REST for final HGVS | Authoritative source for HGVS on specific RefSeq NM_. Free, well-known in clinical workflows. |
| `hg38_chr/pos` in schema | One-time cost during this load; saves re-extracting dbNSFP for every future cross-build annotation. |
| Standalone viz SQLite | Decoupled from the production API DB; can be handed off / shared / queried without project setup. |

## What I'd do differently next time

- **Position-aware VCF splitting for SpliceAI**: line-based split caused
  410 chunk-boundary misses. SpliceAI needs ±5kb context around each
  variant; the split should ensure no variant is within 5kb of a chunk
  boundary.
- **Install `tensorflow-metal`** for M1 Ultra. Likely 5-10× speedup on
  the SpliceAI step. Not done because I didn't want to risk breaking the
  TF install mid-run.
- **Multi-transcript preservation at extract time**: `cp_new.py` dedupes
  by variant_id and silently picks one transcript. For panels where the
  BED-anchored transcript might differ from dbNSFP's chosen one, this
  creates the bug we hit with BRCA1. Better: keep all dbNSFP rows, dedupe
  later with knowledge of the BED's preferred transcript.
- **Resolve via VV during classify**, not after. The VV step happens
  after classification, which means the classifier sees stale c. when
  detecting splice positions. The 255→13 reduction is a clean example —
  if VV had resolved before the splice-position check, the 255 phantom
  count would never have appeared.
- **Add a `--no-canonical-splice` flag** to classify_cp_new.py to
  exclude variants whose c. on the BED transcript hits canonical
  splice positions. Currently those 6 CASR cases come out as W2 Benign
  and require manual review.

## Open items / next steps

1. **6 CASR canonical-splice rows**: clinical decision needed. Either
   filter out of SeqNext (treat as "needs manual review") or accept as-is.
2. **70 intergenic rows**: drop from SeqNext, or re-resolve against an
   alternative FANCD2/SMARCA4 isoform.
3. **`TSC1 c.1334-2A>G`** in W1 — splice-acceptor at FAF >5%. Validate
   the FAF isn't artifact before upload.
4. **The ~49,000 in-BED, unclassified variants** — these need traditional
   ACMG/AMP classification via `src/acmg_scoring.py`. Out of scope for
   this pipeline.

## File inventory

Everything produced by this work:

| Path | Role |
|---|---|
| `src/panels.py` | added cp_new panel definition |
| `src/chunker.py` | added hg38_chr/hg38_pos preservation |
| `src/variantdb.py` | schema migrations + 4 new columns |
| `scripts/cp_new.py` | dbNSFP + gnomAD extractor |
| `scripts/annotate_spliceai.py` | join SpliceAI back into TSVs |
| `scripts/classify_cp_new.py` | apply the three workflows |
| `scripts/resolve_hgvs_vv.py` | re-anchor c. via VariantValidator |
| `scripts/extract_review_lists.py` | dump review TSVs |
| `scripts/build_cp_new_exports.py` | minimal SeqNext + all-annotations + viz DB |
| `scripts/load_cp_new_sqlite.py` | load into production SQLite |
| `docs/cp_new_pipeline.md` | reference documentation |
| `docs/cp_new_process_journal.md` | this file — process narrative |
| `cp_new_bundle/` | self-contained hand-off bundle |

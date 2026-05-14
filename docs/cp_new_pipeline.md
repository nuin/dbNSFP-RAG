# cp_new pipeline

End-to-end annotation + classification pipeline for the 165-gene hereditary
cancer panel `cp_new`. Pulls dbNSFP scores + gnomAD v4.1 FAF + (locally-run)
SpliceAI masked scores, applies three lab-specific benign/LB rules, and emits
per-gene TSVs ready for SeqNext import.

## Scope

- **Panel**: `cp_new` — 165 genes (98 overlap with `NGSgenes`, 67 are new).
  Defined in `src/panels.py`.
- **BED**: `/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed`
  (165 unique genes, 2379 regions). The BED commits one RefSeq transcript per
  region; APC and RAD51D carry two each (alternative first exon /
  alternative E3 respectively). The previous OCT2025 BED had a stray "CHR2"
  entry which the APR2026 version removes.

## Deliverables

Built by `scripts/build_cp_new_exports.py` after the full pipeline runs.

| File | Purpose | Size | Rows |
|---|---|---|---|
| `data/exports/cp_new/cp_new_seqnext_minimal.tsv` | **Main SeqNext upload file.** 4 columns: `gene, transcript, hgvs_c, classification`. Header on row 1. | 347 KB | 4,567 |
| `data/exports/cp_new/seqnext/cp_new_seqnext.tsv` | Combined SeqNext with QA columns appended (chr/pos/ref/alt/vv_status) — for audit | ~600 KB | 4,567 |
| `data/exports/cp_new/seqnext/{GENE}_seqnext.tsv` (100 files) | Per-gene SeqNext, same layout as combined | — | varies |
| `data/exports/cp_new/seqnext/_review_canonical_splice.tsv` | 13 rows at canonical splice positions needing manual review | — | 13 |
| `data/exports/cp_new/seqnext/_review_intergenic.tsv` | 70 rows VV flagged as non-coding on BED transcript | — | 70 |
| `data/exports/cp_new/cp_new_all_annotations.tsv.gz` | **Complete annotation dump** — every dbNSFP/gnomAD/SpliceAI column for every variant in the 165 cp_new genes (69 columns). For analysis / QA / re-classification | 6.2 MB | 100,034 |
| `data/exports/cp_new/cp_new_viz.db` | **Standalone SQLite for visualization tools** (Datasette, DB Browser, Metabase). Two indexed tables: `cp_new_variants` (100,034 rows wide layout) + `cp_new_classifications` (4,567 SeqNext rows) | 38 MB | — |
| `data/sqlite/grch37-all-panels.db` | Production API database (NGSgenes + cp_new + others, panel-tagged via `variant_panels` table) | ~430 MB | 303,628 |

Build all three of the new exports:

```bash
uv run python scripts/build_cp_new_exports.py            # builds all three
uv run python scripts/build_cp_new_exports.py --drop-flagged   # excludes vv_status != ok from SeqNext-minimal
uv run python scripts/build_cp_new_exports.py --only seqnext   # just the minimal TSV
```

### Visualization DB quick examples

```bash
# Browse with Datasette
pip install datasette
datasette data/exports/cp_new/cp_new_viz.db

# Or sqlite3 CLI
sqlite3 data/exports/cp_new/cp_new_viz.db
```

```sql
-- top genes by classified-variant count
SELECT gene, COUNT(*) AS n
FROM cp_new_classifications
GROUP BY gene ORDER BY n DESC LIMIT 10;

-- variants in BRCA1 with their dbNSFP scores
SELECT genename, hg19_chr || ':' || hg19_pos_1_based AS pos,
       ref, alt, CADD_phred, REVEL_score, AlphaMissense_pred,
       gnomad_v41_faf95_grpmax, spliceai_ds_max_masked
FROM cp_new_variants
WHERE genename = 'BRCA1' AND CADD_phred > 20
ORDER BY CADD_phred DESC LIMIT 20;

-- everything flagged by VariantValidator
SELECT gene, transcript, hgvs_c, classification, vv_status
FROM cp_new_classifications
WHERE vv_status != 'ok';
```

## The three workflows

Priority order, first match wins. All applied per (variant, BED-region)
combination.

| # | Class | Rule |
|---|-------|------|
| 1 | Benign | `gnomad_v41_faf95_grpmax > 0.05` |
| 2 | Benign | synonymous AND `spliceai_ds_max_masked <= 0.1` AND `phastCons100way_vertebrate < 1.0` AND (if intronic) `phyloP100way_vertebrate < 0.1` |
| 3 | Likely_benign | `gnomad_v41_faf95_grpmax > 0.001` AND `REVEL_score < 0.290` AND `spliceai_ds_max_masked <= 0.1` |

Variants outside the BED are silently dropped (the assay doesn't cover them).
Synonymous detection: `HGVSp` shows `p.X=`, or `aaref == aaalt`, or
`codon_degeneracy in {2, 4}`. Intronic detection: HGVSc carries a `+/-` offset
and no protein change.

**SpliceAI handling**: a missing `spliceai_ds_max_masked` is treated as **0**
(no detected splice signal). SpliceAI v1.3 silently skips variants outside
its GENCODE V24 canonical gene model — notably the BED's alternative-
transcript regions (e.g. APC E01 on `NM_001127511.3`, RAD51D alt-E3 on
`NM_001142571.2`). Without this rule, those variants would never reach W2/W3
because they have no SpliceAI value. Variants where SpliceAI returned a
non-zero score keep their measured value. **Side effect**: canonical splice
positions (`-1`/`-2`/`+1`/`+2`) get auto-classified by W2 if synonymous-
flagged, because SpliceAI masked is intrinsically near-zero there — see
"Manual-review lists" for the 265-row review pile this generates.

## Variant funnel — what happens to every variant

This is the full accounting from raw dbNSFP rows down to the SeqNext export.
Numbers from the 2026-05-12 production run.

```
  STAGE                                                       COUNT
  ─────────────────────────────────────────────────────────────────
  Raw dbNSFP rows across 165 cp_new genes                   100,034
                                                                 │
                                                                 │ (1)
                                                                 ▼
  Unique variants (dedup on chr,pos,ref,alt)                 99,436
                                                                 │
                                                                 │ (2)
                                                                 ▼
  Unique variants within at least one BED region             53,984
                                                                 │
                                                                 │ (3)
                                                                 ▼
  Successfully scored by SpliceAI (-M 1, hg38)               53,574
                                                                 │
                                                                 │ (4)
                                                                 ▼
  Classifier emitted a Benign / Likely_benign call            4,567
                                                                 │
                                                                 │ (5)
                                                                 ▼
  vv_status = ok in SeqNext output                            4,497
  vv_status = flagged:intergenic                                 70
                                                                 │
                                                                 │ (6)
                                                                 ▼
  CASR canonical-splice subset that warrants W2 re-review         6
  (other 7 splice-pattern matches drop under intergenic filter)
```

What each transition does:

**(1) dbNSFP → unique variants.** dbNSFP records each variant once *per
transcript* it's annotated against; same `chr,pos,ref,alt` can appear on
multiple rows with different `Ensembl_transcriptid` values. `cp_new.py`
keeps the first row per `variant_id`. Loses **598 duplicate transcript
rows**. This dedup is also the source of the "wrong c. notation" bug
that VV resolves later: the transcript that survives dedup is arbitrary,
and often isn't MANE Select.

**(2) Unique variants → in-BED.** The dbNSFP per-gene extraction captures
*every nsSNV in the gene* including intronic and UTR positions. The BED
encodes the lab's actual coverage — typically coding exons + flanking
splice region. **45,452 variants fall outside the assay** and can't be
reported on, so they're dropped at classify time. Those positions are
real biological variants but the assay literally doesn't sequence them.

**(3) In-BED → SpliceAI-scored.** Locally-run SpliceAI v1.3 (`-M 1`) on
the 53,984 in-BED variants. **410 missed** because they sat on chunk-split
boundaries in the 4-way parallel run (the split was naive line-based, not
position-aware; a variant straddling two chunks gets seen by one chunk's
context window but its score lookup ends up empty). Acceptable error rate
(<1%); those variants stay in W2/W3 "pending" and don't get classified.

**(4) SpliceAI-scored → classified.** The three rules (FAF >5%, synonymous
+ low splice + low conservation, rare + low REVEL + low splice) are
*negative-prediction* filters — they identify variants we can confidently
exclude as causes of disease. Most variants don't match any of them
because they aren't common-population (FAF >5%), aren't synonymous, or do
have REVEL/SpliceAI signal. **~49,000 in-BED variants don't trigger any
benign/LB call** and stay in the "need separate classification" bucket
(see "What's NOT in the SeqNext export" below).

**(5) Classified → VV-resolved.** Every classified row gets its c.
notation re-anchored to the BED's RefSeq NM_ via VariantValidator REST.
**4,497** resolved cleanly; **70** came back `flagged:intergenic` (variant
non-coding on the BED's NM_, c. shown is from a different transcript and
shouldn't be uploaded as-is). No errors.

**(6) Manual-review subsets.** Two flags surface for human review:
- 13 rows whose c. on the BED's transcript sits at canonical splice
  `-1`/`-2`/`+1`/`+2` positions. SpliceAI masked is zero by construction
  at canonical sites, so W2 fires inappropriately. Of these:
    - **6 CASR rows** at `c.1609-1G>X` / `c.1609-2A>X` — real review cases
    - 7 dual-flagged as intergenic (drop under the intergenic filter)
- 70 intergenic rows total (62 FANCD2, 8 SMARCA4) — drop or re-resolve
  against an alternative isoform.

## Workflow rule rationale

Each of the three workflows targets a *negative* prediction — "we have
enough evidence to exclude pathogenicity." None of these calls are
positive evidence for benignity in the ACMG sense; they're lab-specific
filtering rules that pre-classify obvious benign/LB before traditional
ACMG criteria evaluation.

### Workflow 1 — Benign (FAF >5%)
```
gnomad_v41_faf95_grpmax > 0.05
```
**Equivalent ACMG criterion**: BA1 (Benign Stand-Alone). ACMG/AMP 2015 sets
BA1 at MAF >5% in any population. The gnomAD v4.1 **grpmax FAF95** (filtering
allele frequency at 95% confidence interval lower bound, max across
genetic ancestry groups) is what ClinGen recommends — more conservative
than raw POPMAX_AF and accounts for sampling noise.

A variant present in >5% of *any* well-sampled population is essentially
incompatible with a highly penetrant Mendelian disease (the disease would
be too common). This is the cleanest benign call you can make from
population data alone, with one caveat: founder populations or
late-onset / incomplete-penetrance disease can violate it. Examples
already flagged (e.g. `TSC1 c.1334-2A>G`) deserve manual sanity-check.

**Hits in production**: 20 variants. Most concentrated in SLX4 (11) — a
known highly-polymorphic locus.

### Workflow 2 — Benign (synonymous + low splice + low conservation)
```
synonymous AND
spliceai_ds_max_masked <= 0.1 AND
phastCons100way_vertebrate < 1.0 AND
(if intronic) phyloP100way_vertebrate < 0.1
```
**Equivalent ACMG criterion**: BP4 (computational evidence supports a
benign effect) + BP7 (synonymous nucleotide change with no predicted
splice impact and no high conservation). This rule combines them into a
single auto-call.

Logic:
- **Synonymous**: variant doesn't change the amino acid (or is intronic
  near a splice site — see "synonymous detection" caveat below). No
  protein-level consequence.
- **SpliceAI masked DS_MAX ≤ 0.1**: no predicted splice impact. The
  masked threshold of 0.1 is more permissive than the unmasked 0.2 that
  Illumina recommends, because masked scores skew lower (they suppress
  signal at canonical splice sites — see "Known issues" below).
- **PhastCons < 1.0**: not in a perfectly-conserved position. PhastCons
  ranges 0-1 with 1 = invariant across 100 vertebrate species; <1 means
  the position has at least some tolerated variation across evolution.
- **PhyloP < 0.1 if intronic**: extra conservation filter for
  intronic/flanking positions, since "synonymous" loses meaning there.
  PhyloP <0.1 means the position evolves at ~neutral rate.

The "if intronic" branch is critical — without it, intronic variants
incorrectly flagged as "synonymous" by the heuristic would pass W2 too
easily. The PhyloP check adds rigor at those positions.

**Hits in production**: 4,460 variants. Dominant rule.

### Workflow 3 — Likely_benign (rare + low REVEL + low splice)
```
gnomad_v41_faf95_grpmax > 0.001 AND
REVEL_score < 0.290 AND
spliceai_ds_max_masked <= 0.1
```
**Equivalent ACMG criterion**: BS1 (MAF higher than expected for disorder)
+ BP4 (computational evidence benign). The threshold combination is
calibrated for the lab's expected disease prevalence — at FAF >0.1% the
variant is rare enough to be plausibly pathogenic in principle, but
combined with low REVEL (no missense damage prediction) and low SpliceAI
(no splice damage prediction) the overall evidence points away from
pathogenicity.

REVEL < 0.290 cutoff: ClinGen SVI work has calibrated REVEL thresholds;
<0.290 is in the "supporting benign" band. Above 0.644 is "supporting
pathogenic"; the 0.290-0.644 middle is non-informative.

**Hits in production**: 87 variants.

## What's NOT in the SeqNext export

The export is **only** the variants that match W1/W2/W3. Everything else
is unclassified by this pipeline:

| Bucket | Count | Why excluded | Where to look |
|---|---|---|---|
| Outside BED | 45,452 | Assay doesn't sequence those positions | n/a (can't report) |
| In BED but no benign rule fires | ~49,000 | Need separate ACMG workflow (VUS, possibly pathogenic) | `src/acmg_scoring.py`, SQLite DB |
| Pending SpliceAI (chunk boundaries) | 58 | Lost in 4-way parallel SpliceAI split | re-run SpliceAI on missing subset |
| Multi-transcript dbNSFP duplicates | 598 | First-row dedup discards them | n/a (collapsed at extraction) |

The **~49,000 in-BED, unclassified variants** are the meat of clinical
interpretation work. They need the full ACMG/AMP criteria evaluation —
not the lab's negative-prediction filters. Those are queryable in the
SQLite DB:

```sql
SELECT v.variant_id, v.gene, v.cadd_phred, v.revel_score, v.clinvar_sig,
       v.gnomad_v41_faf95_grpmax, v.spliceai_ds_max_masked
FROM variants v
JOIN variant_panels p ON v.variant_id = p.variant_id
WHERE p.panel_name = 'cp_new'
  AND v.variant_id NOT IN (
    SELECT variant_id FROM ...  -- the 4,567 classified
  );
```

The cp_new pipeline pre-classifies the easy-benign subset so analysts
focus on the harder ~49k. It doesn't *replace* ACMG — it short-circuits
the low-hanging cases.

## Data provenance per output column

| Column | Source | Computed in |
|---|---|---|
| `gene` | dbNSFP `genename` | cp_new.py |
| `transcript` | BED column 6 (RefSeq NM_) | classify_cp_new.py (interval lookup against BED) |
| `hgvs_c` | VariantValidator REST resolution against `transcript` | resolve_hgvs_vv.py |
| `classification` | rule application | classify_cp_new.py |
| `chr_grch37` | dbNSFP `hg19_chr` | cp_new.py |
| `pos_grch37` | dbNSFP `hg19_pos(1-based)` | cp_new.py |
| `ref`, `alt` | dbNSFP `ref` / `alt` | cp_new.py |
| `vv_status` | VariantValidator response category | resolve_hgvs_vv.py |

Internal annotations used by the classifier (present in per-gene TSVs but
not in SeqNext export):

| Column | Source | Used by |
|---|---|---|
| `gnomad_v41_faf95_grpmax` | gnomAD v4.1 joint sites VCF, `fafmax_faf95_max_joint` INFO field | W1, W3 |
| `gnomad_v41_af_joint` | gnomAD v4.1 joint sites VCF, `AF_joint` INFO field | (not used in rules; reference) |
| `gnomad_v41_filter` | gnomAD v4.1 joint sites VCF, FILTER column | (not used; reference) |
| `spliceai_ds_max_masked` | SpliceAI v1.3 with `-M 1` flag | W2, W3 |
| `phastCons100way_vertebrate` | dbNSFP | W2 |
| `phyloP100way_vertebrate` | dbNSFP | W2 (intronic branch) |
| `REVEL_score` | dbNSFP | W3 |
| `HGVSp_snpEff`, `aaref`, `aaalt`, `codon_degeneracy` | dbNSFP | W2 (synonymous detection) |

## Synonymous detection — edge case

The classifier's `is_synonymous(row)` heuristic returns True when *any*
of these holds:

1. `HGVSp_snpEff` contains `p.=` (explicit synonymous marker), **or**
2. `aaref == aaalt` (same amino acid before/after) and not `X`, **or**
3. `codon_degeneracy in {2, 4}` (synonymous degeneracy classes per dbNSFP)

There's a known quirk: when a variant is **intronic** (no protein change),
dbNSFP records empty/NaN values for `aaref` and `aaalt`. pandas reads
these as `NaN`, and `str(NaN).upper() == "NAN"` — so `aaref == aaalt`
returns True (both "NAN"), and `is_synonymous` returns True for intronic
variants too. This is *not* a bug per se: the user's workflow 2 spec
explicitly mentions "intronic variants" via the PhyloP branch, suggesting
the intent was to filter both synonymous coding AND intronic positions
through the same rule. The intronic-aware PhyloP check (line 3 of the W2
rule) is the safety net.

The downstream consequence: variants at canonical splice acceptor/donor
positions (`-1`, `-2`, `+1`, `+2` HGVS offsets) qualify as "intronic"
under the heuristic, pass the PhyloP check (typically `< 0.1` because
SpliceAI training data may not cover these well), pass SpliceAI masked
≤ 0.1 (zeroed by design at canonical sites), and end up classified as
Benign. The 13 such calls in production are documented under
"Manual-review lists" — they need clinical review, not auto-classification.

## Pipeline stages

```
dbNSFP chr*.gz  ─┐
gnomAD v4.1 (remote tabix over HTTPS) ─┤
                                       ▼
                       scripts/cp_new.py
                                       │
                                       ▼
              data/exports/cp_new/{GENE}.tsv (165 files)
                       │            +
                       │   data/exports/cp_new/cp_new.vcf (99,436 unique variants)
                       │
              spliceai -M 1 (external, 4-way parallel) ─► chunk*.spliceai.vcf
                       │                                     │
                       ▼                                     ▼
          scripts/annotate_spliceai.py joins back into per-gene TSVs
                       │
                       ▼
          scripts/classify_cp_new.py  ─►  seqnext/{GENE}_seqnext.tsv
                       │
                       ▼
          scripts/resolve_hgvs_vv.py  ─► (overwrites c. via VariantValidator)
                       │
                       ▼
          scripts/load_cp_new_sqlite.py
                       │
                       ▼
          data/sqlite/grch37-all-panels.db
          (cp_new variants + variant_panels(panel='cp_new'))
```

## Scripts

### `scripts/cp_new.py`
Pulls dbNSFP + gnomAD v4.1 joint sites for cp_new genes.

- Reads each `dbNSFP5.3.1a_variant.chr{N}.gz` **once**, filters all cp_new
  genes on that chromosome in one streaming pass (chromosome-grouped — adding
  more genes is free).
- Uses dbNSFP gene file (`dbNSFP5.3_gene.gz`) to build `gene→chromosome` map
  at startup, so each gene maps to exactly one chr file.
- For each gene's variant block, runs **one remote tabix call** against
  `gs://gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{N}.vcf.bgz`
  to pull all variants in that gene's hg38 interval, then joins on
  `(#chr, pos, ref, alt)`.
- Emits a consolidated VCF with `--emit-vcf` for the SpliceAI step.
- Resumable: skips per-gene TSVs that already exist non-empty (use `--force`
  to redo).

```bash
uv run python scripts/cp_new.py --emit-vcf data/exports/cp_new/cp_new.vcf
```

Output columns include both GRCh38 (`#chr`, `pos(1-based)`) and GRCh37
(`hg19_chr`, `hg19_pos(1-based)`) — no liftover ever needed. Appended at the
end: `gnomad_v41_faf95_grpmax`, `gnomad_v41_af_joint`, `gnomad_v41_grpmax_anc`,
`gnomad_v41_nhomalt`, `gnomad_v41_filter`, `spliceai_ds_max_masked` (empty
until step 3).

### `scripts/annotate_spliceai.py`
Reads the SpliceAI-output VCF and joins masked DS_MAX into the per-gene TSVs
in place.

```bash
uv run python scripts/annotate_spliceai.py \
    --spliceai-vcf data/exports/cp_new/cp_new.spliceai.vcf \
    --tsv-dir data/exports/cp_new
```

### `scripts/classify_cp_new.py`
Applies the three workflows. Looks up each variant against the BED to find
its region's transcript, then emits SeqNext rows.

```bash
uv run python scripts/classify_cp_new.py
```

`workflow_pending` counters tell you how many variants are candidates for
W2/W3 but waiting on SpliceAI — i.e., expected delta after re-running once
SpliceAI annotation is filled in.

### `scripts/resolve_hgvs_vv.py`
**Critical correction step.** The classifier emits the BED's RefSeq NM_ as
the transcript column and dbNSFP's `HGVSc_snpEff` value as `hgvs_c`, but
dbNSFP's c. notation is on whatever Ensembl transcript dbNSFP chose for that
row — which is *not necessarily* the BED's transcript. For many genes (BRCA1
being the obvious case) dbNSFP's chosen transcript isn't even MANE Select.
Without this step, the SeqNext rows would pair the BED's NM_ with a c. on
a different transcript — wrong c. coordinates.

This script queries `rest.variantvalidator.org` per row:

```
GET /VariantValidator/variantvalidator/GRCh37/17-41197801-T-A/NM_007294.4
   -> top-level key "NM_007294.4:c.5486A>T"
```

Replaces `hgvs_c` with the resolved value. Adds a `vv_status` column:
`ok`, `flagged:intergenic` (variant non-coding on requested transcript),
`flagged:http_429` (rate-limited), etc. Flagged rows keep their original
c. with a `(VV_FLAGGED:reason)` suffix.

Cache at `data/exports/cp_new/.vv_cache.json` makes re-runs incremental.

```bash
# combined file first (one tqdm tick per file -- silent during processing)
uv run python scripts/resolve_hgvs_vv.py --in-place --only-combined

# then propagate to per-gene files (mostly cache hits, seconds)
uv run python scripts/resolve_hgvs_vv.py --in-place
```

**Rate-limiting**: VV is hosted by University of Leicester. ~750 ms/call
typical latency. Without throttling, ~20% of calls came back HTTP 429
(rate-limited). The script sleeps **1.2 s after every API hit** (cache hits
skip). First full pass for the combined file takes ~52 min for 4,567
variants. Retry pass for 429s adds ~30 min. Future runs are fully cached
and finish in seconds.

### `scripts/load_cp_new_sqlite.py`
Loads all 165 per-gene TSVs into `data/sqlite/grch37-all-panels.db` via the
chunker (same code path as `build-sqlite`), then registers `cp_new` panel
membership. `INSERT OR REPLACE` overwrites existing rows for the 96 genes
that overlap with NGSgenes — picking up the new gnomAD/SpliceAI/hg38 columns.

```bash
uv run python scripts/load_cp_new_sqlite.py
```

## Schema changes to `src/variantdb.py`

Added four columns to the `variants` table:

| Column | Type | Purpose |
|---|---|---|
| `gnomad_v41_faf95_grpmax` | REAL | gnomAD v4.1 joint `fafmax_faf95_max_joint` |
| `spliceai_ds_max_masked` | REAL | max delta score from `spliceai -M 1` |
| `hg38_chr` | TEXT | GRCh38 chromosome — preserved going forward so cross-build annotations don't have to re-touch dbNSFP |
| `hg38_pos` | INTEGER | GRCh38 position |

A `MIGRATIONS` list at the top of `VariantDatabase` defines these for forward
migration — `_init_schema` runs `PRAGMA table_info(variants)` and adds any
column that isn't already present. The existing 272k-row DB was migrated in
place without rebuilding.

Indexes added: `idx_faf95` on `gnomad_v41_faf95_grpmax`, `idx_spliceai` on
`spliceai_ds_max_masked`.

## SpliceAI setup

SpliceAI is not in the project's `.venv`. Installed in a **separate** venv
to keep TensorFlow out of the main project deps:

```bash
mkdir -p ~/data/spliceai
uv venv --python python3.11 ~/data/spliceai/venv-spliceai
~/data/spliceai/venv-spliceai/bin/python -m ensurepip
~/data/spliceai/venv-spliceai/bin/python -m pip install tensorflow spliceai 'setuptools<80'
```

Three gotchas worth noting:

1. **`setuptools<80`** — SpliceAI 1.3's `__init__.py` imports `pkg_resources`,
   which setuptools 81+ removed/relocated. Pin setuptools below 80.
2. **`np.fromstring` patch** — SpliceAI 1.3 uses `np.fromstring(seq, np.int8)`,
   which was removed in modern numpy. One-line patch in
   `<venv>/lib/python3.11/site-packages/spliceai/utils.py`:
   ```python
   # was: return map[np.fromstring(seq, np.int8) % 5]
   return map[np.frombuffer(seq.encode("ascii") if isinstance(seq, str) else seq, np.int8) % 5]
   ```
3. **No Metal acceleration** in stock TF install on M1 Ultra here — runs
   ~110 ms/variant on CPU, ~3 hr for 99k variants. `tensorflow-metal` plugin
   not installed; not blocking but a future speedup target.

Reference FASTA at `~/data/hg38/hg38.fa` (3.0 GB uncompressed + `.fai`),
downloaded from UCSC `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`.

Run:
```bash
~/data/spliceai/venv-spliceai/bin/spliceai \
    -I data/exports/cp_new/cp_new.vcf \
    -O data/exports/cp_new/cp_new.spliceai.vcf \
    -R ~/data/hg38/hg38.fa \
    -A grch38 \
    -M 1
```

## Why the BED's transcript matters

For 164 of 166 genes, the BED's RefSeq NM_ is the **MANE Select** transcript
and pairs 1:1 with an Ensembl ENST. CDS coordinates *typically* match, but
**not always**: dbNSFP records each variant on whichever transcript snpEff
chose (often the longest coding isoform or first ENST hit) — which is **not
necessarily MANE Select**. Two failure modes:

1. dbNSFP picked an alternative isoform: e.g. for BRCA1 hg19:41197801 T>A,
   dbNSFP shows `Ensembl_transcriptid = ENST00000468300`, `HGVSc_snpEff =
   c.2100A>T`. MANE Select (`ENST00000357654` / `NM_007294.4`) gives c.5486A>T.
   **MANE Select isn't even in dbNSFP's transcript list for this position**
   — probably the variant is intronic on the canonical and the alternative
   isoform is the only one where it's coding.
2. APC/RAD51D dual-tx regions: the BED explicitly uses non-canonical
   transcripts for specific regions (APC E01 on `NM_001127511.3`, RAD51D
   alt-E3 on `NM_001142571.2`). These differ structurally from canonical.

**Fix**: `scripts/resolve_hgvs_vv.py` queries VariantValidator per row to
get the authoritative c. on the BED's RefSeq NM_. Rows where VV says the
variant is non-coding on the requested transcript get `vv_status =
flagged:intergenic` and keep their original c. with a `(VV_FLAGGED:intergenic)`
suffix — manual review or drop before SeqNext upload.

## Data sources and what's missing

| Source | Has | Doesn't have |
|---|---|---|
| dbNSFP 5.3.1a | REVEL, CADD, AlphaMissense, ClinVar, gnomAD AF (joint), PhastCons, PhyloP, conservation, gene constraint, hg19+hg38 coords per variant | SpliceAI scores (variant-level); FAF grpmax |
| gnomAD v4.1 joint sites VCF | population freqs, FAF (`fafmax_faf95_max_joint`), filter status, popmax | **no VEP CSQ, no SpliceAI** — v4 split annotations out of the sites VCF |
| Illumina SpliceAI v1.3 (local run) | masked DS_MAX per variant | requires local install + 3 GB FASTA |

The gnomAD v4 VCF carrying no SpliceAI / no VEP CSQ was a late discovery —
earlier discussions assumed v4 retained v2/v3's VEP-annotated form. It
doesn't. SpliceAI must be sourced separately (we run it locally with `-M 1`).

## Re-running from scratch

```bash
# 0. one-time setup (already done):
#    - SpliceAI venv at ~/data/spliceai/venv-spliceai
#    - hg38.fa at ~/data/hg38/hg38.fa
#    - dbNSFP 5.3.1a at ~/Downloads/dbNSFP5.3.1a/

# 1. dbNSFP + gnomAD (~30-60 min)
uv run python scripts/cp_new.py \
    --emit-vcf data/exports/cp_new/cp_new.vcf \
    --force

# 2. SpliceAI -M 1 (~3 hr on CPU)
~/data/spliceai/venv-spliceai/bin/spliceai \
    -I data/exports/cp_new/cp_new.vcf \
    -O data/exports/cp_new/cp_new.spliceai.vcf \
    -R ~/data/hg38/hg38.fa \
    -A grch38 -M 1

# 3. join SpliceAI back (~seconds)
uv run python scripts/annotate_spliceai.py \
    --spliceai-vcf data/exports/cp_new/cp_new.spliceai.vcf \
    --tsv-dir data/exports/cp_new

# 4. classify and emit SeqNext rows (~seconds)
uv run python scripts/classify_cp_new.py

# 5. resolve c. on the BED transcript via VariantValidator (~50-80 min
#    first run, seconds on cached re-runs)
uv run python scripts/resolve_hgvs_vv.py --in-place --only-combined
uv run python scripts/resolve_hgvs_vv.py --in-place

# 6. load into SQLite (~15 sec)
uv run python scripts/load_cp_new_sqlite.py
```

## Performance notes for SpliceAI

SpliceAI runs at **~1.4 variants/sec single-process** on M1 Ultra CPU
(no Metal). For the 99k full-VCF run that's ~20 hr — impractical. Two
mitigations applied during the first production run:

1. **BED-filter the input VCF** before SpliceAI. Variants outside the BED
   regions are silently dropped by the classifier anyway, so scoring them
   wastes time. Filter reduces 99,436 → 53,984 variants. The BED is hg19
   coords; the cp_new.vcf emitted by cp_new.py is hg38 — filter using the
   per-gene TSVs' `hg19_chr` / `hg19_pos(1-based)` columns and emit hg38
   coords back into the VCF.
2. **4-way parallel SpliceAI** on chunks of the BED-filtered VCF. Each
   process with `TF_NUM_INTRAOP_THREADS=4` `TF_NUM_INTEROP_THREADS=2`. Wall
   time for the full 53,984 variants ≈ 3.5 hr (each chunk ~13.5k variants,
   chunks 00/01 finish in ~3 hr, chunks 02/03 take 30-60 min longer because
   the OS schedules them onto efficiency cores).

Concatenate the 4 chunk outputs (header from chunk 00, data from all four):

```bash
cd data/exports/cp_new
grep '^#' cp_new.bed.chunk_00.spliceai.vcf > cp_new.spliceai.vcf
for i in 00 01 02 03; do
    grep -v '^#' cp_new.bed.chunk_$i.spliceai.vcf >> cp_new.spliceai.vcf
done
```

## Snapshot — production run (2026-05-13, with missing-SpliceAI-as-zero)

- **dbNSFP scan**: 100,034 variant rows across 165 genes, 8,253 with gnomAD
  FAF (8%). 99,436 unique variant_ids emitted to cp_new.vcf.
- **BED-filter for SpliceAI**: 53,984 BED-covered unique variants.
- **SpliceAI run**: 53,574 of 53,984 (99.2%) got DS_MAX scores. Remaining
  410 either fell on chunk-split boundaries or were silently skipped by
  SpliceAI (outside its GENCODE V24 gene model — e.g. APC alt-E01).
- **Classifier output**: **4,610 SeqNext rows** across 130 genes.
  - W1 (Benign FAF >5%): **20**
  - W2 (Benign synonymous + SpliceAI≤0.1 + conservation low): **4,501**
  - W3 (Likely_benign FAF >0.1% + REVEL<0.29 + SpliceAI≤0.1): **89**
- **VV resolver**: 4,540 OK, 70 flagged:intergenic, 0 errors.
- **Review piles**:
  - 265 canonical-splice review rows (consequence of missing-SpliceAI-as-zero)
  - 70 intergenic flags (variants non-coding on BED's NM_)
- **SQLite DB**: 272,875 → 303,628 rows. `cp_new` panel = 99,433 variants.

### Earlier snapshot (2026-05-12, before missing-SpliceAI-as-zero)
- 4,567 classified. 58 candidates "pending" SpliceAI. 13 canonical-splice
  review rows. Fixing the missing-SpliceAI gap moved +43 net to classified
  (APC alt-E01 alone contributed 32) and expanded the canonical-splice
  review list 13 → 265.

## Manual-review lists

`scripts/extract_review_lists.py` derives two TSVs from the combined SeqNext
output for items needing human review before SeqNext upload:

```bash
uv run python scripts/extract_review_lists.py
```

### `_review_canonical_splice.tsv` — 265 rows

W2 'Benign synonymous' calls whose c. (on the BED's RefSeq NM_, after VV
resolution) sits at a canonical splice acceptor or donor position
(`-1` / `-2` / `+1` / `+2`).

Spread across many genes; top contributors:

| Gene | Count | | Gene | Count |
|---|---|---|---|---|
| MSH2 | 6 | | RAD51D | 3 |
| TRIM28 | 6 | | ANKRD26 | 3 |
| FANCL | 6 | | VHL | 3 |
| CASR | 6 | | UBE2T | 3 |
| SMARCA4 | 5 | | BMPR1A | 3 |
| WRN | 4 | | TP53 | 3 |
| ACD | 4 | | SMAD4 | 3 |
| BRCA2 | 4 | | (and 25 more, mostly 1-3 per gene) | |

**Why they're flagged**: SpliceAI `-M 1` (masked) is designed to **zero out
scores at canonical splice sites** — the whole point of masking is to find
disruption *outside* the canonical positions. So variants at the canonical
sites themselves automatically pass the `SpliceAI <= 0.1` filter regardless
of their true effect, and the W2 rule fires Benign. These calls are filter
artifacts, not evidence of safety.

**History**: this count moved twice.
- Pre-VV resolver: **255** rows. Most were c. on dbNSFP's alternative
  transcripts where positions falsely looked splice-like.
- Post-VV resolver, pre missing-SpliceAI-as-zero fix: **13** rows
  (the genuine canonical-splice cases on the BED's NM_).
- Post missing-SpliceAI-as-zero fix (current): **265** rows. With missing
  SpliceAI treated as 0, every canonical-splice position automatically
  passes the SpliceAI ≤ 0.1 gate — masked SpliceAI is near-zero by design
  at canonical sites.

This expansion is a known trade-off of the missing-SpliceAI rule that
recovers ~43 legitimate Benign calls at alt-transcript regions (APC E01,
RAD51D alt-E3, etc.) where SpliceAI has no gene model. Practical recipe
for SeqNext upload: exclude this review file (`_review_canonical_splice.tsv`)
from the upload, hand-review separately.

Options: filter from SeqNext upload, or re-run classifier with an
`--exclude-canonical-splice` flag (not yet implemented).

### `_review_intergenic.tsv` — 70 rows

Rows where VariantValidator returned `vv_status=flagged:intergenic` —
i.e. the variant is non-coding on the BED's specified RefSeq NM_.

| Gene | Count |
|---|---|
| FANCD2 (NM_033084.5) | 62 |
| SMARCA4 | 8 |

The classification probably still holds biologically (the SpliceAI / FAF /
REVEL filters don't depend on transcript choice), but the c. notation in
these rows is from a different transcript than the BED target and
shouldn't be uploaded to SeqNext as-is. Drop or re-resolve against an
alternative isoform.

### Other one-off flags

- **`TSC1 c.1334-2A>G`** in W1 (FAF >5%). Splice-acceptor position;
  population frequency overrides mechanistic prediction at >5% but worth
  confirming the FAF isn't artifact (e.g. mismapping at this position).

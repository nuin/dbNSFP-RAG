# cp_new pipeline

End-to-end annotation + classification pipeline for the 165-gene hereditary
cancer panel `cp_new`. Pulls dbNSFP scores + gnomAD v4.1 FAF + (locally-run)
SpliceAI masked scores, applies three lab-specific benign/LB rules, and emits
per-gene TSVs ready for SeqNext import.

## Scope

- **Panel**: `cp_new` — 165 genes (98 overlap with `NGSgenes`, 67 are new).
  Defined in `src/panels.py`.
- **BED**: `/Users/nuin/Projects/ahs/BED/CP_new/C+_ALL_IDPE_OCT2025.bed`
  (166 unique genes, 2380 regions). The BED commits one RefSeq transcript per
  region. 164 of 166 genes use a single transcript; **APC** and **RAD51D**
  carry two each (alternative first exon / alternative E3 respectively).
- **Output**: per-gene TSV under `data/exports/cp_new/seqnext/{GENE}_seqnext.tsv`
  + combined `cp_new_seqnext.tsv`. Columns: `gene, transcript, hgvs_c,
  classification`. `transcript` is the BED's RefSeq NM_ for the region the
  variant lies in.

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

## Snapshot — first full production run (2026-05-12)

- **dbNSFP scan**: 100,034 variant rows across 165 genes, 8,253 with gnomAD
  FAF (8%). 99,436 unique variant_ids emitted to cp_new.vcf.
- **BED-filter for SpliceAI**: 53,984 BED-covered unique variants.
- **SpliceAI run**: 53,574 of 53,984 (99.2%) got DS_MAX scores. Boundary
  variants on chunk splits accounted for the small miss.
- **Classifier output**: **4,567 SeqNext rows** across 100 genes.
  - W1 (Benign FAF >5%): **20**
  - W2 (Benign synonymous + SpliceAI≤0.1 + conservation low): **4,460**
  - W3 (Likely_benign FAF >0.1% + REVEL<0.29 + SpliceAI≤0.1): **87**
  - Still "pending": 58 (chunk-boundary missed-SpliceAI variants)
- **VV resolver**: 4,497 OK, 70 flagged:intergenic. Per-gene SeqNext files
  now carry the correct MANE/RefSeq c. notation (e.g. BRCA1 17:41197802 C>G
  → `NM_007294.4:c.5485G>C`, not the `c.2099G>C` dbNSFP would have given).
- **SQLite DB**: 272,875 → 303,628 rows. `cp_new` panel = 99,433 variants.

## Known issues & manual-review items

1. **255 W2 calls at canonical splice positions** (`-1`/`-2`/`+1`/`+2`).
   SpliceAI masked (`-M 1`) is *designed* to zero out scores at canonical
   splice sites — the whole point of masking is to find disruption *outside*
   the canonical positions. So variants at the canonical sites themselves
   automatically pass the `SpliceAI <= 0.1` filter regardless of their true
   effect, and the W2 rule fires Benign. Examples: `AIP c.469-2A>C`,
   `ANKRD26 c.639-1G>T`, `ACD c.986+1G>T`. These need manual review or to be
   filtered out of the classifier (add `--exclude-canonical-splice` flag).
2. **70 intergenic flags** in the SeqNext output. VariantValidator says
   these positions aren't coding on the BED's NM_ — mostly FANCD2
   `NM_033084.5`. Drop from SeqNext upload or hand-review against an
   alternative isoform.
3. **`TSC1 c.1334-2A>G`** appears in W1 (FAF >5%). Splice-acceptor position;
   population frequency overrides mechanistic prediction at >5% but worth
   confirming the FAF isn't artifact (e.g. mismapping at this position).

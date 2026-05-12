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
              spliceai -M 1 (external) ─► cp_new.spliceai.vcf
                       │                      │
                       ▼                      ▼
          scripts/annotate_spliceai.py joins back into per-gene TSVs
                       │
                       ▼
          scripts/classify_cp_new.py  ─►  seqnext/{GENE}_seqnext.tsv
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
and pairs 1:1 with an Ensembl ENST. CDS coordinates match, so dbNSFP's
`HGVSc_snpEff` (annotated against Ensembl) is the same as HGVS c. on the
RefSeq NM_ — we just relabel.

Two exceptions:

- **APC** — 15 of 16 BED regions on `NM_000038.6` (canonical). 1 region (E01
  at chr5:112,043,145–112,043,585) on `NM_001127511.3`, an isoform with an
  alternative first exon. Variants in that region get HGVS on the alternative.
- **RAD51D** — 10 regions on `NM_002878.3`. 1 extra E03 region (chr17:33,443,872–
  33,444,071, 1.6 kb upstream of canonical E03) on `NM_001142571.2`, capturing
  the alternative exon 3.

For these, `HGVSc_snpEff` from dbNSFP may not match what HGVS c. would be on
the BED-listed transcript. Current `classify_cp_new.py` emits dbNSFP's value
as-is; for production SeqNext upload you'll want to validate the dual-tx
regions against VariantValidator (or VEP --refseq) before signing off. Small
volume — the alternative regions cover ~10-50 variants combined.

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

# 5. load into SQLite (~15 sec)
uv run python scripts/load_cp_new_sqlite.py
```

## Snapshot at first full run

Run from steps 1-5 (step 2 still in progress as of writing):

- dbNSFP scan: 100,034 variants across 165 genes, 8,253 with gnomAD FAF (8%)
- After SQLite load: DB grew from 272,875 → 303,628 rows. `cp_new` panel
  registered with 99,433 variants.
- Workflow 1 (FAF >5%): **20 BENIGN** calls across 9 genes (SLX4 dominant
  with 11; the rest in EPCAM, FANCM, GATA2, GCM2, RPS20, SDHA, TRPV6, TSC1).
- Pending SpliceAI: 10,620 W2 candidates + 94 W3 candidates.

Notable flag: `TSC1 c.1334-2A>G` at FAF >5% is a canonical splice acceptor
position. Rule fires Benign on population freq, but normally ACMG BA1 doesn't
override a splice consequence. Manual review before SeqNext upload.

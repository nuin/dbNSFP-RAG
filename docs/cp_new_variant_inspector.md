# cp_new variant inspector

Alamut-style single-variant lookup tool. Pulls together everything we
have locally (cp_new pipeline cache, dbNSFP per-gene TSVs, ClinVar bulk,
VariantValidator cache, reference FASTAs) and augments it with free
public APIs (myvariant.info, VariantValidator REST, Ensembl REST). Two
front-ends share the same data layer:

- **CLI**: `scripts/inspect_variant.py`
- **Web UI**: `scripts/inspect_app.py` (Streamlit)

## CLI usage

```bash
# Single variant (gene + HGVS c.)
uv run python scripts/inspect_variant.py lookup --gene GOT2 --hgvs "c.816C>T"

# Single variant (genomic coords, hg19)
uv run python scripts/inspect_variant.py lookup --gene GOT2 \
    --coord 16:58750604:G:A --build hg19

# Single variant (dbSNP rsID)
uv run python scripts/inspect_variant.py lookup --rsid rs1058192

# Skip external APIs (instant, local cache only)
... --no-external

# JSON output (for scripting)
... --json

# Dump all classified variants for a gene
uv run python scripts/inspect_variant.py gene GOT2
uv run python scripts/inspect_variant.py gene GOT2 --format tsv

# Bulk -- process a TSV list, write per-variant Markdown files
uv run python scripts/inspect_variant.py bulk variants.tsv --out reports/
uv run python scripts/inspect_variant.py bulk variants.tsv --out reports/ --external  # include API calls
```

Bulk TSV input format: at minimum `gene` + `hgvs_c` columns; if
`chr_grch37/pos_grch37/ref/alt` are also present they're used to enrich
the lookup.

## Web UI

```bash
uv pip install streamlit  # one-time
uv run streamlit run scripts/inspect_app.py
# -> opens http://localhost:8501 in browser
```

Layout:
- **Sidebar**: variant input (gene+HGVS / coord / rsID), gene browser.
- **Main**: collapsible panels — Identity, Frequencies, Predictors, Splice
  & Conservation, ClinVar, our pipeline classification, Cross-references.
- **Gene browser**: type a gene symbol, see all classified variants for
  it in a sortable table, download as TSV.

Runs locally on the analyst's workstation. No auth, no network egress
beyond the configured external APIs.

## Data sources

### Local cache (instant lookup)

| Path | Contents |
|---|---|
| `data/exports/cp_new/{GENE}.tsv` | dbNSFP + gnomAD v4.1 per-gene annotations |
| `data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv` | merged classifications (pipeline + ClinVar + catalog) |
| `data/exports/cp_new/.vv_cache.json` | VariantValidator HGVS resolutions |
| `~/data/clinvar/variant_summary.txt.gz` | ClinVar bulk (refreshed weekly) |
| `~/data/hg19/hg19.fa`, `~/data/hg38/hg38.fa` | Reference FASTAs for codon context |

### Public APIs (cached on disk; rate-limited)

| Source | What it adds |
|---|---|
| `rest.variantvalidator.org` | HGVS resolution on any specified transcript |
| `myvariant.info` | Aggregator: dbSNP, ClinVar, COSMIC, EVS, ExAC, gnomAD, CADD |
| `grch37.rest.ensembl.org`, `rest.ensembl.org` | Gene/transcript metadata in both builds |
| `eutils.ncbi.nlm.nih.gov` | ClinVar variation summary + submitter details |

HTTP responses are cached at `data/exports/cp_new/.inspector_cache/` so
repeat lookups are instant. Cache files are JSON, gitignored.

## Example output

```
$ uv run python scripts/inspect_variant.py lookup --gene GOT2 --hgvs "c.816C>T"

# Variant report — GOT2 c.816C>T

## Identity
- gene: GOT2
- Ensembl tx (dbNSFP): (n/a -- synonymous, not in dbNSFP)
- HGVSc (dbNSFP/snpEff): (n/a)
- hg19: chr16:58750604 G>A
- hg38: chr16:58716700 G>A

## Frequencies (gnomAD)
- gnomAD v4.1 joint AF: 0.853
- v4.1 grpmax FAF95: 0.83
- ...

## Our pipeline classification
- Benign ClinVar(2-star, criteria provided)
  - source: clinvar    vv_status: clinvar_assertion
  - transcript: NM_002080.4
  - c.: c.816C>T

## Cross-references
- dbSNP rsID: rs1058192
- myvariant.cosmic: present (run --raw-mv to see full JSON)
- myvariant.evs: present
- ...
```

## What it shows vs Alamut

| Capability | Alamut | This tool |
|---|---|---|
| HGVS on all transcripts | ✅ | ✅ via VariantValidator REST |
| gnomAD frequencies | ✅ | ✅ (local cache + myvariant.info) |
| ClinVar assertions | ✅ | ✅ (bulk + esummary) |
| Conservation (PhastCons/PhyloP/GERP) | ✅ | ✅ from dbNSFP |
| In-silico predictors (REVEL/CADD/AlphaMissense/SIFT/PolyPhen) | ✅ | ✅ from dbNSFP |
| SpliceAI raw + masked | ✅ | ✅ (we cache masked from local run) |
| Alternative splice predictors (MaxEnt, NNSplice, SSF, ESEFinder) | ✅ | ❌ — Alamut's proprietary bundle |
| LOVD entries | ✅ | partial (myvariant.info passthrough only) |
| HGMD entries | ✅ (licensed) | ❌ — requires HGMD license |
| OMIM links | ✅ | ✅ via gene lookup |
| Reference context (codon, surrounding sequence) | ✅ | ✅ from local FASTA |
| Interactive splice diagram | ✅ | ❌ — needs custom render |
| Browser GUI | ✅ | ✅ via Streamlit (`inspect_app.py`) |
| **Per-variant cost** | $$ per query / seat license | **free** |

## Re-build local cache

The tool reads whatever's currently in `data/exports/cp_new/`. To refresh:
- `scripts/cp_new.py` for dbNSFP + gnomAD extraction
- `scripts/annotate_spliceai.py` after running SpliceAI
- `scripts/classify_cp_new.py` to apply W1/W2/W3
- `scripts/resolve_hgvs_vv.py` for VariantValidator cache
- `scripts/pull_clinvar_benign.py` for ClinVar bulk
- `scripts/merge_all_sources.py` to rebuild `cp_new_seqnext_FINAL.tsv`

## Limitations

- The local cache covers the **cp_new panel only**. For genes outside
  cp_new the inspector falls back to API-only mode (slower, less rich).
- ClinVar bulk is refreshed manually (`curl` per script docstring). For
  the latest assertions, re-download `variant_summary.txt.gz`.
- myvariant.info passes through whatever the upstream sources provided
  at its last refresh (~monthly); not strictly real-time.
- No HGMD coverage (license required).
- Alternative splice predictors (MaxEnt etc.) not included; SpliceAI is
  the splice-impact signal we use.

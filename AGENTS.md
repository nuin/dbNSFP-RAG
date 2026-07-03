# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project Overview

This is an ACMG Variant Classification API that uses a fine-tuned LLM to classify genetic variants according to ACMG/AMP guidelines. It processes dbNSFP annotation data, builds a FAISS vector database, and exposes a REST API for variant classification. **FOR RESEARCH USE ONLY** - not validated for clinical/diagnostic use.

## Common Commands

### Development Setup
```bash
uv venv && uv pip install -e .
uv pip install fastapi uvicorn mlx-lm  # for API
uv pip install -e ".[dev]"             # for testing
```

### Start the API Server
```bash
uvicorn api.server:app --host 0.0.0.0 --port 8000
uvicorn api.server:app --reload  # development mode

# Production (RHEL 8): PM2 manages the api on port 8029 + optional ollama
pm2 start ecosystem.config.js --only acmg-api
```

### Build SQLite Database (preferred, production backend)
```bash
# Builds one .db file covering ALL panels (no embeddings, no ML deps)
uv run python -m src.main build-sqlite --build grch37

# Resume / fresh / chr filter all supported
uv run python -m src.main build-sqlite --fresh --chr 17 --limit 50000
```
Output: `data/sqlite/grch37-all-panels.db` (~400 MB). Used by API in production
(see PM2 `ACMG_DB_PATH`).

### Build FAISS Vector Database (legacy)
```bash
# Build for specific gene panel with GRCh37 coordinates
uv run python -m src.main build-panel --panel NGSgenes --build grch37

# Build full index (all chromosomes)
uv run python -m src.main build

# Resume interrupted build (uses checkpoint)
uv run python -m src.main build-panel --panel NGSgenes

# Start fresh, ignore checkpoint
uv run python -m src.main build-panel --panel NGSgenes --fresh

# Limit variants for testing
uv run python -m src.main build --limit 50000
```

### Training Pipeline
```bash
# Generate ACMG training data from vector database
uv run python training/acmg_training.py --db data/vectordb/grch37-ngsgenes --output data/training/acmg_training

# Fine-tune on Apple Silicon (MLX)
uv run python training/finetune_acmg.py --method mlx --iters 500

# Fine-tune on NVIDIA GPU
uv run python training/finetune_acmg.py --method transformers

# Test existing model
uv run python training/finetune_acmg.py --test-only
```

### Utility Commands
```bash
uv run python -m src.main stats           # Show index statistics
uv run python -m src.main list            # List chromosomes in source
uv run python -m src.main list-panels     # List available gene panels
uv run python -m src.main clear           # Clear index and checkpoint
```

### Testing
```bash
pytest                                     # Run all tests
pytest validation/test_suite/test_api.py  # Run single test file
pytest -k "test_validate_chromosome"      # Run tests matching name
pytest -m "not slow"                       # Skip slow tests
pytest -m "not requires_model"             # Skip tests requiring ML model
```

### Linting
```bash
ruff check .
ruff format .
```

### Docker Deployment
```bash
# Build bundled image (includes model + database)
docker build -f Dockerfile.bundled -t acmg-api:bundled .

# Run container
docker run -d -p 8000:8000 --name acmg-api acmg-api:bundled

# Check logs
docker logs acmg-api --tail 50
```

## Architecture

### Data Pipeline Flow
1. **Ingestion** (`src/ingest.py`): Reads dbNSFP chromosome files (gzipped TSV) from zip or directory. Streams data in chunks to handle large files. Filters to keep only ~50 key columns from 600+.

2. **Gene Data** (`src/gene_data.py`): Loads gene-level constraint scores (pLI, LOEUF, mis_z) from dbNSFP gene file. Cached in memory for fast access during variant processing.

3. **Chunking** (`src/chunker.py`): Converts variant rows to structured text for embedding. Formats scores with interpretations (e.g., CADD phred thresholds). Derives consequence type from HGVSp/amino acid changes. Integrates gene constraint scores. Supports GRCh37/GRCh38 coordinate systems.

4. **Storage backends — pick one:**
   - **`src/variantdb.py` (SQLite, production default)**: Indexed coordinate/gene/ClinVar lookups, no ML deps, ~400 MB for all-panels. Schema in `VariantDatabase.SCHEMA`. The API auto-selects this when `ACMG_DB_PATH` points to a `.db` file.
   - **`src/vectorstore.py` (FAISS, legacy)**: sentence-transformer embeddings (all-MiniLM-L6-v2) for semantic search. Still buildable via `build`/`build-panel`. API falls back to it when `ACMG_DB_PATH` is a directory.

   Both backends return the same `{"id", "document", "metadata"}` shape, so downstream ACMG/API code is backend-agnostic.

5. **ACMG Scoring** (`src/acmg_scoring.py`): Rule-based evaluation of 28 ACMG/AMP criteria. 17 criteria can be automated from dbNSFP data. 11 criteria require external clinical/family data. Uses ClinVar API for PS1/PM5 evaluation when configured.

6. **Validation** (`src/validation.py`): Input validation for variants, HGVS notation, and ACMG codes. Returns `ValidationResult` dataclass with normalized values or detailed errors.

7. **Training** (`training/acmg_training.py`): Generates instruction-tuning data from ClinVar annotations. Maps ClinVar significance to ACMG 5-tier classification. Outputs Alpaca JSON format for mlx-lm.

8. **Fine-tuning** (`training/finetune_acmg.py`): LoRA fine-tuning using mlx-lm (Apple Silicon) or transformers (NVIDIA). Uses Llama-3.2-3B-Instruct-4bit as base model.

9. **API** (`api/server.py`): FastAPI server. Endpoints:
   - `GET|POST /classify` — accepts `chr/pos/ref/alt` **or** `hgvs=` (e.g. `17:g.41197801T>A`), plus `genome_build` (GRCh37/GRCh38, default GRCh37). Returns ACMG classification from rule-based scoring on stored dbNSFP data (LLM fallback was removed in commit 8707ed8).
   - `GET /variant/{id}`, `GET /gene/{symbol}` — direct lookups
   - `GET /chromosomes`, `/chromosome/{c}/genes`, `/panels`, `/panel/{name}/genes` — browse navigation (SQLite backend only; returns 501 on FAISS)
   - `GET /lookup` — dbNSFP raw record lookup

   `get_database(genome_build)` lazy-loads and caches per-build databases in `_databases` dict.

### Key Data Structures

**Variant ID format**: `{chr}_{pos}_{ref}_{alt}` (e.g., `17_41197801_T_A`)

**Metadata stored per variant**:
- chr, pos, ref, alt, gene, transcript
- cadd_phred, revel_score, gnomad_af, gnomad_hom
- sift_pred, polyphen_pred, alphamissense_pred
- clinvar_sig, clinvar_id, clinvar_review
- consequence (derived from HGVSp: missense_variant, stop_gained, etc.)
- gnomad_pli, loeuf, gnomad_mis_oe (gene constraint scores)
- interpro_domain (functional domain annotation)

**ValidationResult** (from `src/validation.py`):
- `valid`: bool - whether validation passed
- `error`: optional error message
- `warning`: optional warning message
- `normalized`: normalized value
- `variant_id`: constructed variant ID (for `validate_variant`)

**ACMGCriterion** (from `src/acmg_scoring.py`):
- `code`: ACMG code (e.g., "PVS1", "PM2")
- `status`: CriterionStatus (MET, NOT_MET, NOT_EVALUATED)
- `evidence`: Human-readable explanation
- `strength`: EvidenceStrength (very_strong, strong, moderate, supporting)

### ACMG Criteria Implementation

**Automated from dbNSFP (17 criteria)**:
- PVS1: LOF variants with LOEUF/pLI constraint
- PM1: Functional domain (InterPro)
- PM2: Rare variant (gnomAD AF)
- PM4: In-frame protein length change
- PS1/PM5: ClinVar same/different AA change (requires ClinVar API)
- PP2: Missense in constrained gene (mis_z)
- PP3/BP4: Computational predictor consensus
- PP5/BP6: ClinVar assertions
- BA1/BS1: Common variants
- BS2: Homozygotes in gnomAD
- BP1: Missense in truncating-variant gene
- BP7: Synonymous variant

**Require external data (11 criteria)**:
- PS2, PM6: De novo status (trio data)
- PS3, BS3: Functional studies (literature)
- PS4: Case-control prevalence
- PM3, BP2: Phasing (cis/trans)
- PP1, BS4: Segregation
- PP4: Patient phenotype
- BP5: Alternate diagnosis

### Gene Panels
Defined in `src/panels.py`. Main panel is `NGSgenes` (314 genes covering cardiac + cancer). Other panels: `hereditary_cancer`, `cardiac`, `neurological`.

### Environment Variables
- `ACMG_MODEL_PATH`: Path to fine-tuned model (default: `models/acmg-classifier/model`)
- `ACMG_DB_PATH`: Path to GRCh37 variant database. **Suffix decides backend**: `.db` → SQLite (`VariantDatabase`), directory → FAISS (`VariantVectorStore`). Default: `data/sqlite/grch37-all-panels.db`.
- `ACMG_DB_PATH_GRCH38`: Same, for GRCh38 (default: `data/vectordb/grch38-ngsgenes`)
- `DBNSFP_DEVICE`: Force device for embeddings (`cpu`, `cuda`, `mps`)
- `DBNSFP_DIR`: Path to dbNSFP data directory
- `OLLAMA_HOST`: Ollama server URL for fallback inference

### Coordinate System
Default is GRCh37 (hg19) for clinical compatibility. The `use_grch37` flag in chunking/metadata functions controls which coordinates are used. Variants without GRCh37 liftover are skipped when building GRCh37 databases.

### Test Infrastructure
Tests are in `validation/test_suite/`. Pytest markers: `slow`, `requires_model`, `requires_database`, `integration`. Test fixtures in `conftest.py`.

## Important Implementation Details

### Backend Selection (SQLite vs FAISS)
`api/server.py:get_database()` inspects the path suffix:
- `.db` and file exists → loads `VariantDatabase`
- otherwise → loads `VariantVectorStore` (FAISS)

If `ACMG_DB_PATH` is set to a `.db` path that doesn't exist, the API silently falls back to the default FAISS directory. This is intentional for backwards compatibility but can mask a missing/misplaced SQLite file — check startup logs for "Loaded GRCh37 SQLite/FAISS database with N variants".

### Database Rebuilding
When modifying metadata fields in `src/chunker.py` or `src/config.py`:
1. The database must be rebuilt for changes to take effect
2. Use `--fresh` flag to ignore checkpoints
3. For SQLite: rebuild via `build-sqlite`; the schema in `VariantDatabase.SCHEMA` must be updated if new metadata fields are added (otherwise they'll only appear in `full_annotation` JSON, not be queryable).
4. For FAISS: stored in `data/vectordb/grch37-ngsgenes/`
5. Docker images bundle the database and must be rebuilt

### ACMG Scoring Logic
The `evaluate_all_criteria()` function in `src/acmg_scoring.py`:
- Returns list of `ACMGCriterion` objects, one per ACMG code
- Each criterion has status: MET, NOT_MET, or NOT_EVALUATED
- `combine_criteria()` applies official ACMG combination rules
- PVS1 requires consequence type + gene constraint (pLI/LOEUF)
- PP2/BP1 use gnomAD mis_z scores for missense constraint

### Consequence Derivation
In `src/chunker.py`, `derive_consequence()` function:
- Parses HGVSp notation to determine variant consequence
- Maps patterns like `p.X123*` → `stop_gained`
- Maps `p.A123B` → `missense_variant`
- Returns "unknown" if HGVSp cannot be parsed
- Critical for PVS1, PM4, BP1, BP7 criteria evaluation

### Gene Constraint Integration
`src/gene_data.py` loads gene file once at startup:
- Reads `dbNSFP5.3_gene.gz` from DBNSFP_DIR
- Caches pLI, LOEUF, mis_z scores in memory
- `get_gene_scores()` returns dict or None
- Used by chunker to add metadata, by ACMG scoring for PVS1/PP2/BP1

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an ACMG Variant Classification API that uses a fine-tuned LLM to classify genetic variants according to ACMG/AMP guidelines. It processes dbNSFP annotation data, builds a FAISS vector database, and exposes a REST API for variant classification.

## Common Commands

### Development Setup
```bash
uv venv && uv pip install -e .
uv pip install fastapi uvicorn mlx-lm  # for API
```

### Start the API Server
```bash
uvicorn api.server:app --host 0.0.0.0 --port 8000
uvicorn api.server:app --reload  # development mode
```

### Build Vector Database
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

### Linting
```bash
ruff check .
ruff format .
```

## Architecture

### Data Pipeline Flow
1. **Ingestion** (`src/ingest.py`): Reads dbNSFP chromosome files (gzipped TSV) from zip or directory. Streams data in chunks to handle large files. Filters to keep only ~50 key columns from 600+.

2. **Chunking** (`src/chunker.py`): Converts variant rows to structured text for embedding. Formats scores with interpretations (e.g., CADD phred thresholds). Supports GRCh37/GRCh38 coordinate systems.

3. **Vector Store** (`src/vectorstore.py`): FAISS-backed store with sentence-transformer embeddings (all-MiniLM-L6-v2). Stores variant ID, text document, and metadata. Supports semantic search, gene lookup, and region queries.

4. **Training** (`training/acmg_training.py`): Generates instruction-tuning data from ClinVar annotations. Maps ClinVar significance to ACMG 5-tier classification. Outputs Alpaca JSON format for mlx-lm.

5. **Fine-tuning** (`training/finetune_acmg.py`): LoRA fine-tuning using mlx-lm (Apple Silicon) or transformers (NVIDIA). Uses Llama-3.2-3B-Instruct-4bit as base model.

6. **API** (`api/server.py`): FastAPI server with GET/POST `/classify` endpoints. Lazy-loads model (MLX > transformers > Ollama fallback). Returns ACMG classification, criteria, and confidence.

### Key Data Structures

**Variant ID format**: `{chr}_{pos}_{ref}_{alt}` (e.g., `17_41197801_T_A`)

**Metadata stored per variant**:
- chr, pos, ref, alt, gene, transcript
- cadd_phred, revel_score, gnomad_af
- sift_pred, polyphen_pred, alphamissense_pred
- clinvar_sig

### Gene Panels
Defined in `src/panels.py`. Main panel is `NGSgenes` (314 genes covering cardiac + cancer). Other panels: `hereditary_cancer`, `cardiac`, `neurological`.

### Environment Variables
- `ACMG_MODEL_PATH`: Path to fine-tuned model (default: `models/acmg-classifier/model`)
- `ACMG_DB_PATH`: Path to variant database (default: `data/vectordb/grch37-ngsgenes`)
- `DBNSFP_DEVICE`: Force device for embeddings (`cpu`, `cuda`, `mps`)
- `DBNSFP_DIR`: Path to dbNSFP data directory
- `OLLAMA_HOST`: Ollama server URL for fallback inference

### Coordinate System
Default is GRCh37 (hg19) for clinical compatibility. The `use_grch37` flag in chunking/metadata functions controls which coordinates are used. Variants without GRCh37 liftover are skipped when building GRCh37 databases.

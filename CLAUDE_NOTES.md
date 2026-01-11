# Notes for Next Claude Instance

Hey there! Here's the context on this project:

## What This Is
A RAG (Retrieval Augmented Generation) pipeline for dbNSFP variant annotation at Alberta Precision Labs. It indexes variants into a FAISS vector store and uses Ollama for LLM-based clinical interpretation.

## Current State (2026-01-10)

### Completed
- GRCh37 NGSgenes database built (268,954 variants)
- Training pipeline fully operational
- All three training types tested and working
- Full documentation created (README.md)
- CPU-only support added via environment variables

### Trained Models
| Model | Location | Performance |
|-------|----------|-------------|
| XGBoost Classifier | `models/classifier.xgb` | 97% accuracy, ROC-AUC 0.9965 |
| Fine-tuned Embeddings | `models/variant-embeddings/` | Separation: 0.1349 |
| LLM Training Data | `data/exports/NGSgenes_llm_training.alpaca.json` | 232,838 samples |

### Database
- **Location:** `data/vectordb/grch37-ngsgenes/`
- **Size:** ~695 MB (413 MB faiss.index + 282 MB metadata.pkl)
- **Variants:** 268,954
- **Genes:** 314 (NGSgenes panel)

## Source Data
Extracted dbNSFP files at: `~/Downloads/dbNSFP5.3.1a/`

## Key Commands

```bash
# Build database
uv run python -m src.main build-panel --panel NGSgenes --build grch37

# Export training data
uv run python -m src.export --panel NGSgenes --type all --db data/vectordb/grch37-ngsgenes

# Train classifier
uv run python training/train_classifier.py --input data/exports/NGSgenes_classifier_features.jsonl

# Train embeddings
uv run python training/train_embeddings.py --input data/exports/NGSgenes_embedding_pairs.jsonl

# Prepare LLM data
uv run python training/train_llm.py --input data/exports/NGSgenes_llm_training.jsonl --prepare-only

# RAG interpretation
uv run python -c "
from src.rag import VariantRAG
from src.config import VECTORDB_GRCH37_NGSGENES
rag = VariantRAG(db_path=VECTORDB_GRCH37_NGSGENES, model='llama3.2:3b')
print(rag.interpret('17_41197801_T_A'))
"
```

## CPU-Only Usage

Set environment variables for non-GPU machines:

```bash
export DBNSFP_DEVICE=cpu
export DBNSFP_DIR=/path/to/dbNSFP5.3.1a  # Optional
```

## Gene Panels (`src/panels.py`)
- `NGSgenes`: 314 genes (main clinical panel - cardiac + cancer)
- `hereditary_cancer`: 25 genes
- `cardiac`: 19 genes
- `neurological`: 14 genes

## Key Files
- `src/config.py` - Configuration (supports env vars: DBNSFP_DIR, DBNSFP_DEVICE)
- `src/vectorstore.py` - FAISS wrapper with CPU/GPU auto-detection
- `src/rag.py` - LLM integration with Ollama
- `src/export.py` - Training data export (supports --db parameter)
- `src/main.py` - CLI with build-panel command
- `training/` - Training scripts for classifier, embeddings, LLM
- `README.md` - Full documentation

## Things That Work
- Streaming from extracted directory
- Gene panel filtering during ingestion
- GRCh37 coordinate support
- Checkpoint/resume for long indexing jobs
- Semantic search by query or gene
- LLM interpretation via Ollama (llama3.2:3b tested)
- XGBoost pathogenicity classifier
- Fine-tuned embeddings for better similarity search
- CPU-only mode via DBNSFP_DEVICE=cpu

## Tested Variants
```
17_41197801_T_A  # BRCA1
13_32953652_G_A  # BRCA2
10_90701009_C_T  # ACTA2 (Pathogenic/Likely_pathogenic)
```

## GitHub
https://github.com/nuin/dbNSFP-RAG

— Claude Instance (2026-01-10)

# Notes for Next Claude Instance

## What This Is
ACMG Variant Classification API - a standalone fine-tuned LLM for classifying genetic variants according to ACMG/AMP guidelines. Trained on 314 clinical genes (NGSgenes panel) from dbNSFP.

## Current State (2026-01-10)

### Completed
- GRCh37 NGSgenes database: 268,954 variants
- ACMG training data: 9,743 labeled examples
- Fine-tuned Llama-3.2-3B model via LoRA (500 iterations)
- FastAPI server with /classify endpoint
- Full documentation

### Key Components

| Component | Location | Description |
|-----------|----------|-------------|
| ACMG Model | `models/acmg-classifier/model/` | Fine-tuned Llama-3.2-3B |
| API Server | `api/server.py` | FastAPI with /classify endpoint |
| Training Data | `data/training/acmg_training.json` | 9,743 ACMG examples |
| Database | `data/vectordb/grch37-ngsgenes/` | 268,954 variants |

### Model Training
```
Base: mlx-community/Llama-3.2-3B-Instruct-4bit
Method: LoRA (16 layers, 500 iterations)
Initial loss: 4.033
Final loss: 0.203
Training time: ~10 min on M1 Ultra
```

## Quick Start

```bash
# Start API
uvicorn api.server:app --host 0.0.0.0 --port 8000

# Classify variant
curl "http://localhost:8000/classify?chr=17&pos=41197801&ref=T&alt=A"
```

## API Endpoints

- `GET /classify?chr=X&pos=Y&ref=R&alt=A` - Classify variant
- `POST /classify` - Classify with JSON body
- `GET /gene/{symbol}` - Get variants for gene
- `GET /health` - API status

## Training Commands

```bash
# Generate ACMG training data
uv run python training/acmg_training.py --db data/vectordb/grch37-ngsgenes

# Fine-tune model (Apple Silicon)
uv run python training/finetune_acmg.py --method mlx --iters 500

# Fine-tune model (NVIDIA)
uv run python training/finetune_acmg.py --method transformers
```

## Key Files

- `api/server.py` - FastAPI server
- `training/acmg_training.py` - Generate ACMG training data
- `training/finetune_acmg.py` - Fine-tune LLM
- `src/vectorstore.py` - FAISS database (supports DBNSFP_DEVICE env var)
- `src/panels.py` - Gene panels (NGSgenes: 314 genes)

## Environment Variables

- `ACMG_MODEL_PATH` - Path to fine-tuned model
- `ACMG_DB_PATH` - Path to variant database
- `DBNSFP_DEVICE` - Force cpu/cuda/mps

## Non-GPU Deployment

```bash
export DBNSFP_DEVICE=cpu
export ACMG_MODEL_PATH=  # Empty = use Ollama
ollama pull llama3.2:3b
uvicorn api.server:app --port 8000
```

## GitHub
https://github.com/nuin/dbNSFP-RAG

— Claude Instance (2026-01-10)

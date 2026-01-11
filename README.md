# ACMG Variant Classification API

A standalone LLM-powered API for ACMG variant classification trained on 314 clinical genes (NGSgenes panel) using dbNSFP annotations.

## Features

- **Fine-tuned LLM** for ACMG/AMP variant classification
- **REST API** - GET/POST endpoints for variant classification
- **No source file dependency** - Model has learned patterns from 268,954 variants
- **GRCh37 coordinates** - Uses hg19 positions
- **Apple Silicon optimized** - MLX-based inference on M1/M2/M3/M4

## Quick Start

```bash
# Install
git clone https://github.com/nuin/dbNSFP-RAG.git
cd dbNSFP-RAG
uv venv && uv pip install -e .

# Start API
uvicorn api.server:app --host 0.0.0.0 --port 8000

# Classify a variant
curl "http://localhost:8000/classify?chr=17&pos=41197801&ref=T&alt=A"
```

## API Endpoints

### GET /classify

Classify a variant by genomic coordinates.

```bash
curl "http://localhost:8000/classify?chr=17&pos=41197801&ref=T&alt=A"
```

**Parameters:**
- `chr` - Chromosome (1-22, X, Y)
- `pos` - Position (1-based, GRCh37)
- `ref` - Reference allele
- `alt` - Alternate allele

**Response:**
```json
{
  "variant_id": "17_41197801_T_A",
  "gene": "BRCA1",
  "acmg_classification": "Likely_pathogenic",
  "criteria": [],
  "confidence": 0.75,
  "interpretation": "...",
  "scores": {
    "cadd_phred": 19.09,
    "sift_pred": "",
    "polyphen_pred": ""
  }
}
```

### POST /classify

Same as GET but with JSON body.

```bash
curl -X POST "http://localhost:8000/classify" \
  -H "Content-Type: application/json" \
  -d '{"chr": "17", "pos": 41197801, "ref": "T", "alt": "A"}'
```

### GET /gene/{symbol}

Get variants for a gene.

```bash
curl "http://localhost:8000/gene/BRCA1?limit=10"
```

### GET /health

Check API status.

```bash
curl "http://localhost:8000/health"
```

## Installation

### Requirements

- Python 3.11+
- 8GB+ RAM
- Apple Silicon (M1/M2/M3/M4) or NVIDIA GPU recommended

### Install Dependencies

```bash
# Using uv (recommended)
uv venv
uv pip install -e .
uv pip install fastapi uvicorn mlx-lm

# Or pip
pip install -e .
pip install fastapi uvicorn mlx-lm
```

### CPU-Only Mode

```bash
export DBNSFP_DEVICE=cpu
```

## Training the Model

### 1. Build Database (if not provided)

```bash
# Download dbNSFP to ~/Downloads/dbNSFP5.3.1a/

# Build GRCh37 database for NGSgenes panel
uv run python -m src.main build-panel --panel NGSgenes --build grch37
```

### 2. Generate Training Data

```bash
uv run python training/acmg_training.py \
  --db data/vectordb/grch37-ngsgenes \
  --output data/training/acmg_training
```

Output:
- `data/training/acmg_training.jsonl` - Training examples
- `data/training/acmg_training.json` - Alpaca format

### 3. Fine-tune Model

```bash
# Apple Silicon (MLX)
uv run python training/finetune_acmg.py --method mlx --iters 500

# NVIDIA GPU
uv run python training/finetune_acmg.py --method transformers
```

Output: `models/acmg-classifier/model/`

### 4. Start API

```bash
export ACMG_MODEL_PATH=models/acmg-classifier/model
uvicorn api.server:app --host 0.0.0.0 --port 8000
```

## Gene Panels

| Panel | Genes | Description |
|-------|-------|-------------|
| NGSgenes | 314 | Main clinical panel (cardiac + cancer) |
| hereditary_cancer | 25 | Cancer predisposition genes |
| cardiac | 19 | Cardiac disease genes |
| neurological | 14 | Neurological disease genes |

## ACMG Classifications

The model outputs one of 5 classifications:

| Classification | Description |
|----------------|-------------|
| Pathogenic | Strong evidence of disease causation |
| Likely_pathogenic | Moderate evidence of disease causation |
| Uncertain_significance | Insufficient evidence to classify |
| Likely_benign | Moderate evidence against pathogenicity |
| Benign | Strong evidence against pathogenicity |

## Training Data Distribution

From NGSgenes panel (9,743 labeled variants):

| Class | Count |
|-------|-------|
| Pathogenic | 1,011 |
| Likely_pathogenic | 997 |
| Uncertain_significance | 5,000 |
| Likely_benign | 2,484 |
| Benign | 251 |

## Model Performance

Training metrics (500 iterations on Llama-3.2-3B):
- Initial validation loss: 4.033
- Final validation loss: 0.203
- Training time: ~10 minutes on M1 Ultra

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `ACMG_MODEL_PATH` | Path to fine-tuned model | `models/acmg-classifier/model` |
| `ACMG_DB_PATH` | Path to variant database | `data/vectordb/grch37-ngsgenes` |
| `DBNSFP_DEVICE` | Device for embeddings | Auto-detect |

## Project Structure

```
dbNSFP-RAG/
├── api/
│   └── server.py          # FastAPI server
├── training/
│   ├── acmg_training.py   # Generate ACMG training data
│   └── finetune_acmg.py   # Fine-tune LLM
├── src/
│   ├── config.py          # Configuration
│   ├── vectorstore.py     # FAISS database
│   ├── panels.py          # Gene panel definitions
│   └── main.py            # CLI
├── models/
│   └── acmg-classifier/   # Fine-tuned model
└── data/
    ├── vectordb/          # Variant databases
    └── training/          # Training data
```

## Non-GPU Deployment

For servers without GPU:

```bash
# Set CPU mode
export DBNSFP_DEVICE=cpu

# Use Ollama instead of MLX
export ACMG_MODEL_PATH=  # Leave empty to use Ollama

# Install and start Ollama
ollama serve &
ollama pull llama3.2:3b

# Start API
uvicorn api.server:app --host 0.0.0.0 --port 8000
```

## License

MIT License

## Citation

If using dbNSFP data:
> Liu X, et al. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Med. 2020.

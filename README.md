# dbNSFP RAG Pipeline

A Retrieval-Augmented Generation (RAG) pipeline for clinical variant interpretation using dbNSFP annotations. Developed for Alberta Precision Labs.

## Features

- **Vector database** of variant annotations (FAISS + sentence-transformers)
- **Gene panel filtering** during ingestion for focused databases
- **GRCh37/GRCh38 coordinate support**
- **LLM-powered interpretation** via Ollama
- **Training pipeline** for pathogenicity classifier, embeddings, and LLM fine-tuning

## Quick Start

```bash
# Clone and install
git clone https://github.com/nuin/dbNSFP-RAG.git
cd dbNSFP-RAG
uv venv && uv pip install -e .

# Interpret a variant (requires pre-built database and Ollama)
uv run python -c "
from src.rag import VariantRAG
rag = VariantRAG()
print(rag.interpret('17_41197801_T_A'))
"
```

## Installation

### Requirements

- Python 3.11+
- 8GB+ RAM (16GB+ recommended for large databases)
- ~50GB disk space for full dbNSFP data

### Standard Installation (with GPU support)

```bash
# Using uv (recommended)
uv venv
uv pip install -e .

# Or using pip
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### CPU-Only Installation

For machines without GPU (CUDA/MPS), set the environment variable before running:

```bash
# Force CPU mode
export DBNSFP_DEVICE=cpu

# Or set in your shell profile (~/.bashrc, ~/.zshrc)
echo 'export DBNSFP_DEVICE=cpu' >> ~/.zshrc
```

**Performance notes for CPU-only:**
- Embedding generation: ~3-5x slower than GPU
- Database building: ~2-4 hours for NGSgenes panel (vs ~30 min on GPU)
- RAG queries: ~1-2 seconds per query (vs ~200ms on GPU)
- Training: Significantly slower, consider using smaller batch sizes

### Install Ollama (for LLM interpretation)

```bash
# macOS
brew install ollama

# Linux
curl -fsSL https://ollama.com/install.sh | sh

# Start Ollama and pull a model
ollama serve &
ollama pull llama3.2:3b
```

## Data Setup

### Download dbNSFP

1. Download dbNSFP from https://sites.google.com/site/jpaborern/dbnsfp
2. Extract to a folder (e.g., `~/Downloads/dbNSFP5.3.1a/`)

### Configure Data Path

Edit `src/config.py`:

```python
DBNSFP_DIR = Path.home() / "Downloads" / "dbNSFP5.3.1a"
```

Or set environment variable:

```bash
export DBNSFP_DIR=/path/to/dbNSFP5.3.1a
```

## Building the Database

### Build Panel-Specific Database (Recommended)

Build a focused database for specific gene panels:

```bash
# List available panels
uv run python -m src.main list-panels

# Build GRCh37 database for NGSgenes panel (314 genes)
uv run python -m src.main build-panel --panel NGSgenes --build grch37

# Build GRCh38 database
uv run python -m src.main build-panel --panel NGSgenes --build grch38

# Test with single chromosome first
uv run python -m src.main build-panel --panel NGSgenes --build grch37 --chr 17
```

### Available Gene Panels

| Panel | Genes | Description |
|-------|-------|-------------|
| NGSgenes | 314 | Main clinical panel (cardiac + cancer) |
| hereditary_cancer | 25 | Cancer predisposition genes |
| cardiac | 19 | Cardiac disease genes |
| neurological | 14 | Neurological disease genes |

### Database Locations

- Default (GRCh38): `data/vectordb/`
- GRCh37 NGSgenes: `data/vectordb/grch37-ngsgenes/`

## Usage

### RAG Interpretation

```python
from src.rag import VariantRAG
from src.config import VECTORDB_GRCH37_NGSGENES

# Initialize with GRCh37 database
rag = VariantRAG(db_path=VECTORDB_GRCH37_NGSGENES, model="llama3.2:3b")

# Interpret a variant
result = rag.interpret("17_41197801_T_A")
print(result)

# Also supports underscore format
result = rag.interpret("17_41197801_T_A")
```

### Vector Store Search

```python
from src.vectorstore import VariantVectorStore
from src.config import VECTORDB_GRCH37_NGSGENES

store = VariantVectorStore(db_path=VECTORDB_GRCH37_NGSGENES)

# Search by gene
results = store.search_by_gene("BRCA1", k=10)
for r in results:
    print(f"{r['id']} - {r['metadata'].get('clinvar_sig', 'N/A')}")

# Search by variant ID
result = store.get_variant("17_41197801_T_A")
```

### Command Line Interface

```bash
# Interpret variant
uv run python -m src.main interpret 17_41197801_T_A

# Search by gene
uv run python -m src.main search --gene BRCA1 --limit 10

# Build database
uv run python -m src.main build-panel --panel NGSgenes --build grch37
```

## Training Pipeline

### 1. Export Training Data

```bash
# Export all training data types for a panel
uv run python -m src.export --panel NGSgenes --type all --db data/vectordb/grch37-ngsgenes

# Export specific type
uv run python -m src.export --panel NGSgenes --type classifier
uv run python -m src.export --panel NGSgenes --type llm
uv run python -m src.export --panel NGSgenes --type embeddings
```

**Exported files:**
- `data/exports/NGSgenes_classifier_features.jsonl` - For XGBoost classifier
- `data/exports/NGSgenes_llm_training.jsonl` - For LLM fine-tuning
- `data/exports/NGSgenes_embedding_pairs.jsonl` - For embedding fine-tuning

### 2. Train Pathogenicity Classifier

XGBoost classifier using CADD, REVEL, and population frequency features.

```bash
uv run python training/train_classifier.py \
    --input data/exports/NGSgenes_classifier_features.jsonl \
    --output models/classifier.json
```

**Output:** `models/classifier.xgb`

**CPU-only notes:**
- Training is fast (~1-2 minutes) even on CPU
- No special configuration needed

### 3. Train Embedding Model

Fine-tune sentence-transformers for better variant similarity search.

```bash
uv run python training/train_embeddings.py \
    --input data/exports/NGSgenes_embedding_pairs.jsonl \
    --output models/variant-embeddings
```

**Output:** `models/variant-embeddings/`

**CPU-only notes:**
- Set smaller batch size: Edit `train_embeddings.py` line 57: `batch_size=8` (instead of 16)
- Training takes ~15-30 minutes on CPU (vs ~2 minutes on GPU)

### 4. Prepare LLM Training Data

Convert to Alpaca format for LLM fine-tuning.

```bash
uv run python training/train_llm.py \
    --input data/exports/NGSgenes_llm_training.jsonl \
    --prepare-only
```

**Output:** `data/exports/NGSgenes_llm_training.alpaca.json`

#### LLM Fine-tuning Options

**Option A: Apple Silicon (MLX)**
```bash
pip install mlx-lm

mlx_lm.lora \
    --model meta-llama/Llama-3.2-3B-Instruct \
    --data data/exports \
    --train \
    --batch-size 4 \
    --lora-layers 16 \
    --iters 1000
```

**Option B: NVIDIA GPU (Transformers + LoRA)**
```bash
pip install transformers datasets accelerate peft bitsandbytes

# See training/train_llm.py for full script
```

**Option C: CPU-Only (Not Recommended)**
LLM fine-tuning on CPU is extremely slow. Consider:
- Using pre-trained models with RAG (no fine-tuning)
- Using cloud GPU instances (Colab, Lambda Labs, etc.)
- Using smaller models (TinyLlama, Phi-2)

## Using Trained Models

### Load Fine-tuned Embeddings

```python
from src.vectorstore import VariantVectorStore

# Use custom embedding model
store = VariantVectorStore(
    db_path="data/vectordb/grch37-ngsgenes",
    model_name="models/variant-embeddings"
)
```

### Load Classifier

```python
import xgboost as xgb
import numpy as np

model = xgb.XGBClassifier()
model.load_model("models/classifier.xgb")

# Predict (features: cadd_phred, revel_score, gnomad_af)
features = np.array([[25.0, 0.8, 0.0001]])
prediction = model.predict(features)
probability = model.predict_proba(features)[:, 1]
```

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `DBNSFP_DIR` | Path to extracted dbNSFP data | `~/Downloads/dbNSFP5.3.1a` |
| `DBNSFP_DEVICE` | Force device (`cpu`, `cuda`, `mps`) | Auto-detect |
| `OLLAMA_URL` | Ollama server URL | `http://localhost:11434` |

### Config File (`src/config.py`)

Key settings:
```python
EMBEDDING_MODEL = "all-MiniLM-L6-v2"  # or custom model path
EMBEDDING_BATCH_SIZE = 1024           # Reduce for CPU (256-512)
CHUNK_SIZE = 10000                    # Rows per batch during parsing
DEFAULT_MODEL = "llama3.1:8b"         # Ollama model for RAG
```

## CPU-Only Quick Reference

```bash
# 1. Set CPU mode
export DBNSFP_DEVICE=cpu

# 2. Install with minimal dependencies
uv pip install -e .

# 3. Build database (slower but works)
uv run python -m src.main build-panel --panel hereditary_cancer --build grch37

# 4. Train classifier (fast on CPU)
uv run python training/train_classifier.py --input data/exports/hereditary_cancer_classifier_features.jsonl

# 5. Use RAG (works but slower queries)
uv run python -m src.main interpret 17_41197801_T_A
```

## Project Structure

```
dbNSFP-RAG/
├── src/
│   ├── config.py          # Configuration and paths
│   ├── ingest.py          # Data loading and filtering
│   ├── chunker.py         # Variant text/metadata generation
│   ├── vectorstore.py     # FAISS vector store wrapper
│   ├── rag.py             # LLM integration with Ollama
│   ├── panels.py          # Gene panel definitions
│   ├── export.py          # Training data export
│   └── main.py            # CLI entry point
├── training/
│   ├── train_classifier.py    # XGBoost pathogenicity classifier
│   ├── train_embeddings.py    # Sentence-transformer fine-tuning
│   └── train_llm.py           # LLM fine-tuning preparation
├── data/
│   ├── vectordb/              # FAISS databases
│   └── exports/               # Training data exports
├── models/
│   ├── classifier.xgb         # Trained XGBoost model
│   └── variant-embeddings/    # Fine-tuned embeddings
└── CLAUDE_NOTES.md            # Development notes
```

## Troubleshooting

### "CUDA out of memory"
```bash
export DBNSFP_DEVICE=cpu
# Or reduce batch size in config.py
```

### "Ollama connection refused"
```bash
ollama serve &
ollama list  # Verify running
```

### "Module not found" errors
```bash
# Run as module, not script
uv run python -m src.main ...
uv run python -m src.export ...
```

### Slow database building
- Use panel filtering instead of full database
- Process one chromosome at a time: `--chr 17`
- Resume from checkpoint (automatic)

## License

MIT License - See LICENSE file

## Citation

If using dbNSFP data, please cite:
> Liu X, et al. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Med. 2020.

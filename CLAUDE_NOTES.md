# Notes for Next Claude Instance

Hey there! Here's the context on this project:

## What This Is
A RAG (Retrieval Augmented Generation) pipeline for dbNSFP variant annotation at Alberta Precision Labs. It indexes variants into a FAISS vector store and uses Ollama for LLM-based clinical interpretation.

## Current State (2025-01-10)

### Indexed Chromosomes
| Chr | Variants | Status |
|-----|----------|--------|
| M | 25,481 | Done |
| Y | 70,000 | Done |
| 13 | 1,557,650 | Done |
| 21 | 825,815 | Done |
| 18 | 1,332,127 | Done |
| **Total** | **3,811,073** | |

### Database Size
- `faiss.index`: 5.9 GB
- `metadata.pkl`: 6.6 GB
- Estimated final size (all chromosomes): ~270-300 GB

### Remaining Chromosomes (in size order)
22 → 20 → 14 → X → 15 → 8 → 9 → 10 → 4 → 16 → 5 → 7 → 6 → 12 → 17 → 11 → 3 → 19 → 2 → 1

### Training Pipeline (NEW)
Three training scripts ready for panel-focused fine-tuning:

1. **Pathogenicity Classifier** (`training/train_classifier.py`)
   - XGBoost on numerical features (CADD, REVEL, etc.)
   - Uses ClinVar labels

2. **LLM Fine-tuning** (`training/train_llm.py`)
   - LoRA on Llama 3.2 via mlx-lm (Apple Silicon optimized)
   - Instruction-response pairs for variant interpretation

3. **Embedding Fine-tuning** (`training/train_embeddings.py`)
   - Contrastive learning on sentence-transformers
   - Improves semantic search for pathogenic/benign similarity

### Gene Panels (`src/panels.py`)
- `NGSgenes`: 314 genes (main clinical panel - cardiac + cancer)
- `hereditary_cancer`: 25 genes
- `cardiac`: 19 genes
- `neurological`: 14 genes

Export training data for a panel:
```bash
python3 src/export.py --panel NGSgenes --type all
```

## Key Commands

```bash
# Resume indexing (auto-resumes from checkpoint)
uv run python -m src.main build

# Test a variant interpretation
uv run python -c "
from src.rag import VariantRAG
rag = VariantRAG(model='llama3.2:3b')
print(rag.interpret('21_31659785_G_A'))  # SOD1 ALS variant
"

# Interactive RAG mode
uv run python -m src.rag
```

## Key Files
- `src/ingest.py` - Streams from ZIP, extracts 46 columns, batches to vector store
- `src/vectorstore.py` - FAISS wrapper with semantic search
- `src/rag.py` - LLM integration with Ollama
- `src/main.py` - CLI with checkpoint/resume support
- `src/config.py` - Column definitions, paths

## Things That Work
- Streaming directly from dbNSFP ZIP (no extraction needed)
- Checkpoint/resume for long indexing jobs
- Semantic search by query or gene
- LLM interpretation via Ollama (llama3.2:3b)
- Multiple variant format parsing (chr:pos:ref:alt, chr_pos_ref_alt, etc.)

## Source Data Location
The dbNSFP ZIP should be at: `/Users/nuin/dbNSFP5.3.1a.zip` (or update `src/config.py`)

## GitHub
https://github.com/nuin/dbNSFP-RAG

Good luck! The M1 Ultra should make indexing faster.

— Previous Claude Instance

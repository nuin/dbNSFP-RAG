#!/usr/bin/env python
"""
Fine-tune embedding model for better variant semantic search.

Uses contrastive learning to improve similarity between:
- Variants with same pathogenicity class
- Variants in same gene
- Variants with similar clinical significance

Usage:
    python training/train_embeddings.py --panel hereditary_cancer
    python training/train_embeddings.py --input data/exports/custom_embedding_pairs.jsonl
"""

import json
import argparse
from pathlib import Path


def train_sentence_transformer(data_path: Path, output_dir: Path):
    """
    Fine-tune sentence-transformers model using contrastive learning.

    Uses MultipleNegativesRankingLoss for efficient training.
    """
    try:
        from sentence_transformers import (
            SentenceTransformer,
            InputExample,
            losses,
        )
        from torch.utils.data import DataLoader
    except ImportError:
        print("Install: pip install sentence-transformers")
        return

    # Load base model (same as RAG uses)
    model = SentenceTransformer("all-MiniLM-L6-v2")

    # Load training pairs
    train_examples = []
    with open(data_path) as f:
        for line in f:
            data = json.loads(line)
            # MultipleNegativesRankingLoss uses anchor-positive pairs
            # Negatives are sampled from other pairs in batch
            train_examples.append(InputExample(
                texts=[data["anchor"], data["positive"]]
            ))

    print(f"Loaded {len(train_examples)} training pairs")

    # DataLoader
    train_dataloader = DataLoader(
        train_examples,
        shuffle=True,
        batch_size=16
    )

    # Loss function - uses in-batch negatives
    train_loss = losses.MultipleNegativesRankingLoss(model)

    # Train
    output_dir.mkdir(parents=True, exist_ok=True)

    model.fit(
        train_objectives=[(train_dataloader, train_loss)],
        epochs=3,
        warmup_steps=100,
        output_path=str(output_dir),
        show_progress_bar=True,
    )

    print(f"\nModel saved to {output_dir}")
    return model


def evaluate_embeddings(model_path: Path, test_data_path: Path):
    """Evaluate embedding quality on held-out test set."""
    try:
        from sentence_transformers import SentenceTransformer
        from sklearn.metrics.pairwise import cosine_similarity
        import numpy as np
    except ImportError:
        print("Install: pip install sentence-transformers scikit-learn")
        return

    model = SentenceTransformer(str(model_path))

    # Load test pairs
    pathogenic_texts = []
    benign_texts = []

    with open(test_data_path) as f:
        for line in f:
            data = json.loads(line)
            if data.get("label") == "pathogenic":
                pathogenic_texts.append(data["anchor"])
            else:
                benign_texts.append(data["anchor"])

    if not pathogenic_texts or not benign_texts:
        print("Need both pathogenic and benign examples for evaluation")
        return

    # Encode
    path_emb = model.encode(pathogenic_texts[:100])
    benign_emb = model.encode(benign_texts[:100])

    # Compute similarities
    path_path_sim = cosine_similarity(path_emb).mean()
    benign_benign_sim = cosine_similarity(benign_emb).mean()
    path_benign_sim = cosine_similarity(path_emb, benign_emb).mean()

    print("\n=== Embedding Quality ===")
    print(f"Pathogenic-Pathogenic similarity: {path_path_sim:.4f}")
    print(f"Benign-Benign similarity: {benign_benign_sim:.4f}")
    print(f"Pathogenic-Benign similarity: {path_benign_sim:.4f}")
    print(f"\nSeparation (higher is better): {(path_path_sim + benign_benign_sim) / 2 - path_benign_sim:.4f}")


def main():
    parser = argparse.ArgumentParser(description="Fine-tune embedding model")
    parser.add_argument("--panel", "-p", help="Panel name")
    parser.add_argument("--input", "-i", help="Pre-exported JSONL file")
    parser.add_argument("--output", "-o", default="models/variant-embeddings",
                       help="Output directory")
    parser.add_argument("--evaluate", action="store_true", help="Evaluate existing model")
    args = parser.parse_args()

    output_dir = Path(args.output)

    if args.evaluate:
        if args.input:
            evaluate_embeddings(output_dir, Path(args.input))
        else:
            print("Specify --input for evaluation data")
        return

    # Get training data
    if args.input:
        data_path = Path(args.input)
    elif args.panel:
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from src.export import export_for_embeddings
        data_path = export_for_embeddings(args.panel)
    else:
        print("Specify --panel or --input")
        return

    # Train
    train_sentence_transformer(data_path, output_dir)

    # Evaluate
    print("\nEvaluating trained model...")
    evaluate_embeddings(output_dir, data_path)


if __name__ == "__main__":
    main()

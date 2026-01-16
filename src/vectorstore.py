"""Vector store using FAISS for fast similarity search."""

import pickle
from pathlib import Path

import faiss
import numpy as np
from sentence_transformers import SentenceTransformer

from .config import EMBEDDING_BATCH_SIZE, EMBEDDING_MODEL, EMBEDDING_DIMENSION, VECTORDB_DIR


class VariantVectorStore:
    """FAISS-backed vector store for variant embeddings."""

    def __init__(
        self,
        db_path: Path | None = None,
        embedding_model: str | None = None,
    ):
        self.db_path = Path(db_path or VECTORDB_DIR)
        self.db_path.mkdir(parents=True, exist_ok=True)

        self.index_path = self.db_path / "faiss.index"
        self.meta_path = self.db_path / "metadata.pkl"

        # Initialize embedding model with device selection
        model_name = embedding_model or EMBEDDING_MODEL
        print(f"Loading embedding model: {model_name}")

        import os
        import torch

        # Check for forced device via environment variable
        forced_device = os.environ.get("DBNSFP_DEVICE", "").lower()
        if forced_device in ("cpu", "cuda", "mps"):
            device = forced_device
            print(f"Using {device.upper()} (set via DBNSFP_DEVICE)")
        elif torch.backends.mps.is_available():
            device = "mps"
            print("Using MPS (Apple Silicon GPU)")
        elif torch.cuda.is_available():
            device = "cuda"
            print("Using CUDA GPU")
        else:
            device = "cpu"
            print("Using CPU")

        self.embedder = SentenceTransformer(model_name, device=device)

        # Load or create index
        if self.index_path.exists():
            print(f"Loading existing index from {self.index_path}")
            self.index = faiss.read_index(str(self.index_path))
            with open(self.meta_path, "rb") as f:
                self.metadata = pickle.load(f)
        else:
            print("Creating new FAISS index")
            # Use IndexFlatIP for cosine similarity (normalize vectors first)
            self.index = faiss.IndexFlatIP(EMBEDDING_DIMENSION)
            self.metadata = {
                "ids": [],
                "documents": [],
                "metas": [],
            }

    def count(self) -> int:
        """Return number of variants in the index."""
        return self.index.ntotal

    def save(self):
        """Save index and metadata to disk."""
        faiss.write_index(self.index, str(self.index_path))
        with open(self.meta_path, "wb") as f:
            pickle.dump(self.metadata, f)

    def add_variants(
        self,
        ids: list[str],
        texts: list[str],
        metadatas: list[dict],
        batch_size: int = EMBEDDING_BATCH_SIZE,
        show_progress: bool = True,
    ):
        """
        Embed and add variants to FAISS.
        """
        total = len(ids)

        # Generate embeddings
        if show_progress:
            print(f"  Embedding {total} variants...")

        embeddings = self.embedder.encode(
            texts,
            show_progress_bar=show_progress,
            convert_to_numpy=True,
            batch_size=batch_size,
            normalize_embeddings=True,  # For cosine similarity
        )

        # Add to FAISS index
        self.index.add(embeddings.astype(np.float32))

        # Store metadata
        self.metadata["ids"].extend(ids)
        self.metadata["documents"].extend(texts)
        self.metadata["metas"].extend(metadatas)

        # Save after each batch
        self.save()

    def get_by_id(self, variant_id: str) -> dict | None:
        """Get a variant by its exact ID."""
        try:
            idx = self.metadata["ids"].index(variant_id)
            return {
                "id": variant_id,
                "document": self.metadata["documents"][idx],
                "metadata": self.metadata["metas"][idx],
            }
        except ValueError:
            return None

    def get_by_coordinates(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> dict | None:
        """Get a variant by genomic coordinates."""
        chrom = str(chrom).replace("chr", "")
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        return self.get_by_id(variant_id)

    def semantic_search(
        self,
        query: str,
        k: int = 10,
        gene_filter: str | None = None,
        clinvar_filter: str | None = None,
        min_cadd: float | None = None,
    ) -> list[dict]:
        """
        Semantic search for variants using vector similarity.
        """
        # Embed query
        query_embedding = self.embedder.encode(
            query,
            convert_to_numpy=True,
            normalize_embeddings=True,
        ).reshape(1, -1).astype(np.float32)

        # Search more than k to allow for filtering
        search_k = k * 10 if any([gene_filter, clinvar_filter, min_cadd]) else k
        search_k = min(search_k, self.count())

        if search_k == 0:
            return []

        distances, indices = self.index.search(query_embedding, search_k)

        # Filter and format results
        results = []
        for dist, idx in zip(distances[0], indices[0]):
            if idx < 0:
                continue

            meta = self.metadata["metas"][idx]

            # Apply filters
            if gene_filter and meta.get("gene", "").upper() != gene_filter.upper():
                continue
            if clinvar_filter and clinvar_filter.lower() not in meta.get("clinvar_sig", "").lower():
                continue
            if min_cadd is not None:
                cadd = meta.get("cadd_phred")
                if cadd is None or cadd < min_cadd:
                    continue

            results.append({
                "id": self.metadata["ids"][idx],
                "document": self.metadata["documents"][idx],
                "metadata": meta,
                "distance": float(1 - dist),  # Convert similarity to distance
            })

            if len(results) >= k:
                break

        return results

    def search_by_gene(self, gene: str, k: int = 100) -> list[dict]:
        """Get all variants for a gene."""
        gene_upper = gene.upper()
        results = []

        for i, meta in enumerate(self.metadata["metas"]):
            if meta.get("gene", "").upper() == gene_upper:
                results.append({
                    "id": self.metadata["ids"][i],
                    "document": self.metadata["documents"][i],
                    "metadata": meta,
                })
                if len(results) >= k:
                    break

        return results

    def search_by_region(
        self,
        chrom: str,
        start: int,
        end: int,
        k: int = 1000,
    ) -> list[dict]:
        """Get variants in a genomic region."""
        chrom = str(chrom).replace("chr", "")
        results = []

        for i, meta in enumerate(self.metadata["metas"]):
            if meta.get("chr") == chrom:
                pos = meta.get("pos", 0)
                if start <= pos <= end:
                    results.append({
                        "id": self.metadata["ids"][i],
                        "document": self.metadata["documents"][i],
                        "metadata": meta,
                    })
                    if len(results) >= k:
                        break

        return sorted(results, key=lambda x: x["metadata"].get("pos", 0))

    def delete_all(self):
        """Delete all variants (use with caution)."""
        self.index = faiss.IndexFlatIP(EMBEDDING_DIMENSION)
        self.metadata = {"ids": [], "documents": [], "metas": []}
        self.save()
        print("Deleted all variants")

    def optimize(self):
        """No-op for FAISS (already optimized)."""
        pass

"""RAG module - retrieval and LLM integration."""

import json
import re
from pathlib import Path

import ollama
import requests

from .config import DEFAULT_MODEL, OLLAMA_URL, VECTORDB_DIR
from .vectorstore import VariantVectorStore


SYSTEM_PROMPT = """You are a clinical genomics expert assistant at Alberta Precision Labs.
Your role is to help interpret genetic variants using evidence from dbNSFP annotations.

When analyzing variants:
1. Cite specific prediction scores and their interpretations
2. Consider the concordance between different predictors
3. Note population frequencies and their clinical significance
4. Reference ClinVar classifications when available
5. Consider conservation scores for evolutionary context
6. Be clear about uncertainty and conflicting evidence

Always base your interpretations on the provided variant data. If information is missing,
acknowledge this rather than speculating."""


class VariantRAG:
    """RAG system for variant annotation queries."""

    def __init__(
        self,
        db_path: Path | None = None,
        ollama_url: str | None = None,
        model: str | None = None,
    ):
        self.vectorstore = VariantVectorStore(db_path=db_path)
        self.ollama_url = ollama_url or OLLAMA_URL
        self.model = model or DEFAULT_MODEL

    def _check_ollama(self) -> bool:
        """Check if Ollama is running."""
        try:
            response = requests.get(f"{self.ollama_url}/api/tags", timeout=5)
            return response.status_code == 200
        except requests.exceptions.ConnectionError:
            return False

    def _query_ollama(self, prompt: str, system: str = SYSTEM_PROMPT) -> str:
        """Query local Ollama instance."""
        try:
            response = ollama.chat(
                model=self.model,
                messages=[
                    {"role": "system", "content": system},
                    {"role": "user", "content": prompt},
                ],
            )
            return response["message"]["content"]
        except Exception as e:
            return f"Error querying LLM: {e}"

    def exact_lookup(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> dict | None:
        """
        Look up a specific variant by coordinates.

        Args:
            chrom: Chromosome (e.g., "17" or "chr17")
            pos: Position (1-based)
            ref: Reference allele
            alt: Alternate allele

        Returns:
            Variant data or None
        """
        return self.vectorstore.get_by_coordinates(chrom, pos, ref, alt)

    def parse_variant_string(self, variant_str: str) -> dict | None:
        """
        Parse variant string in common formats.

        Supports:
        - chr17:43092919:G:A
        - chr17:43092919 G>A
        - 17-43092919-G-A
        - chr17:g.43092919G>A

        Returns:
            Dict with chr, pos, ref, alt or None if parsing fails
        """
        # Normalize
        s = variant_str.strip().upper()

        # Try different patterns
        patterns = [
            # chr17:43092919:G:A or 17:43092919:G:A
            r"(?:CHR)?(\w+):(\d+):([ACGT]+):([ACGT]+)",
            # chr17:43092919 G>A
            r"(?:CHR)?(\w+):(\d+)\s*([ACGT]+)>([ACGT]+)",
            # 17-43092919-G-A
            r"(?:CHR)?(\w+)-(\d+)-([ACGT]+)-([ACGT]+)",
            # 13_32326240_A_C (internal ID format with underscores)
            r"(?:CHR)?(\w+)_(\d+)_([ACGT]+)_([ACGT]+)",
            # chr17:g.43092919G>A (HGVS-like)
            r"(?:CHR)?(\w+):G\.(\d+)([ACGT]+)>([ACGT]+)",
        ]

        for pattern in patterns:
            match = re.match(pattern, s)
            if match:
                return {
                    "chr": match.group(1).replace("CHR", ""),
                    "pos": int(match.group(2)),
                    "ref": match.group(3),
                    "alt": match.group(4),
                }

        return None

    def lookup_variant(self, variant_str: str) -> dict | None:
        """Look up variant from string representation."""
        parsed = self.parse_variant_string(variant_str)
        if not parsed:
            return None
        return self.exact_lookup(
            chrom=parsed["chr"],
            pos=parsed["pos"],
            ref=parsed["ref"],
            alt=parsed["alt"],
        )

    def search(
        self,
        query: str,
        k: int = 5,
        gene: str | None = None,
    ) -> list[dict]:
        """
        Semantic search for variants.

        Args:
            query: Natural language query
            k: Number of results
            gene: Optional gene filter

        Returns:
            List of matching variants
        """
        return self.vectorstore.semantic_search(query, k=k, gene_filter=gene)

    def annotate(
        self,
        query: str,
        k: int = 3,
        gene: str | None = None,
    ) -> str:
        """
        Full RAG pipeline: retrieve relevant variants and generate response.

        Args:
            query: User question about variant(s)
            k: Number of context variants to retrieve
            gene: Optional gene filter

        Returns:
            LLM-generated response with citations
        """
        # Check if query contains a specific variant
        parsed = self.parse_variant_string(query)
        if parsed:
            # Direct lookup
            result = self.exact_lookup(**parsed)
            if result:
                context = result["document"]
            else:
                return f"Variant {query} not found in database."
        else:
            # Semantic search
            results = self.search(query, k=k, gene=gene)
            if not results:
                return "No relevant variants found for your query."
            context = "\n\n---\n\n".join([r["document"] for r in results])

        # Build prompt
        prompt = f"""Based on the following variant annotation data from dbNSFP, please answer the question.

VARIANT DATA:
{context}

QUESTION: {query}

Provide a clear, evidence-based response. Cite specific scores and predictions from the data."""

        # Query LLM
        if not self._check_ollama():
            return (
                "Ollama is not running. Start it with: ollama serve\n\n"
                f"Retrieved context:\n{context}"
            )

        return self._query_ollama(prompt)

    def interpret_variant(self, variant_str: str) -> str:
        """
        Get clinical interpretation of a specific variant.

        Args:
            variant_str: Variant in any supported format

        Returns:
            Clinical interpretation
        """
        result = self.lookup_variant(variant_str)
        if not result:
            return f"Variant {variant_str} not found in database."

        prompt = f"""Please provide a clinical interpretation of this variant.

VARIANT DATA:
{result['document']}

Provide:
1. Summary of pathogenicity predictions (concordance between predictors)
2. Population frequency interpretation
3. Clinical significance based on ClinVar (if available)
4. Conservation analysis
5. Overall assessment with confidence level"""

        if not self._check_ollama():
            return (
                "Ollama is not running. Start it with: ollama serve\n\n"
                f"Variant data:\n{result['document']}"
            )

        return self._query_ollama(prompt)

    def interpret(self, variant_str: str) -> str:
        """Alias for interpret_variant()."""
        return self.interpret_variant(variant_str)

    def compare_variants(self, variant_strs: list[str]) -> str:
        """
        Compare multiple variants.

        Args:
            variant_strs: List of variant strings

        Returns:
            Comparative analysis
        """
        variants = []
        for vs in variant_strs:
            result = self.lookup_variant(vs)
            if result:
                variants.append(result["document"])
            else:
                variants.append(f"Variant {vs}: NOT FOUND")

        context = "\n\n===\n\n".join(variants)

        prompt = f"""Compare the following variants:

{context}

Provide a comparative analysis covering:
1. Pathogenicity prediction comparison
2. Population frequency differences
3. Clinical significance comparison
4. Which variant(s) appear more likely to be pathogenic and why"""

        if not self._check_ollama():
            return (
                "Ollama is not running. Start it with: ollama serve\n\n"
                f"Variant data:\n{context}"
            )

        return self._query_ollama(prompt)


def main():
    """CLI for RAG queries."""
    import argparse

    parser = argparse.ArgumentParser(description="Query dbNSFP RAG system")
    parser.add_argument("query", nargs="?", help="Query or variant to look up")
    parser.add_argument("--lookup", "-l", help="Look up specific variant")
    parser.add_argument("--interpret", "-i", help="Interpret specific variant")
    parser.add_argument("--search", "-s", help="Semantic search query")
    parser.add_argument("--gene", "-g", help="Filter by gene")
    parser.add_argument("--k", type=int, default=5, help="Number of results")
    parser.add_argument("--model", default=DEFAULT_MODEL, help="Ollama model")
    parser.add_argument(
        "--compare", "-c", nargs="+", help="Compare multiple variants"
    )

    args = parser.parse_args()

    rag = VariantRAG(model=args.model)

    if args.lookup:
        result = rag.lookup_variant(args.lookup)
        if result:
            print(result["document"])
        else:
            print(f"Variant {args.lookup} not found")

    elif args.interpret:
        print(rag.interpret_variant(args.interpret))

    elif args.search:
        results = rag.search(args.search, k=args.k, gene=args.gene)
        for i, r in enumerate(results, 1):
            print(f"\n--- Result {i} (distance: {r['distance']:.4f}) ---")
            print(r["document"][:500] + "..." if len(r["document"]) > 500 else r["document"])

    elif args.compare:
        print(rag.compare_variants(args.compare))

    elif args.query:
        print(rag.annotate(args.query, k=args.k, gene=args.gene))

    else:
        # Interactive mode
        print("dbNSFP RAG System")
        print(f"Loaded {rag.vectorstore.count()} variants")
        print("Commands: /lookup <variant>, /interpret <variant>, /gene <name>, /quit")
        print()

        while True:
            try:
                query = input("Query> ").strip()
            except (EOFError, KeyboardInterrupt):
                break

            if not query:
                continue

            if query.lower() in ("/quit", "/exit", "/q"):
                break

            if query.startswith("/lookup "):
                variant = query[8:].strip()
                result = rag.lookup_variant(variant)
                if result:
                    print(result["document"])
                else:
                    print(f"Not found: {variant}")

            elif query.startswith("/interpret "):
                variant = query[11:].strip()
                print(rag.interpret_variant(variant))

            elif query.startswith("/gene "):
                gene = query[6:].strip()
                results = rag.vectorstore.search_by_gene(gene, k=10)
                print(f"Found {len(results)} variants in {gene}")
                for r in results[:5]:
                    print(f"  - {r['id']}")

            else:
                print(rag.annotate(query))

            print()


if __name__ == "__main__":
    main()

"""Export utilities for training data extraction."""

import json
from pathlib import Path

from .config import VECTORDB_DIR
from .panels import get_panel, get_all_panel_genes, PANELS
from .vectorstore import VariantVectorStore


def export_panel_variants(
    panel_name: str | None = None,
    genes: set[str] | None = None,
    output_dir: Path | None = None,
    format: str = "jsonl",
) -> Path:
    """
    Export variants for a gene panel to training-ready format.

    Args:
        panel_name: Name of predefined panel (from panels.py)
        genes: Custom set of genes (overrides panel_name)
        output_dir: Output directory (default: data/exports/)
        format: Output format - 'jsonl', 'csv', or 'parquet'

    Returns:
        Path to exported file
    """
    output_dir = output_dir or Path("data/exports")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get genes
    if genes:
        target_genes = genes
        name = "custom"
    elif panel_name:
        target_genes = get_panel(panel_name)
        name = panel_name
    else:
        target_genes = get_all_panel_genes()
        name = "all_panels"

    print(f"Exporting variants for {len(target_genes)} genes...")

    # Load vector store
    store = VariantVectorStore()

    # Collect variants
    variants = []
    for gene in sorted(target_genes):
        results = store.search_by_gene(gene, k=10000)  # Get all variants
        for r in results:
            variants.append({
                "id": r["id"],
                "gene": gene,
                "document": r["document"],
                "metadata": r["metadata"],
            })
        if results:
            print(f"  {gene}: {len(results)} variants")

    print(f"Total: {len(variants)} variants")

    # Export
    if format == "jsonl":
        output_path = output_dir / f"{name}_variants.jsonl"
        with open(output_path, "w") as f:
            for v in variants:
                f.write(json.dumps(v) + "\n")

    elif format == "csv":
        import csv
        output_path = output_dir / f"{name}_variants.csv"
        if variants:
            fieldnames = ["id", "gene"] + list(variants[0]["metadata"].keys())
            with open(output_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for v in variants:
                    row = {"id": v["id"], "gene": v["gene"]}
                    row.update(v["metadata"])
                    writer.writerow(row)

    else:
        raise ValueError(f"Unknown format: {format}")

    print(f"Exported to: {output_path}")
    return output_path


def export_for_classifier(
    panel_name: str | None = None,
    output_dir: Path | None = None,
) -> Path:
    """
    Export features for pathogenicity classifier training.

    Extracts numerical features (CADD, REVEL, etc.) and labels (ClinVar).
    """
    output_dir = output_dir or Path("data/exports")
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = get_panel(panel_name) if panel_name else get_all_panel_genes()
    store = VariantVectorStore()

    rows = []
    for gene in sorted(genes):
        results = store.search_by_gene(gene, k=10000)
        for r in results:
            meta = r["metadata"]
            # Extract numerical features
            row = {
                "id": r["id"],
                "gene": gene,
                "sift_score": meta.get("sift_score"),
                "polyphen2_hdiv": meta.get("polyphen2_hdiv"),
                "polyphen2_hvar": meta.get("polyphen2_hvar"),
                "cadd_phred": meta.get("cadd_phred"),
                "revel": meta.get("revel"),
                "alphamissense": meta.get("alphamissense"),
                "clinpred": meta.get("clinpred"),
                "dann": meta.get("dann"),
                "phylop100": meta.get("phylop100"),
                "phastcons100": meta.get("phastcons100"),
                "gerp": meta.get("gerp"),
                "gnomad_af": meta.get("gnomad_exome_af"),
                # Label
                "clinvar_sig": meta.get("clinvar_sig"),
            }
            rows.append(row)

    output_path = output_dir / f"{panel_name or 'all'}_classifier_features.jsonl"
    with open(output_path, "w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")

    print(f"Exported {len(rows)} variants for classifier training")
    print(f"Output: {output_path}")
    return output_path


def export_for_llm_finetuning(
    panel_name: str | None = None,
    output_dir: Path | None = None,
) -> Path:
    """
    Export instruction-response pairs for LLM fine-tuning.

    Format: {"instruction": "Interpret variant X", "response": "Clinical interpretation..."}
    """
    output_dir = output_dir or Path("data/exports")
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = get_panel(panel_name) if panel_name else get_all_panel_genes()
    store = VariantVectorStore()

    pairs = []
    for gene in sorted(genes):
        results = store.search_by_gene(gene, k=10000)
        for r in results:
            # Create instruction-response pair
            instruction = f"Interpret the clinical significance of variant {r['id']} in gene {gene}."

            # Build response from metadata
            meta = r["metadata"]
            response_parts = [f"Variant {r['id']} in {gene}:"]

            if meta.get("clinvar_sig"):
                response_parts.append(f"ClinVar classification: {meta['clinvar_sig']}")

            if meta.get("cadd_phred"):
                score = meta["cadd_phred"]
                interp = "high" if float(score) > 20 else "moderate" if float(score) > 15 else "low"
                response_parts.append(f"CADD score: {score} ({interp} predicted deleteriousness)")

            if meta.get("gnomad_exome_af"):
                af = meta["gnomad_exome_af"]
                response_parts.append(f"Population frequency (gnomAD): {af}")

            pairs.append({
                "instruction": instruction,
                "input": r["document"],  # Full variant context
                "output": "\n".join(response_parts),
            })

    output_path = output_dir / f"{panel_name or 'all'}_llm_training.jsonl"
    with open(output_path, "w") as f:
        for pair in pairs:
            f.write(json.dumps(pair) + "\n")

    print(f"Exported {len(pairs)} instruction-response pairs for LLM fine-tuning")
    print(f"Output: {output_path}")
    return output_path


def export_for_embeddings(
    panel_name: str | None = None,
    output_dir: Path | None = None,
) -> Path:
    """
    Export text pairs for embedding model fine-tuning.

    Format: {"anchor": "variant text", "positive": "similar variant", "negative": "different variant"}
    """
    output_dir = output_dir or Path("data/exports")
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = get_panel(panel_name) if panel_name else get_all_panel_genes()
    store = VariantVectorStore()

    # Collect variants by classification for contrastive pairs
    by_class = {"pathogenic": [], "benign": [], "uncertain": []}

    for gene in sorted(genes):
        results = store.search_by_gene(gene, k=10000)
        for r in results:
            sig = (r["metadata"].get("clinvar_sig") or "").lower()
            if "pathogenic" in sig and "conflict" not in sig:
                by_class["pathogenic"].append(r)
            elif "benign" in sig and "conflict" not in sig:
                by_class["benign"].append(r)
            else:
                by_class["uncertain"].append(r)

    # Create contrastive pairs
    import random
    pairs = []

    for path_var in by_class["pathogenic"][:1000]:  # Limit for training
        if by_class["benign"]:
            neg = random.choice(by_class["benign"])
            pairs.append({
                "anchor": path_var["document"],
                "positive": random.choice(by_class["pathogenic"])["document"],
                "negative": neg["document"],
                "label": "pathogenic",
            })

    for benign_var in by_class["benign"][:1000]:
        if by_class["pathogenic"]:
            neg = random.choice(by_class["pathogenic"])
            pairs.append({
                "anchor": benign_var["document"],
                "positive": random.choice(by_class["benign"])["document"],
                "negative": neg["document"],
                "label": "benign",
            })

    output_path = output_dir / f"{panel_name or 'all'}_embedding_pairs.jsonl"
    with open(output_path, "w") as f:
        for pair in pairs:
            f.write(json.dumps(pair) + "\n")

    print(f"Exported {len(pairs)} contrastive pairs for embedding fine-tuning")
    print(f"  Pathogenic variants: {len(by_class['pathogenic'])}")
    print(f"  Benign variants: {len(by_class['benign'])}")
    print(f"  Uncertain variants: {len(by_class['uncertain'])}")
    print(f"Output: {output_path}")
    return output_path


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Export training data")
    parser.add_argument("--panel", "-p", help="Panel name")
    parser.add_argument(
        "--type", "-t",
        choices=["all", "classifier", "llm", "embeddings"],
        default="all",
        help="Export type"
    )
    args = parser.parse_args()

    if args.type == "classifier" or args.type == "all":
        export_for_classifier(args.panel)

    if args.type == "llm" or args.type == "all":
        export_for_llm_finetuning(args.panel)

    if args.type == "embeddings" or args.type == "all":
        export_for_embeddings(args.panel)

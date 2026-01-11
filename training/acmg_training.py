#!/usr/bin/env python
"""
Generate ACMG classification training data for LLM fine-tuning.

Creates instruction-tuning examples that teach the model to:
1. Interpret variant features (CADD, REVEL, population frequency, etc.)
2. Apply ACMG/AMP criteria (PVS, PS, PM, PP, BS, BP, BA)
3. Output 5-tier classification with reasoning

Usage:
    python training/acmg_training.py --db data/vectordb/grch37-ngsgenes
"""

import json
import argparse
from pathlib import Path
from typing import Literal

# ACMG 5-tier classification
ACMGClass = Literal["Pathogenic", "Likely_pathogenic", "Uncertain_significance", "Likely_benign", "Benign"]


def map_clinvar_to_acmg(clinvar_sig: str) -> ACMGClass | None:
    """Map ClinVar significance to ACMG 5-tier classification."""
    if not clinvar_sig:
        return None

    sig = clinvar_sig.lower()

    # Skip conflicting/uncertain for training (want clean labels)
    if "conflict" in sig:
        return None
    if "not_provided" in sig:
        return None

    if "pathogenic/likely_pathogenic" in sig:
        return "Likely_pathogenic"  # Conservative
    if "benign/likely_benign" in sig:
        return "Likely_benign"  # Conservative
    if "pathogenic" in sig and "likely" not in sig:
        return "Pathogenic"
    if "likely_pathogenic" in sig:
        return "Likely_pathogenic"
    if "benign" in sig and "likely" not in sig:
        return "Benign"
    if "likely_benign" in sig:
        return "Likely_benign"
    if "uncertain" in sig:
        return "Uncertain_significance"

    return None


def generate_acmg_criteria(meta: dict, acmg_class: ACMGClass) -> list[str]:
    """Generate ACMG criteria based on variant evidence."""
    criteria = []

    # Population frequency criteria
    gnomad_af = meta.get("gnomad_af")
    if gnomad_af:
        try:
            af = float(gnomad_af)
            if af > 0.05:
                criteria.append("BA1: Allele frequency >5% in population databases")
            elif af > 0.01:
                criteria.append("BS1: Allele frequency greater than expected for disorder")
            elif af < 0.0001:
                criteria.append("PM2: Absent or extremely low frequency in population databases")
        except (ValueError, TypeError):
            pass
    else:
        if acmg_class in ("Pathogenic", "Likely_pathogenic"):
            criteria.append("PM2: Absent from population databases")

    # Computational evidence
    cadd = meta.get("cadd_phred")
    if cadd:
        try:
            cadd_score = float(cadd)
            if cadd_score >= 25:
                criteria.append(f"PP3: Computational evidence supports deleterious effect (CADD={cadd_score:.1f})")
            elif cadd_score < 15:
                criteria.append(f"BP4: Computational evidence suggests no impact (CADD={cadd_score:.1f})")
        except (ValueError, TypeError):
            pass

    revel = meta.get("revel_score")
    if revel:
        try:
            revel_score = float(revel)
            if revel_score >= 0.7:
                criteria.append(f"PP3: REVEL score {revel_score:.3f} predicts pathogenic")
            elif revel_score < 0.3:
                criteria.append(f"BP4: REVEL score {revel_score:.3f} predicts benign")
        except (ValueError, TypeError):
            pass

    # Prediction algorithms
    sift = meta.get("sift_pred", "")
    if "D" in str(sift):
        criteria.append("PP3: SIFT predicts damaging")
    elif "T" in str(sift):
        criteria.append("BP4: SIFT predicts tolerated")

    polyphen = meta.get("polyphen_pred", "")
    if "D" in str(polyphen):
        criteria.append("PP3: PolyPhen-2 predicts probably damaging")
    elif "B" in str(polyphen):
        criteria.append("BP4: PolyPhen-2 predicts benign")

    alphamissense = meta.get("alphamissense_pred", "")
    if "pathogenic" in str(alphamissense).lower():
        criteria.append("PP3: AlphaMissense predicts pathogenic")
    elif "benign" in str(alphamissense).lower():
        criteria.append("BP4: AlphaMissense predicts benign")

    return criteria


def format_variant_input(variant_id: str, gene: str, meta: dict) -> str:
    """Format variant features as model input."""
    parts = variant_id.split("_")
    if len(parts) >= 4:
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
    else:
        chrom, pos, ref, alt = "?", "?", "?", "?"

    lines = [
        f"Variant: chr{chrom}:{pos} {ref}>{alt}",
        f"Gene: {gene}",
    ]

    # Add available scores
    if meta.get("cadd_phred"):
        lines.append(f"CADD phred: {meta['cadd_phred']}")
    if meta.get("revel_score"):
        lines.append(f"REVEL: {meta['revel_score']}")
    if meta.get("gnomad_af"):
        lines.append(f"gnomAD AF: {meta['gnomad_af']}")
    if meta.get("sift_pred"):
        lines.append(f"SIFT: {meta['sift_pred']}")
    if meta.get("polyphen_pred"):
        lines.append(f"PolyPhen-2: {meta['polyphen_pred']}")
    if meta.get("alphamissense_pred"):
        lines.append(f"AlphaMissense: {meta['alphamissense_pred']}")

    return "\n".join(lines)


def format_acmg_output(acmg_class: ACMGClass, criteria: list[str], gene: str) -> str:
    """Format ACMG classification as model output."""
    lines = [
        f"ACMG Classification: {acmg_class.replace('_', ' ')}",
        "",
        "Evidence:",
    ]

    if criteria:
        for c in criteria:
            lines.append(f"- {c}")
    else:
        lines.append("- Limited computational evidence available")

    # Add classification reasoning
    lines.append("")
    if acmg_class == "Pathogenic":
        lines.append(f"Conclusion: This variant in {gene} meets criteria for pathogenic classification based on strong evidence of deleteriousness.")
    elif acmg_class == "Likely_pathogenic":
        lines.append(f"Conclusion: This variant in {gene} meets criteria for likely pathogenic classification based on moderate evidence of deleteriousness.")
    elif acmg_class == "Likely_benign":
        lines.append(f"Conclusion: This variant in {gene} meets criteria for likely benign classification based on evidence suggesting no functional impact.")
    elif acmg_class == "Benign":
        lines.append(f"Conclusion: This variant in {gene} meets criteria for benign classification based on strong evidence of no pathogenicity.")
    else:
        lines.append(f"Conclusion: This variant in {gene} has insufficient or conflicting evidence for definitive classification.")

    return "\n".join(lines)


def generate_training_data(db_path: Path, output_path: Path, max_per_class: int = 5000):
    """Generate ACMG training data from vector database."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from src.vectorstore import VariantVectorStore
    from src.panels import get_panel

    store = VariantVectorStore(db_path=db_path)
    genes = get_panel("NGSgenes")

    # Collect by class for balanced training
    by_class: dict[ACMGClass, list] = {
        "Pathogenic": [],
        "Likely_pathogenic": [],
        "Uncertain_significance": [],
        "Likely_benign": [],
        "Benign": [],
    }

    print(f"Processing {len(genes)} genes...")

    for gene in sorted(genes):
        results = store.search_by_gene(gene, k=10000)
        for r in results:
            meta = r["metadata"]
            clinvar = meta.get("clinvar_sig", "")

            acmg_class = map_clinvar_to_acmg(clinvar)
            if not acmg_class:
                continue

            if len(by_class[acmg_class]) >= max_per_class:
                continue

            criteria = generate_acmg_criteria(meta, acmg_class)

            example = {
                "instruction": "Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria.",
                "input": format_variant_input(r["id"], gene, meta),
                "output": format_acmg_output(acmg_class, criteria, gene),
                "acmg_class": acmg_class,
                "variant_id": r["id"],
                "gene": gene,
            }

            by_class[acmg_class].append(example)

    # Report distribution
    print("\nACMG class distribution:")
    total = 0
    for cls, examples in by_class.items():
        print(f"  {cls}: {len(examples)}")
        total += len(examples)
    print(f"  Total: {total}")

    # Combine and shuffle
    import random
    all_examples = []
    for examples in by_class.values():
        all_examples.extend(examples)
    random.shuffle(all_examples)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # JSONL format (for training)
    jsonl_path = output_path.with_suffix(".jsonl")
    with open(jsonl_path, "w") as f:
        for ex in all_examples:
            f.write(json.dumps(ex) + "\n")
    print(f"\nSaved JSONL: {jsonl_path}")

    # Alpaca JSON format (for mlx-lm)
    alpaca_path = output_path.with_suffix(".json")
    alpaca_data = [
        {"instruction": ex["instruction"], "input": ex["input"], "output": ex["output"]}
        for ex in all_examples
    ]
    with open(alpaca_path, "w") as f:
        json.dump(alpaca_data, f, indent=2)
    print(f"Saved Alpaca JSON: {alpaca_path}")

    return all_examples


def main():
    parser = argparse.ArgumentParser(description="Generate ACMG training data")
    parser.add_argument("--db", type=Path, default=Path("data/vectordb/grch37-ngsgenes"),
                       help="Path to vector database")
    parser.add_argument("--output", "-o", type=Path, default=Path("data/training/acmg_training"),
                       help="Output path (without extension)")
    parser.add_argument("--max-per-class", type=int, default=5000,
                       help="Maximum examples per ACMG class")
    args = parser.parse_args()

    generate_training_data(args.db, args.output, args.max_per_class)


if __name__ == "__main__":
    main()

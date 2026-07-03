#!/usr/bin/env python3
"""Cross-check the rule-based classifier against the fine-tuned ACMG LLM.

For each variant in the upload set, build the prompt the model was trained
on (gene + scores), run inference, parse the output for the 5-tier ACMG
classification, and compare against our rule-based call.

Reports:
  - agreement matrix (our class x LLM class)
  - rows where LLM disagrees with us (these are the ones to manually review)
  - per-source agreement (pipeline vs ClinVar vs synonymous_catalog)

The model is base Llama-3.2-3B-Instruct-4bit at models/acmg-classifier/model
plus LoRA adapters at models/acmg-classifier/adapters (1000 training iters
on NGSgenes ACMG examples).

Usage:
  uv run python scripts/llm_crosscheck.py                          # default: 5,796 upload
  uv run python scripts/llm_crosscheck.py --limit 50                # quick sanity
  uv run python scripts/llm_crosscheck.py --source pipeline         # filter
  uv run python scripts/llm_crosscheck.py --input some.tsv          # other input
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
from tqdm import tqdm

DEFAULT_INPUT = Path("cp_new_bundle/outputs/new_genes_only/cp_new_seqnext_UPLOAD.tsv")
DEFAULT_OUT = Path("data/exports/cp_new/seqnext/llm_crosscheck.tsv")
PER_GENE_DIR = Path("data/exports/cp_new")
MODEL_PATH = "models/acmg-classifier/model"
ADAPTER_PATH = "models/acmg-classifier/adapters"

INSTRUCTION = "Classify this variant according to ACMG/AMP guidelines and provide the evidence criteria."

CLASS_RE = re.compile(
    r"(Pathogenic|Likely[\s_]+pathogenic|Uncertain[\s_]+significance|Likely[\s_]+benign|Benign)",
    re.IGNORECASE,
)

ACMG_NORM = {
    "pathogenic": "Pathogenic",
    "likely pathogenic": "Likely_pathogenic", "likely_pathogenic": "Likely_pathogenic",
    "uncertain significance": "Uncertain_significance", "uncertain_significance": "Uncertain_significance",
    "likely benign": "Likely_benign", "likely_benign": "Likely_benign",
    "benign": "Benign",
}


def parse_llm_class(text: str) -> str | None:
    m = CLASS_RE.search(text)
    if not m: return None
    return ACMG_NORM.get(re.sub(r"\s+", " ", m.group(1).lower()), None)


def format_variant_input(gene: str, chr_: str, pos: str, ref: str, alt: str, meta: dict) -> str:
    """Match the training format used in training/acmg_training.py."""
    lines = [
        f"Variant: chr{chr_}:{pos} {ref}>{alt}",
        f"Gene: {gene}",
    ]
    fields = [
        ("CADD phred", meta.get("CADD_phred") or meta.get("CADD_phred_dbnsfp")),
        ("REVEL", meta.get("REVEL_score") or meta.get("REVEL")),
        ("gnomAD AF", meta.get("gnomad_v41_af_joint")),
        ("FAF95 grpmax", meta.get("gnomad_v41_faf95_grpmax") or meta.get("FAF95_grpmax")),
        ("SIFT", meta.get("SIFT_pred")),
        ("PolyPhen-2", meta.get("Polyphen2_HDIV_pred")),
        ("AlphaMissense", meta.get("AlphaMissense_pred") or meta.get("AlphaMissense_pred_dbnsfp")),
        ("MutationTaster", meta.get("MutationTaster_pred")),
        ("BayesDel", meta.get("BayesDel_addAF_pred")),
        ("SpliceAI masked", meta.get("spliceai_ds_max_masked") or meta.get("SpliceAI_masked")),
        ("PhastCons", meta.get("phastCons100way_vertebrate") or meta.get("PhastCons100way")),
        ("PhyloP", meta.get("phyloP100way_vertebrate") or meta.get("PhyloP100way")),
        ("ClinVar", meta.get("clinvar_clnsig") or meta.get("ClinVar_sig")),
    ]
    for label, val in fields:
        if val and str(val) not in ("nan", "", "."):
            lines.append(f"{label}: {val}")
    return "\n".join(lines)


def load_per_gene_cache(gene: str) -> dict[tuple, dict]:
    """Index per-gene dbNSFP TSV by (hg19_chr, hg19_pos, ref, alt)."""
    tsv = PER_GENE_DIR / f"{gene}.tsv"
    if not tsv.exists(): return {}
    try:
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".",""], low_memory=False)
    except Exception:
        return {}
    out = {}
    for _, r in df.iterrows():
        try:
            key = (str(r["hg19_chr"]),
                   str(int(float(r["hg19_pos(1-based)"]))),
                   r["ref"], r["alt"])
        except (TypeError, ValueError, KeyError):
            continue
        if key not in out:
            out[key] = r.to_dict()
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    p.add_argument("--limit", type=int, help="cap rows for quick test")
    p.add_argument("--source", help="filter by source column (pipeline/clinvar/synonymous_catalog)")
    p.add_argument("--max-tokens", type=int, default=80)
    args = p.parse_args()

    if not args.input.exists():
        sys.exit(f"Missing input: {args.input}")
    df = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    if args.source and "source" in df.columns:
        df = df[df["source"] == args.source]
    if args.limit:
        df = df.head(args.limit)
    print(f"Cross-checking {len(df):,} variants against the fine-tuned ACMG LLM\n")

    print(f"Loading model {MODEL_PATH} + adapters {ADAPTER_PATH} ...")
    from mlx_lm import load, generate
    t0 = time.time()
    model, tokenizer = load(MODEL_PATH, adapter_path=ADAPTER_PATH)
    print(f"  loaded in {time.time()-t0:.1f}s\n")

    # Per-gene cache loader (lazy)
    gene_cache: dict[str, dict] = {}
    def cache(gene):
        if gene not in gene_cache:
            gene_cache[gene] = load_per_gene_cache(gene)
        return gene_cache[gene]

    out_rows = []
    n_dbnsfp = 0  # how many had dbNSFP context for the prompt
    matrix = Counter()
    t0 = time.time()
    for i, r in tqdm(df.iterrows(), total=len(df), desc="llm"):
        gene = r["gene"]
        # Look up rich dbNSFP row if available
        key = (str(r.get("chr_grch37","")),
               str(r.get("pos_grch37","")),
               r.get("ref",""), r.get("alt",""))
        meta = cache(gene).get(key, {})
        # Merge with row's own scores so we still have something when dbNSFP misses
        meta = {**r.to_dict(), **meta}
        if any(meta.get(k) for k in ("REVEL_score","CADD_phred","phastCons100way_vertebrate")):
            n_dbnsfp += 1

        prompt_input = format_variant_input(
            gene, str(r.get("chr_grch37","")), str(r.get("pos_grch37","")),
            r.get("ref",""), r.get("alt",""), meta,
        )
        prompt = f"{INSTRUCTION}\n\n{prompt_input}\n\nACMG Classification:"

        text = generate(model, tokenizer, prompt=prompt, max_tokens=args.max_tokens, verbose=False)
        llm_class = parse_llm_class(text)

        # Our class (strip rule string from the full classification cell)
        our_full = r.get("classification","")
        our_class = "Benign" if our_full.startswith("Benign") else \
                    "Likely_benign" if our_full.startswith("Likely_benign") else \
                    "Pathogenic" if our_full.startswith("Pathogenic") else \
                    "Likely_pathogenic" if our_full.startswith("Likely_pathogenic") else "Uncertain_significance"
        agree = (llm_class == our_class) if llm_class else None
        matrix[(our_class, llm_class or "no_parse")] += 1

        out_rows.append({
            "gene": gene,
            "transcript": r.get("transcript",""),
            "hgvs_c": r.get("hgvs_c",""),
            "our_classification": our_full,
            "our_class": our_class,
            "llm_class": llm_class or "",
            "agree": agree,
            "source": r.get("source",""),
            "llm_response_head": text.split("\n",1)[0].strip()[:200],
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    dt = time.time() - t0
    print(f"\nWrote {len(out_df):,} rows -> {args.output}  ({dt:.0f}s total, {dt/max(len(out_df),1):.2f}s/variant)")
    print(f"With dbNSFP context: {n_dbnsfp:,} / {len(out_df):,}")

    print(f"\n=== Agreement matrix (our class x LLM class) ===")
    rows = sorted({k[0] for k in matrix})
    cols = sorted({k[1] for k in matrix})
    print(f"{'our':<20} | " + " ".join(f"{c:<22}" for c in cols))
    print("-" * (24 + 23 * len(cols)))
    for r_ in rows:
        line = f"{r_:<20} | "
        line += " ".join(f"{matrix.get((r_,c),0):<22}" for c in cols)
        print(line)

    n_agree = sum(1 for r in out_rows if r["agree"] is True)
    n_disagree = sum(1 for r in out_rows if r["agree"] is False)
    n_unparsed = sum(1 for r in out_rows if r["agree"] is None)
    print(f"\nAgreement: {n_agree:,} / {len(out_df):,} ({100*n_agree/max(len(out_df),1):.1f}%)")
    print(f"Disagree:  {n_disagree:,}  (rows where LLM called different than rules)")
    print(f"Unparsed:  {n_unparsed:,}  (LLM output didn't contain a recognized class)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

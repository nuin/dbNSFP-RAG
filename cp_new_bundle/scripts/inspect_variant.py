#!/usr/bin/env python3
"""Alamut-style variant inspector for the cp_new pipeline.

Pulls together everything we know (and everything we can fetch from free
public APIs) about a single variant, formatted as a multi-section report
the analyst can read in one screen.

Subcommands:
  lookup    -- single variant from gene+c. or chr:pos:ref:alt or rsID
  gene      -- dump all classified variants for a gene
  bulk      -- process a TSV list, write per-variant reports

Data sources (cascading local -> external, all free):
  Local cache (instant):
    - data/exports/cp_new/{GENE}.tsv             dbNSFP + gnomAD per-gene
    - data/exports/cp_new/seqnext/cp_new_seqnext_FINAL.tsv  our classifications
    - data/exports/cp_new/.vv_cache.json         VariantValidator cache
    - ~/data/clinvar/variant_summary.txt.gz      ClinVar bulk
    - ~/data/hg19/hg19.fa  ~/data/hg38/hg38.fa   reference FASTAs

  External APIs (rate-limited):
    - rest.variantvalidator.org    HGVS resolution on any transcript
    - myvariant.info               aggregator: dbSNP, ClinVar, CADD, EVS, COSMIC
    - rest.ensembl.org (+ grch37)  gene/transcript metadata
    - eutils.ncbi.nlm.nih.gov      ClinVar variant detail + dbSNP

Usage:
  inspect_variant.py lookup --gene GOT2 --hgvs c.816C>T
  inspect_variant.py lookup --coord 16:58750604:G:A --build hg19
  inspect_variant.py lookup --rsid rs1058192
  inspect_variant.py gene GOT2
  inspect_variant.py bulk variants.tsv --out reports/
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from inspector_report import (
    SECTION_TITLES,
    add_classifications,
    add_external_annotations,
    add_local_annotations,
    add_local_identity,
    render_variant,
)
from inspector_sources import (
    FINAL,
    Variant,
    fetch_myvariant,
    find_in_classifications,
    find_in_per_gene,
    parse_input,
)

try:
    from rich.console import Console
    from rich.markdown import Markdown

    HAS_RICH = True
except ImportError:
    HAS_RICH = False


def resolve_rsid(v: Variant, sections: dict) -> None:
    mv = fetch_myvariant(v)
    if not mv:
        return
    hg19 = mv.get("hg19")
    vcf = mv.get("vcf")
    if not isinstance(vcf, dict):
        sections["Notes"].append("- rsID response has no VCF coordinates.")
        return
    start = hg19.get("start") if isinstance(hg19, dict) else None
    try:
        pos = int(start or vcf.get("position"))
    except (TypeError, ValueError, OverflowError):
        sections["Notes"].append("- rsID response has no valid hg19 position.")
        return
    v.chr_hg19 = str(mv.get("chrom") or "").removeprefix("chr")
    v.pos_hg19 = pos
    v.ref = vcf.get("ref")
    v.alt = vcf.get("alt")


def cmd_lookup(args: argparse.Namespace) -> int:
    v = parse_input(args)
    print("# resolving variant ...", file=sys.stderr)
    sections: dict[str, list[str]] = {title: [] for title in SECTION_TITLES}
    if v.rsid and not (v.chr_hg19 and v.pos_hg19) and not args.no_external:
        resolve_rsid(v, sections)
    per_gene = find_in_per_gene(v) if v.gene else None
    if per_gene:
        add_local_identity(per_gene, sections)
        add_local_annotations(per_gene, sections)
    else:
        sections["Notes"].append("- *Not in our cp_new dbNSFP cache.*")
    classifications = find_in_classifications(v)
    add_classifications(classifications, sections)
    if not args.no_external:
        add_external_annotations(v, classifications, sections)
    report = render_variant(v, sections)
    if args.json:
        print(json.dumps({"variant": asdict(v), "sections": sections}, default=str, indent=2))
    elif HAS_RICH and sys.stdout.isatty():
        Console().print(Markdown(report))
    else:
        print(report)
    return 0


def cmd_gene(args: argparse.Namespace) -> int:
    gene = args.gene_symbol.upper()
    if not FINAL.exists():
        sys.exit(f"Missing classifications: {FINAL}")
    import pandas as pd

    df = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    sub = df[df["gene"] == gene]
    if sub.empty:
        print(f"No classifications for {gene} in {FINAL}", file=sys.stderr)
        return 1
    cols = [
        "gene",
        "transcript",
        "hgvs_c",
        "classification",
        "FAF95_grpmax",
        "REVEL",
        "SpliceAI_masked",
        "PhastCons100way",
        "ClinVar_sig",
        "source",
    ]
    cols = [c for c in cols if c in sub.columns]
    if args.format == "tsv":
        sub[cols].to_csv(sys.stdout, sep="\t", index=False)
    else:
        # markdown table
        print(f"# {gene} — {len(sub)} classified variants\n")
        print("| " + " | ".join(cols) + " |")
        print("|" + "|".join(["---"] * len(cols)) + "|")
        for _, r in sub.iterrows():
            print("| " + " | ".join(str(r.get(c, "")).replace("|", "\\|") for c in cols) + " |")
    return 0


def cmd_bulk(args: argparse.Namespace) -> int:
    inp = Path(args.input)
    if not inp.exists():
        sys.exit(f"Input not found: {inp}")
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    import pandas as pd

    df = pd.read_csv(inp, sep="\t", dtype=str).fillna("")
    print(f"Processing {len(df)} variants -> {out_dir}", file=sys.stderr)
    for i, r in df.iterrows():
        # Build a Variant from each row (expect gene + hgvs_c columns)
        v = Variant(gene=r.get("gene", ""), hgvs_c=r.get("hgvs_c", ""))
        if r.get("chr_grch37"):
            v.chr_hg19 = r["chr_grch37"]
            v.pos_hg19 = int(r["pos_grch37"]) if r.get("pos_grch37") else None
            v.ref = r["ref"]
            v.alt = r["alt"]
        # Render minimal lookup (no external by default in bulk; opt-in)
        ns = argparse.Namespace(
            coord=None,
            rsid=None,
            gene=v.gene,
            hgvs=v.hgvs_c,
            transcript=None,
            build="hg19",
            json=False,
            no_external=not args.external,
        )
        try:
            # Don't actually call cmd_lookup (which prints) -- just write to file
            from io import StringIO

            buf = StringIO()
            _stdout = sys.stdout
            sys.stdout = buf
            try:
                cmd_lookup(ns)
            finally:
                sys.stdout = _stdout
            fname = f"{v.gene}_{(v.hgvs_c or '').replace('>', '-').replace('/', '_')[:40]}.md"
            (out_dir / fname).write_text(buf.getvalue())
        except SystemExit:
            continue
        if (i + 1) % 50 == 0:
            print(f"  {i + 1}/{len(df)}", file=sys.stderr)
    print("done", file=sys.stderr)
    return 0


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("lookup", help="Single-variant report")
    pl.add_argument("--gene")
    pl.add_argument("--hgvs", help="c. notation, e.g. c.816C>T")
    pl.add_argument("--coord", help="chr:pos:ref:alt")
    pl.add_argument("--build", choices=["hg19", "hg38"], default="hg19")
    pl.add_argument("--rsid")
    pl.add_argument("--transcript", help="Force a specific NM_ for HGVS resolution")
    pl.add_argument("--no-external", action="store_true", help="Skip API calls (local cache only)")
    pl.add_argument("--json", action="store_true", help="Emit JSON instead of Markdown")
    pl.set_defaults(func=cmd_lookup)

    pg = sub.add_parser("gene", help="Dump all classified variants for a gene")
    pg.add_argument("gene_symbol")
    pg.add_argument("--format", choices=["md", "tsv"], default="md")
    pg.set_defaults(func=cmd_gene)

    pb = sub.add_parser("bulk", help="Process a TSV list of variants")
    pb.add_argument(
        "input", help="TSV with at least gene+hgvs_c (or chr_grch37+pos_grch37+ref+alt)"
    )
    pb.add_argument(
        "--out", default="reports/", help="Output directory for per-variant Markdown files"
    )
    pb.add_argument("--external", action="store_true", help="Include external API calls (slower)")
    pb.set_defaults(func=cmd_bulk)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

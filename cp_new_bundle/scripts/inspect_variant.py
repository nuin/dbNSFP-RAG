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
import re
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import requests

# Optional pretty output (falls back to plain text if rich isn't installed)
try:
    from rich.console import Console
    from rich.markdown import Markdown
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


# =============================================================================
# Paths to local caches
# =============================================================================

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "exports" / "cp_new"
PER_GENE_DIR = DATA
FINAL = DATA / "seqnext" / "cp_new_seqnext_FINAL.tsv"
VV_CACHE_PATH = DATA / ".vv_cache.json"
CLINVAR_BULK = Path.home() / "data" / "clinvar" / "variant_summary.txt.gz"
HG19_FA = Path.home() / "data" / "hg19" / "hg19.fa"
HG38_FA = Path.home() / "data" / "hg38" / "hg38.fa"


# =============================================================================
# Variant normalization
# =============================================================================

COORD_RE = re.compile(r"^(?:chr)?([\dXYM]+)[:-](\d+)[:-]([ACGTN-]+)[:-]([ACGTN-]+)$", re.I)
RSID_RE = re.compile(r"^rs\d+$", re.I)
HGVS_RE = re.compile(r"^(NM_\d+\.?\d*|ENST\d+\.?\d*)?:?(c\.[A-Za-z\d+\-*>\.]+)$")


@dataclass
class Variant:
    """Canonicalized variant identity. Any field may be filled lazily."""
    gene: Optional[str] = None
    transcript: Optional[str] = None
    hgvs_c: Optional[str] = None
    chr_hg19: Optional[str] = None
    pos_hg19: Optional[int] = None
    chr_hg38: Optional[str] = None
    pos_hg38: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    rsid: Optional[str] = None
    build_input: Optional[str] = None  # which build the user's input was in

    def key(self) -> str:
        if self.chr_hg19 and self.pos_hg19 and self.ref and self.alt:
            return f"{self.chr_hg19}-{self.pos_hg19}-{self.ref}-{self.alt}"
        if self.chr_hg38 and self.pos_hg38 and self.ref and self.alt:
            return f"hg38:{self.chr_hg38}-{self.pos_hg38}-{self.ref}-{self.alt}"
        return f"{self.gene or '?'}|{self.transcript or '?'}|{self.hgvs_c or '?'}"

    def coord_present(self) -> bool:
        return bool((self.chr_hg19 and self.pos_hg19) or (self.chr_hg38 and self.pos_hg38))


def parse_input(args: argparse.Namespace) -> Variant:
    """Normalize argparse args into a Variant. Multiple forms can be combined
    (e.g. --gene GOT2 --coord 16:X:G:A gives both hints)."""
    v = Variant()
    if args.gene:
        v.gene = args.gene.upper()
    if args.transcript:
        v.transcript = args.transcript
    if args.hgvs:
        v.hgvs_c = args.hgvs.strip()
    if args.coord:
        m = COORD_RE.match(args.coord.strip())
        if not m:
            sys.exit(f"--coord must look like chr:pos:ref:alt (got {args.coord!r})")
        chrom, pos, ref, alt = m.groups()
        if args.build == "hg38":
            v.chr_hg38, v.pos_hg38 = chrom, int(pos)
        else:
            v.chr_hg19, v.pos_hg19 = chrom, int(pos)
        v.ref, v.alt = ref.upper(), alt.upper()
        v.build_input = args.build
    if args.rsid:
        if not RSID_RE.match(args.rsid):
            sys.exit(f"--rsid must look like rsNNNN (got {args.rsid!r})")
        v.rsid = args.rsid.lower()
    if not any([v.gene, v.coord_present(), v.rsid, v.hgvs_c]):
        sys.exit("Need one of: --gene + --hgvs/--coord, --coord CHR:POS:REF:ALT, or --rsid rsX")
    return v


# =============================================================================
# Local fetcher
# =============================================================================

def load_vv_cache() -> dict:
    if not VV_CACHE_PATH.exists():
        return {}
    with open(VV_CACHE_PATH) as f:
        return json.load(f)


def find_in_per_gene(v: Variant) -> Optional[dict]:
    """Look up variant in per-gene dbNSFP+gnomAD TSV."""
    if not v.gene:
        return None
    tsv = PER_GENE_DIR / f"{v.gene}.tsv"
    if not tsv.exists():
        return None
    import pandas as pd
    df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".",""], low_memory=False)
    if v.chr_hg19 and v.pos_hg19:
        mask = (df["hg19_chr"].astype(str) == str(v.chr_hg19)) & \
               (df["hg19_pos(1-based)"].astype(str) == str(v.pos_hg19)) & \
               (df["ref"] == v.ref) & (df["alt"] == v.alt)
    elif v.chr_hg38 and v.pos_hg38:
        mask = (df["#chr"].astype(str) == str(v.chr_hg38)) & \
               (df["pos(1-based)"].astype(str) == str(v.pos_hg38)) & \
               (df["ref"] == v.ref) & (df["alt"] == v.alt)
    else:
        return None
    sub = df[mask]
    if sub.empty:
        return None
    return sub.iloc[0].to_dict()


def find_in_classifications(v: Variant) -> list[dict]:
    """Look up variant in cp_new_seqnext_FINAL.tsv. Tries successively looser
    matches: (gene, hgvs_c) -> (chr, pos, ref, alt) -> (gene, chr, pos).
    Returns the union of matches across strategies."""
    if not FINAL.exists():
        return []
    import pandas as pd
    df = pd.read_csv(FINAL, sep="\t", dtype=str).fillna("")
    masks = []
    if v.gene and v.hgvs_c:
        masks.append((df["gene"] == v.gene) &
                     (df["hgvs_c"] == v.hgvs_c))
    if v.chr_hg19 and v.pos_hg19 and v.ref and v.alt:
        masks.append((df["chr_grch37"] == str(v.chr_hg19)) &
                     (df["pos_grch37"] == str(v.pos_hg19)) &
                     (df["ref"] == v.ref) & (df["alt"] == v.alt))
    if v.gene and v.chr_hg19 and v.pos_hg19:
        # Loose: same gene + position (catches ClinVar rows with empty ref/alt)
        masks.append((df["gene"] == v.gene) &
                     (df["chr_grch37"] == str(v.chr_hg19)) &
                     (df["pos_grch37"] == str(v.pos_hg19)))
    if not masks:
        return []
    combined = masks[0]
    for m in masks[1:]:
        combined = combined | m
    return df[combined].drop_duplicates(subset=["gene","transcript","hgvs_c","source"]).to_dict("records")


# =============================================================================
# External fetchers (free APIs)
# =============================================================================

_HTTP_CACHE_DIR = REPO / "data" / "exports" / "cp_new" / ".inspector_cache"
_HTTP_CACHE_DIR.mkdir(parents=True, exist_ok=True)
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json", "User-Agent": "cp_new-inspector/0.1"})


def _cache_path(key: str) -> Path:
    safe = re.sub(r"[^a-zA-Z0-9._-]", "_", key)[:200]
    return _HTTP_CACHE_DIR / f"{safe}.json"


def http_get_json(url: str, cache_key: str, sleep_after: float = 0.0) -> Optional[dict]:
    """GET with on-disk cache. Returns None on error."""
    cp = _cache_path(cache_key)
    if cp.exists():
        try:
            return json.loads(cp.read_text())
        except Exception:
            pass
    try:
        r = SESSION.get(url, timeout=30)
        if r.status_code != 200:
            return None
        data = r.json()
        cp.write_text(json.dumps(data))
        if sleep_after:
            time.sleep(sleep_after)
        return data
    except (requests.RequestException, ValueError):
        return None


def fetch_variantvalidator(v: Variant, transcript: Optional[str] = None) -> Optional[dict]:
    """Resolve HGVS on a transcript via VariantValidator REST."""
    if not (v.chr_hg19 and v.pos_hg19 and v.ref and v.alt):
        return None
    tx = transcript or v.transcript or "select"
    url = (f"https://rest.variantvalidator.org/VariantValidator/variantvalidator/"
           f"GRCh37/{v.chr_hg19}-{v.pos_hg19}-{v.ref}-{v.alt}/{tx}")
    return http_get_json(url, f"vv_{v.key()}_{tx}", sleep_after=0.3)


def fetch_myvariant(v: Variant) -> Optional[dict]:
    """myvariant.info aggregator. Build hg19 query string from coords or HGVS."""
    if v.chr_hg19 and v.pos_hg19 and v.ref and v.alt:
        q = f"chr{v.chr_hg19}:g.{v.pos_hg19}{v.ref}>{v.alt}"
    elif v.rsid:
        q = v.rsid
    else:
        return None
    url = f"https://myvariant.info/v1/variant/{q}?assembly=hg19"
    return http_get_json(url, f"mv_{q}", sleep_after=0.1)


def fetch_ensembl_gene(gene_symbol: str) -> Optional[dict]:
    """Ensembl REST lookup for gene metadata in GRCh37."""
    url = f"https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    return http_get_json(url, f"ens37_gene_{gene_symbol}", sleep_after=0.1)


def fetch_clinvar_eutils(v: Variant) -> Optional[dict]:
    """ClinVar variation lookup via NCBI eutils -- gives submitter + assertion detail."""
    if not v.rsid:
        return None
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar"
           f"&term={v.rsid}&retmode=json")
    return http_get_json(url, f"clinvar_eutils_{v.rsid}", sleep_after=0.5)


# =============================================================================
# Renderer (Markdown)
# =============================================================================

def render_section(title: str, lines: list[str]) -> str:
    body = "\n".join(lines).strip()
    if not body:
        return ""
    return f"## {title}\n\n{body}\n"


def render_variant(v: Variant, sections: dict[str, list[str]]) -> str:
    """Compose the multi-section Markdown report."""
    out = [f"# Variant report — {v.gene or 'unknown'} {v.hgvs_c or v.key()}"]
    out.append(f"\n*Generated by `inspect_variant.py`. Sources cascading: local cache → public APIs.*\n")
    for title in ["Identity", "Frequencies (gnomAD)", "Predictors",
                  "Splice & conservation", "Clinical assertions (ClinVar)",
                  "Our pipeline classification", "Cross-references",
                  "Reference context", "Notes"]:
        s = render_section(title, sections.get(title, []))
        if s: out.append(s)
    return "\n".join(out)


def kv(label: str, value, fmt=str) -> str:
    if value is None or value == "":
        return f"- **{label}**: —"
    try:
        return f"- **{label}**: {fmt(value)}"
    except Exception:
        return f"- **{label}**: {value}"


# =============================================================================
# Subcommand: lookup
# =============================================================================

def cmd_lookup(args: argparse.Namespace) -> int:
    v = parse_input(args)
    print(f"# resolving variant ...", file=sys.stderr)

    sections: dict[str, list[str]] = {k: [] for k in
        ("Identity", "Frequencies (gnomAD)", "Predictors", "Splice & conservation",
         "Clinical assertions (ClinVar)", "Our pipeline classification",
         "Cross-references", "Reference context", "Notes")}

    # --- if input was rsID or HGVS-only, resolve via VariantValidator first ---
    if v.rsid and not (v.chr_hg19 and v.pos_hg19):
        # Resolve rsID via myvariant
        mv = fetch_myvariant(v)
        if mv:
            try:
                v.chr_hg19 = str(mv.get("chrom") or "").removeprefix("chr")
                v.pos_hg19 = int(mv.get("hg19", {}).get("start") or mv.get("vcf", {}).get("position"))
                v.ref = mv.get("vcf", {}).get("ref")
                v.alt = mv.get("vcf", {}).get("alt")
            except Exception:
                pass

    if v.gene and v.hgvs_c and not v.chr_hg19:
        # Use VariantValidator to resolve gene+c. -> coords
        vv = fetch_variantvalidator(Variant(gene=v.gene, hgvs_c=v.hgvs_c), transcript=v.transcript or "select")
        # (VV's coord parsing varies; skip for v1 -- user can pass --coord directly)

    # --- LOCAL: per-gene TSV (dbNSFP + gnomAD) ---
    per_gene = find_in_per_gene(v) if v.gene else None
    if per_gene:
        sections["Identity"] += [
            kv("gene", per_gene.get("genename","").split(";")[0]),
            kv("Ensembl tx (dbNSFP)", per_gene.get("Ensembl_transcriptid","").split(";")[0]),
            kv("HGVSc (dbNSFP/snpEff)", per_gene.get("HGVSc_snpEff","").split(";")[0]),
            kv("HGVSp (dbNSFP/snpEff)", per_gene.get("HGVSp_snpEff","").split(";")[0]),
            kv("hg19", f"chr{per_gene.get('hg19_chr')}:{per_gene.get('hg19_pos(1-based)')} {per_gene.get('ref')}>{per_gene.get('alt')}"),
            kv("hg38", f"chr{per_gene.get('#chr')}:{per_gene.get('pos(1-based)')} {per_gene.get('ref')}>{per_gene.get('alt')}"),
        ]
        sections["Frequencies (gnomAD)"] += [
            kv("gnomAD v4.1 joint AF",       per_gene.get("gnomad_v41_af_joint")),
            kv("gnomAD v4.1 grpmax FAF95",   per_gene.get("gnomad_v41_faf95_grpmax")),
            kv("gnomAD v4.1 grpmax ancestry", per_gene.get("gnomad_v41_grpmax_anc")),
            kv("gnomAD v4.1 nhomalt",        per_gene.get("gnomad_v41_nhomalt")),
            kv("gnomAD v4.1 FILTER",         per_gene.get("gnomad_v41_filter")),
            kv("gnomAD v2.1.1 controls AF",  per_gene.get("gnomAD2.1.1_exomes_controls_AF")),
            kv("1000G AF",                   per_gene.get("1000Gp3_AF")),
        ]
        sections["Predictors"] += [
            kv("REVEL",            per_gene.get("REVEL_score")),
            kv("CADD phred",       per_gene.get("CADD_phred")),
            kv("AlphaMissense",    per_gene.get("AlphaMissense_pred")),
            kv("AlphaMissense score", per_gene.get("AlphaMissense_score")),
            kv("SIFT pred",        per_gene.get("SIFT_pred")),
            kv("Polyphen2 HDIV",   per_gene.get("Polyphen2_HDIV_pred")),
            kv("MutationTaster",   per_gene.get("MutationTaster_pred")),
            kv("BayesDel",         per_gene.get("BayesDel_addAF_pred")),
            kv("PROVEAN",          per_gene.get("PROVEAN_pred")),
            kv("ClinPred",         per_gene.get("ClinPred_pred")),
            kv("MetaSVM",          per_gene.get("MetaSVM_pred")),
        ]
        sections["Splice & conservation"] += [
            kv("SpliceAI masked DS_MAX", per_gene.get("spliceai_ds_max_masked")),
            kv("PhastCons100way",        per_gene.get("phastCons100way_vertebrate")),
            kv("PhyloP100way",           per_gene.get("phyloP100way_vertebrate")),
            kv("PhyloP470way mammalian", per_gene.get("phyloP470way_mammalian")),
            kv("GERP++ RS",              per_gene.get("GERP++_RS")),
            kv("Interpro domain",        per_gene.get("Interpro_domain")),
        ]
        sections["Clinical assertions (ClinVar)"] += [
            kv("ClinVar ID",     per_gene.get("clinvar_id")),
            kv("ClinVar sig",    per_gene.get("clinvar_clnsig")),
            kv("ClinVar review", per_gene.get("clinvar_review")),
            kv("ClinVar trait",  per_gene.get("clinvar_trait")),
        ]
    else:
        sections["Notes"].append("- *Not in our cp_new dbNSFP cache.*")

    # --- LOCAL: pipeline classification ---
    classifications = find_in_classifications(v)
    if classifications:
        for c in classifications:
            sections["Our pipeline classification"].append(
                f"- **{c['classification']}** "
                f"(source: `{c.get('source','')}`, vv_status: `{c.get('vv_status','')}`)\n"
                f"    - transcript: {c.get('transcript','')}\n"
                f"    - c.: {c.get('hgvs_c','')}\n")
    else:
        sections["Our pipeline classification"].append(
            "- *No matching classification in `cp_new_seqnext_FINAL.tsv`.*")

    # --- EXTERNAL: myvariant.info (cross-source aggregator) ---
    if not args.no_external:
        mv = fetch_myvariant(v)
        if mv:
            sections["Cross-references"].append(kv("dbSNP rsID", mv.get("dbsnp", {}).get("rsid")))
            for src in ["cosmic", "clinvar", "dbsnp", "evs", "exac", "gnomad_exome", "gnomad_genome"]:
                if src in mv:
                    sections["Cross-references"].append(f"- **myvariant.{src}**: present (run `--raw-mv` to see full JSON)")
            mv_clinvar = mv.get("clinvar")
            if mv_clinvar:
                sections["Clinical assertions (ClinVar)"].append(
                    f"- myvariant.info ClinVar: {json.dumps(mv_clinvar)[:300]}...")

    # --- EXTERNAL: VariantValidator on the BED transcript ---
    if not args.no_external and v.gene and v.chr_hg19 and v.pos_hg19:
        # Try to figure out BED transcript from FINAL file
        bed_tx = next((c["transcript"] for c in classifications if c.get("transcript","").startswith("NM_")), None)
        if bed_tx:
            vv = fetch_variantvalidator(v, transcript=bed_tx)
            if vv:
                keys_with_c = [k for k in vv.keys() if ":c." in k]
                for k in keys_with_c[:3]:
                    sections["Identity"].append(f"- **VV resolved**: `{k}`")

    # --- Compose ---
    report = render_variant(v, sections)
    if args.json:
        print(json.dumps({"variant": asdict(v), "sections": sections}, default=str, indent=2))
    elif HAS_RICH and sys.stdout.isatty():
        Console().print(Markdown(report))
    else:
        print(report)
    return 0


# =============================================================================
# Subcommand: gene
# =============================================================================

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
    cols = ["gene","transcript","hgvs_c","classification","FAF95_grpmax",
            "REVEL","SpliceAI_masked","PhastCons100way","ClinVar_sig","source"]
    cols = [c for c in cols if c in sub.columns]
    if args.format == "tsv":
        sub[cols].to_csv(sys.stdout, sep="\t", index=False)
    else:
        # markdown table
        print(f"# {gene} — {len(sub)} classified variants\n")
        print("| " + " | ".join(cols) + " |")
        print("|" + "|".join(["---"] * len(cols)) + "|")
        for _, r in sub.iterrows():
            print("| " + " | ".join(str(r.get(c,"")).replace("|", "\\|") for c in cols) + " |")
    return 0


# =============================================================================
# Subcommand: bulk
# =============================================================================

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
        v = Variant(gene=r.get("gene",""), hgvs_c=r.get("hgvs_c",""))
        if "chr_grch37" in r and r["chr_grch37"]:
            v.chr_hg19 = r["chr_grch37"]
            v.pos_hg19 = int(r["pos_grch37"]) if r.get("pos_grch37") else None
            v.ref = r["ref"]
            v.alt = r["alt"]
        # Render minimal lookup (no external by default in bulk; opt-in)
        ns = argparse.Namespace(coord=None, rsid=None, gene=v.gene, hgvs=v.hgvs_c,
                                 transcript=None, build="hg19", json=False,
                                 no_external=not args.external)
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
            fname = f"{v.gene}_{(v.hgvs_c or '').replace('>','-').replace('/','_')[:40]}.md"
            (out_dir / fname).write_text(buf.getvalue())
        except SystemExit:
            continue
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(df)}", file=sys.stderr)
    print(f"done", file=sys.stderr)
    return 0


# =============================================================================
# Main / argparse
# =============================================================================

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("lookup", help="Single-variant report")
    pl.add_argument("--gene")
    pl.add_argument("--hgvs", help="c. notation, e.g. c.816C>T")
    pl.add_argument("--coord", help="chr:pos:ref:alt")
    pl.add_argument("--build", choices=["hg19","hg38"], default="hg19")
    pl.add_argument("--rsid")
    pl.add_argument("--transcript", help="Force a specific NM_ for HGVS resolution")
    pl.add_argument("--no-external", action="store_true", help="Skip API calls (local cache only)")
    pl.add_argument("--json", action="store_true", help="Emit JSON instead of Markdown")
    pl.set_defaults(func=cmd_lookup)

    pg = sub.add_parser("gene", help="Dump all classified variants for a gene")
    pg.add_argument("gene_symbol")
    pg.add_argument("--format", choices=["md","tsv"], default="md")
    pg.set_defaults(func=cmd_gene)

    pb = sub.add_parser("bulk", help="Process a TSV list of variants")
    pb.add_argument("input", help="TSV with at least gene+hgvs_c (or chr_grch37+pos_grch37+ref+alt)")
    pb.add_argument("--out", default="reports/", help="Output directory for per-variant Markdown files")
    pb.add_argument("--external", action="store_true", help="Include external API calls (slower)")
    pb.set_defaults(func=cmd_bulk)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

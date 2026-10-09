"""dbNSFP panel input and indexed VCF annotation helpers."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

GNOMAD_URL_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"
)
DBNSFP_GENE_FILE = "dbNSFP5.3_gene.gz"


def require_tabix() -> str:
    path = shutil.which("tabix")
    if not path:
        sys.exit("ERROR: tabix not found on PATH. brew install htslib  (or apt install tabix)")
    return path


def tabix_region(url_or_path: str, chrom: str, start: int, end: int) -> list[str]:
    """Run tabix on a region and return raw VCF lines."""
    region = f"{chrom}:{start}-{end}"
    try:
        result = subprocess.run(
            ["tabix", url_or_path, region],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
        )
    except subprocess.TimeoutExpired:
        print(f"  tabix TIMEOUT on {region} -> {url_or_path}", file=sys.stderr)
        return []
    if result.returncode != 0:
        print(
            f"  tabix error ({result.returncode}) {region}: {result.stderr.strip()[:200]}",
            file=sys.stderr,
        )
        return []
    return [line for line in result.stdout.splitlines() if line and not line.startswith("#")]


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def fetch_gnomad_region(chrom: str, start: int, end: int) -> pd.DataFrame:
    """Pull a region from gnomAD v4.1 joint sites and return key fields."""
    url = GNOMAD_URL_TEMPLATE.format(chrom=chrom)
    lines = tabix_region(url, f"chr{chrom}", start, end)
    rows = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        c, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
        # Multiallelic VCF rows are pre-split in gnomAD sites; ALT is single.
        ifields = parse_info(info)
        rows.append(
            {
                "gnomad_chr": c.removeprefix("chr"),
                "gnomad_pos": int(pos),
                "gnomad_ref": ref,
                "gnomad_alt": alt,
                # gnomAD v4.1 joint VCF uses _joint suffix; older releases dropped it
                "gnomad_v41_faf95_grpmax": ifields.get("fafmax_faf95_max_joint")
                or ifields.get("fafmax_faf95_max"),
                "gnomad_v41_af_joint": ifields.get("AF_joint") or ifields.get("AF"),
                "gnomad_v41_grpmax_anc": ifields.get("fafmax_faf95_max_gen_anc_joint")
                or ifields.get("grpmax_joint"),
                "gnomad_v41_nhomalt": ifields.get("nhomalt_joint") or ifields.get("nhomalt"),
                "gnomad_v41_filter": parts[6],  # PASS / EXOMES_FILTERED / etc.
            }
        )
    return pd.DataFrame(rows)


def fetch_spliceai_region(
    snv_vcf: Path | None,
    indel_vcf: Path | None,
    chrom: str,
    start: int,
    end: int,
) -> pd.DataFrame:
    """Pull a region from Illumina SpliceAI masked VCFs. Returns DS_MAX per site."""
    rows = []
    for src in (snv_vcf, indel_vcf):
        if src is None:
            continue
        lines = tabix_region(str(src), chrom, start, end)
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            c, pos, _id, ref, alt, _q, _f, info = parts[:8]
            # SpliceAI INFO: SpliceAI=ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            spliceai = parse_info(info).get("SpliceAI")
            if not spliceai:
                continue
            best = 0.0
            for entry in spliceai.split(","):
                fields = entry.split("|")
                if len(fields) < 6:
                    continue
                try:
                    ds = max(float(x) for x in fields[2:6] if x not in ("", "."))
                except ValueError:
                    continue
                best = max(best, ds)
            rows.append(
                {
                    "spliceai_chr": c.removeprefix("chr"),
                    "spliceai_pos": int(pos),
                    "spliceai_ref": ref,
                    "spliceai_alt": alt,
                    "spliceai_ds_max_masked": best,
                }
            )
    return pd.DataFrame(rows)


def load_gene_to_chr(data_dir: Path) -> dict[str, str]:
    """Read dbNSFP gene file and return {GENE_UPPER: chr_str}."""
    gf = data_dir / DBNSFP_GENE_FILE
    if not gf.exists():
        sys.exit(f"Gene file not found: {gf}")
    df = pd.read_csv(gf, sep="\t", usecols=["Gene_name", "chr"], dtype=str, na_values=".")
    df = df.dropna(subset=["Gene_name", "chr"])
    return dict(zip(df["Gene_name"].str.upper(), df["chr"].astype(str)))


def find_chr_file(chr_files: list[Path], chrom: str) -> Path | None:
    """Pick the dbNSFP chr file for a given chromosome label."""
    needle = f"chr{chrom}."
    for f in chr_files:
        if needle in f.name:
            return f
    return None


def read_dbnsfp_chr_for_genes(
    chr_file: Path,
    genes_upper: set[str],
    usecols: list[str],
) -> pd.DataFrame:
    """Single pass through one chr file, keep rows for any gene in genes_upper."""
    keep = []
    with gzip.open(chr_file, "rt") as f:
        reader = pd.read_csv(
            f,
            sep="\t",
            usecols=usecols,
            chunksize=50000,
            na_values=".",
            low_memory=False,
            dtype={
                "#chr": str,
                "ref": str,
                "alt": str,
                "aaref": str,
                "aaalt": str,
                "genename": str,
                "hg19_chr": str,
            },
        )
        for chunk in reader:
            # dbNSFP `genename` is sometimes a ';'-delimited multi-transcript
            # string (e.g. "POLE;POLE;POLE"). Match if ANY part is in the panel.
            up_raw = chunk["genename"].fillna("").str.upper()
            up_any = up_raw.apply(
                lambda g: next((p.strip() for p in g.split(";") if p.strip() in genes_upper), None)
            )
            mask = up_any.notna()
            if mask.any():
                hit = chunk[mask].copy()
                hit["_gene_upper"] = up_any[mask].values
                keep.append(hit)
    if not keep:
        return pd.DataFrame(columns=usecols + ["_gene_upper"])
    return pd.concat(keep, ignore_index=True)


def emit_consolidated_vcf(tsv_dir: Path, vcf_path: Path) -> None:
    """Write a minimal VCF of unique (chr, pos, ref, alt) across all per-gene TSVs.

    Output is for feeding to `spliceai -M 1` externally. Variants are written
    in chr,pos sort order. Chromosomes get a `chr` prefix to match hg38 builds.
    """
    seen: set[tuple[str, int, str, str]] = set()
    for tsv in sorted(tsv_dir.glob("*.tsv")):
        # Skip placeholder "# no variants" files; they have no header
        first = tsv.read_text().splitlines()[:1]
        if not first or not first[0].startswith("#chr"):
            continue
        try:
            df = pd.read_csv(
                tsv, sep="\t", usecols=["#chr", "pos(1-based)", "ref", "alt"], dtype=str
            )
        except (ValueError, pd.errors.EmptyDataError):
            continue
        for _, r in df.iterrows():
            try:
                pos = int(float(r["pos(1-based)"]))
            except (TypeError, ValueError):
                continue
            seen.add((str(r["#chr"]), pos, str(r["ref"]), str(r["alt"])))

    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    with open(vcf_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=.,Number=0,Type=Flag,Description="placeholder">\n')
        # contig lines for chr1..22, X, Y -- SpliceAI doesn't strictly need them
        # but downstream tools (bcftools) complain without them.
        out.writelines(
            f"##contig=<ID=chr{c}>\n" for c in [str(i) for i in range(1, 23)] + ["X", "Y"]
        )
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in sorted(seen, key=lambda x: (x[0], x[1])):
            chrom_out = chrom if chrom.startswith("chr") else f"chr{chrom}"
            out.write(f"{chrom_out}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"Emitted {len(seen):,} unique variants to {vcf_path}")

#!/usr/bin/env python3
"""Enumerate every possible synonymous SNV in cp_new BED CDS regions.

Solves dbNSFP's structural blind spot: dbNSFP only carries non-synonymous
+ splice-site variants. Synonymous coding variants (which the lab's W2 rule
specifically targets) are invisible to the dbNSFP-based pipeline. This
script generates an independent catalog by walking the BED-anchored RefSeq
NM_'s CDS exons against the hg19 reference FASTA and enumerating every
possible single-nucleotide substitution at each CDS position, keeping only
those that don't change the encoded amino acid.

Sources:
  - UCSC ncbiRefSeq hg19 table for CDS exon structure per NM_
  - UCSC hg19 reference FASTA for codon sequence
  - BioPython codon table for translation (handles selenocysteine/stop edge cases)

Output: data/exports/cp_new/cp_new_synonymous_catalog.tsv
  One row per (transcript, genomic position, alt allele) where alt is
  synonymous. Columns: gene, transcript, chr_grch37, pos_grch37, ref, alt,
  hgvs_c, codon_ref, codon_alt, aa, strand.

Downstream: annotate this catalog with gnomAD FAF, PhastCons, PhyloP,
SpliceAI; then apply the W2 classifier rule. See annotate_synonymous_catalog.py.
"""

from __future__ import annotations

import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam
from Bio.Seq import Seq
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.panels import get_panel

REFSEQ_TXT = Path.home() / "data" / "refseq_hg19" / "ncbiRefSeq.txt.gz"
HG19_FA    = Path.home() / "data" / "hg19" / "hg19.fa"
# Prefer the lab's source, fall back to the snapshot we keep in the bundle.
_BED_CANDIDATES = [
    Path("/Users/nuin/Projects/ahs/new_bed/CP_new/C+_ALL_IDPE_APR2026.bed"),
    Path("cp_new_bundle/bed/C+_ALL_IDPE_APR2026.bed"),
]
BED = next((p for p in _BED_CANDIDATES if p.exists()), _BED_CANDIDATES[0])
OUT        = Path("data/exports/cp_new/cp_new_synonymous_catalog.tsv")

COMP = str.maketrans("ACGT", "TGCA")


def load_bed_transcripts(bed: Path, panel: set[str]) -> dict[str, str]:
    """Return {gene: most-common NM_} from BED for the panel."""
    counts: dict[str, Counter] = defaultdict(Counter)
    with open(bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6: continue
            if p[4] not in panel: continue
            counts[p[4]][p[5]] += 1
    return {g: c.most_common(1)[0][0] for g, c in counts.items()}


def parse_refseq(refseq_path: Path, wanted_nm: set[str]) -> dict[str, dict]:
    """For each wanted NM_, return CDS structure.

    Match by NM_ stem (NM_002080) -- versions (.4 .5) may differ slightly
    between BED and ncbiRefSeq, but the CDS exon coordinates are stable for
    most transcripts. Prefer exact version match, fall back to stem match.
    """
    wanted_stems = {n.split(".")[0]: n for n in wanted_nm}
    by_exact: dict[str, dict] = {}
    by_stem: dict[str, dict] = {}
    with gzip.open(refseq_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 16: continue
            nm = fields[1]
            stem = nm.split(".")[0]
            if nm not in wanted_nm and stem not in wanted_stems:
                continue
            try:
                row = {
                    "name": nm,
                    "chrom": fields[2],
                    "strand": fields[3],
                    "cds_start": int(fields[6]),  # 0-based half-open
                    "cds_end": int(fields[7]),
                    "exon_starts": [int(x) for x in fields[9].rstrip(",").split(",")],
                    "exon_ends":   [int(x) for x in fields[10].rstrip(",").split(",")],
                    "gene": fields[12],
                }
            except (ValueError, IndexError):
                continue
            if row["cds_end"] <= row["cds_start"]:
                continue  # non-coding
            # Only consider main chr (skip alt contigs, hap, etc.)
            if "_" in row["chrom"]:
                continue
            if nm in wanted_nm:
                by_exact[nm] = row
            if stem in wanted_stems:
                # Keep the longest transcript per stem if multiple
                prev = by_stem.get(stem)
                if not prev or (row["cds_end"]-row["cds_start"]) > (prev["cds_end"]-prev["cds_start"]):
                    by_stem[stem] = row

    out = {}
    for nm in wanted_nm:
        if nm in by_exact:
            out[nm] = by_exact[nm]
        else:
            stem = nm.split(".")[0]
            if stem in by_stem:
                out[nm] = by_stem[stem]
    return out


def build_cds_seq(tx: dict, fa: pysam.FastaFile) -> tuple[str, list[int]]:
    """Return (cds_seq, [genomic_positions_0based]) on the coding strand."""
    positions: list[int] = []
    for s, e in zip(tx["exon_starts"], tx["exon_ends"]):
        cs = max(s, tx["cds_start"])
        ce = min(e, tx["cds_end"])
        if cs >= ce: continue
        positions.extend(range(cs, ce))
    if not positions:
        return "", []
    chrom = tx["chrom"]
    # Fetch in one call (much faster than per-base)
    spans = []
    run_s = positions[0]
    run_p = positions[0]
    for p in positions[1:]:
        if p == run_p + 1:
            run_p = p
        else:
            spans.append((run_s, run_p + 1))
            run_s = p
            run_p = p
    spans.append((run_s, run_p + 1))
    seq = "".join(fa.fetch(chrom, s, e).upper() for s, e in spans)
    if tx["strand"] == "-":
        seq = str(Seq(seq).reverse_complement())
        positions = list(reversed(positions))
    return seq, positions


def enumerate_synonymous(gene: str, nm_bed: str, tx: dict, fa: pysam.FastaFile):
    """Yield dicts for each synonymous SNV in this transcript's CDS."""
    cds, positions = build_cds_seq(tx, fa)
    if not cds:
        return
    strand = tx["strand"]
    chrom = tx["chrom"].removeprefix("chr")
    # Walk full codons only
    last_codon_start = (len(cds) // 3) * 3
    for i in range(last_codon_start):
        codon_start = (i // 3) * 3
        codon = cds[codon_start:codon_start + 3]
        if "N" in codon: continue
        ref_aa = str(Seq(codon).translate())
        cds_ref = cds[i]
        if cds_ref not in "ACGT": continue
        codon_pos = i - codon_start  # 0,1,2
        # Genomic ref: + strand same as CDS, - strand complement
        genomic_pos_1b = positions[i] + 1
        genomic_ref = cds_ref if strand == "+" else cds_ref.translate(COMP)
        for genomic_alt in "ACGT":
            if genomic_alt == genomic_ref: continue
            cds_alt = genomic_alt if strand == "+" else genomic_alt.translate(COMP)
            new_codon = codon[:codon_pos] + cds_alt + codon[codon_pos + 1:]
            alt_aa = str(Seq(new_codon).translate())
            if alt_aa != ref_aa:
                continue  # not synonymous
            c_pos = i + 1  # 1-based
            yield {
                "gene": gene,
                "transcript": nm_bed,
                "chr_grch37": chrom,
                "pos_grch37": genomic_pos_1b,
                "ref": genomic_ref,
                "alt": genomic_alt,
                "hgvs_c": f"c.{c_pos}{cds_ref}>{cds_alt}",
                "codon_ref": codon,
                "codon_alt": new_codon,
                "aa": ref_aa,
                "strand": strand,
                "consequence": "synonymous_variant",
            }


def main() -> int:
    if not HG19_FA.exists():
        sys.exit(f"hg19 FASTA not found: {HG19_FA}")
    if not REFSEQ_TXT.exists():
        sys.exit(f"RefSeq table not found: {REFSEQ_TXT}")

    panel = set(get_panel("cp_new"))
    gene_tx = load_bed_transcripts(BED, panel)
    wanted_nm = set(gene_tx.values())
    print(f"cp_new genes: {len(panel)}, with BED transcript: {len(gene_tx)}")

    refseq = parse_refseq(REFSEQ_TXT, wanted_nm)
    print(f"Resolved {len(refseq)} of {len(wanted_nm)} NM_ in ncbiRefSeq hg19")

    missing = [nm for nm in wanted_nm if nm not in refseq]
    if missing:
        print(f"  missing NM_ (likely non-coding or hg38-only versions): {missing[:10]}")

    fa = pysam.FastaFile(str(HG19_FA))
    OUT.parent.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    n_genes = 0
    with open(OUT, "w") as out:
        header = ["gene","transcript","chr_grch37","pos_grch37","ref","alt",
                  "hgvs_c","codon_ref","codon_alt","aa","strand","consequence"]
        out.write("\t".join(header) + "\n")
        for gene, nm in tqdm(sorted(gene_tx.items()), desc="enumerate"):
            tx = refseq.get(nm)
            if tx is None: continue
            gene_n = 0
            for row in enumerate_synonymous(gene, nm, tx, fa):
                out.write("\t".join(str(row[k]) for k in header) + "\n")
                gene_n += 1
                n_rows += 1
            if gene_n: n_genes += 1

    print(f"\nWrote {n_rows:,} synonymous candidates across {n_genes} genes -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

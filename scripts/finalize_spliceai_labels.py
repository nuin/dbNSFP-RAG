#!/usr/bin/env python3
"""Finalize SpliceAI honesty in the upload (lab request 2026-07):

1. Backfill any missing missense masked/unmasked from ALL SpliceAI VCFs.
2. Fix STALE labels: rows with a real measured SpliceAI whose classification
   text still says '(assumed 0)' -> '(measured)'. (The score was measured; only
   the label was stale, which made the lab think it was assumed.)
3. EXCLUDE rows that are GENUINELY assumed — SpliceAI_masked == '0 (assumed)',
   i.e. SpliceAI could not score them (ref mismatch / contig edge). Per the lab:
   "if these scores are actually just assumed and not based on SpliceAI data, we
   would want these rows excluded."

Updates v3_3 / v3_4 (full + roi), regenerates v3_5, refreshes per-gene.
"""
from __future__ import annotations
import glob
from pathlib import Path
import pandas as pd
from pyliftover import LiftOver

VCF_DIR = Path("data/exports/cp_new")
MASKED_VCFS = ["cp_new.spliceai.vcf", "cp_new.bed.chunk_*.spliceai.vcf",
               "cp_new_v2_needed.spliceai.vcf", "missense_spliceai.masked.vcf",
               "syn_assumed0_spliceai.masked.vcf", "missense_gap_spliceai.masked.vcf"]
UNMASKED_VCFS = ["cp_new_v2_needed.spliceai_unmasked.vcf",
                 "missense_spliceai.unmasked.vcf",
                 "syn_assumed0_spliceai.unmasked.vcf", "missense_gap_spliceai.unmasked.vcf"]
TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]


def load(patterns):
    m = {}
    for pat in patterns:
        for path in glob.glob(str(VCF_DIR / pat)):
            for line in open(path):
                if line.startswith("#"): continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 8: continue
                c, pos, _i, ref, alt = p[0].removeprefix("chr"), p[1], p[2], p[3], p[4]
                sa = next((kv[9:] for kv in p[7].split(";") if kv.startswith("SpliceAI=")), None)
                if not sa: continue
                best = 0.0
                for e in sa.split(","):
                    f = e.split("|")
                    if len(f) < 6: continue
                    try: best = max(best, max(float(x) for x in f[2:6] if x not in ("", ".")))
                    except ValueError: pass
                m.setdefault((c, int(pos), ref, alt), best)
    return m


def main():
    masked, unmasked = load(MASKED_VCFS), load(UNMASKED_VCFS)
    print(f"SpliceAI VCF index: masked {len(masked):,}, unmasked {len(unmasked):,}")
    lo = LiftOver("hg19", "hg38"); _c = {}
    def h38(c, p):
        k = (c, p)
        if k in _c: return _c[k]
        r = lo.convert_coordinate(f"chr{c}", p - 1)
        v = r[0][1] + 1 if r else None
        _c[k] = v; return v

    for path, pgname, tag in TABLES:
        p = Path(path)
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        n_fill_u = n_fill_m = n_relabel = 0
        for i in df.index:
            r = df.loc[i]
            if not r["chr_grch37"] or not r["pos_grch37"]: continue
            pos38 = h38(r["chr_grch37"], int(float(r["pos_grch37"])))
            k = (r["chr_grch37"], pos38, r["ref"], r["alt"]) if pos38 else None
            # backfill any missing masked/unmasked (esp. missense unmasked gap)
            if k and r["SpliceAI_masked"] in ("", "0 (assumed)") and masked.get(k) is not None:
                df.at[i, "SpliceAI_masked"] = f"{masked[k]:.4f}"; n_fill_m += 1
            if k and r["SpliceAI_unmasked"] == "" and unmasked.get(k) is not None:
                df.at[i, "SpliceAI_unmasked"] = f"{unmasked[k]:.4f}"; n_fill_u += 1
            # fix stale '(assumed 0)' label where SpliceAI is actually measured
            r = df.loc[i]
            if df.at[i, "SpliceAI_masked"] not in ("", "0 (assumed)") and "(assumed 0)" in r["classification"]:
                df.at[i, "classification"] = r["classification"].replace("(assumed 0)", "(measured)")
                n_relabel += 1

        # EXCLUDE genuinely-assumed rows (SpliceAI could not score them)
        before = len(df)
        excl = df["SpliceAI_masked"].eq("0 (assumed)")
        excluded = df[excl]
        df = df[~excl]
        df.to_csv(p, sep="\t", index=False)
        excluded.to_csv(p.with_name(p.stem + "_assumed_excluded.tsv"), sep="\t", index=False)
        print(f"  {p.name}: filled masked {n_fill_m}, unmasked {n_fill_u} | "
              f"relabeled stale {n_relabel} | excluded genuinely-assumed {before-len(df)} -> {len(df):,}")

        pg = p.parent / pgname; pg.mkdir(exist_ok=True)
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{tag}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))

    # regenerate v3.5
    src = "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv"
    v5 = Path("cp_new_bundle/outputs/new_genes_68_v3_5")
    d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    d.to_csv(v5/"UPLOAD_v3_5.tsv", sep="\t", index=False)
    allg = sorted(set(open("cp_new_bundle/outputs/_shared/new_genes.txt").read().split()))
    hdr = "\t".join(d.columns)
    for f in (v5/"per_gene").glob("*.tsv"): f.unlink()
    for g in allg:
        sub = d[d["gene"]==g]
        (v5/"per_gene"/f"{g}__v3_5.tsv").write_text(sub.to_csv(sep="\t",index=False) if len(sub) else hdr+"\n")
    print(f"\nv3.5: {len(d):,} rows")
    print(f"  still 'assumed': {d['SpliceAI_masked'].eq('0 (assumed)').sum()}")
    print(f"  label says 'assumed': {d['classification'].str.contains('assumed',na=False).sum()}")
    mis = d[d["classification"].str.startswith("Likely_benign")]
    print(f"  missense missing unmasked: {(mis['SpliceAI_unmasked']=='').sum()} / {len(mis)}")


if __name__ == "__main__":
    main()

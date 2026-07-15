#!/usr/bin/env python3
"""Backfill real SpliceAI (masked + unmasked) onto the synonymous rows that had
'0 (assumed)' placeholders, and apply the agreed W2 splice gate (SpliceAI <= 0.1
on the unmasked/spliceai.org value) with the real numbers.

Rows whose real SpliceAI > 0.1 were WRONGLY eligible for benign (e.g. HOXB13
c.600A>G, real 0.27). They are relabeled to a splice-signal review flag and kept
in the list (per the keep+filter preference), with real scores populated.

Updates v3_3 / v3_4 (full + roi) in place; regenerates v3_5. Per-gene refreshed.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from pyliftover import LiftOver

MASKED = Path("data/exports/cp_new/syn_assumed0_spliceai.masked.vcf")
UNMASKED = Path("data/exports/cp_new/syn_assumed0_spliceai.unmasked.vcf")
SPLICE_MAX = 0.1

TABLES = [
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",     "per_gene",     "v3_3"),
    ("cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv", "per_gene_roi", "v3_3_roi"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",     "per_gene",     "v3_4"),
    ("cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv", "per_gene_roi", "v3_4_roi"),
]


def parse_vcf(path):
    out = {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 8: continue
        c, pos, _i, ref, alt = p[0].removeprefix("chr"), p[1], p[2], p[3], p[4]
        sa = None
        for kv in p[7].split(";"):
            if kv.startswith("SpliceAI="): sa = kv[9:]; break
        if not sa: continue
        best = 0.0
        for e in sa.split(","):
            f = e.split("|")
            if len(f) < 6: continue
            try: best = max(best, max(float(x) for x in f[2:6] if x not in ("", ".")))
            except ValueError: pass
        out[(c, int(pos), ref, alt)] = best
    return out


def main():
    masked = parse_vcf(MASKED)
    unmasked = parse_vcf(UNMASKED)
    print(f"parsed masked {len(masked):,}, unmasked {len(unmasked):,}")
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
        assumed = df["SpliceAI_masked"].eq("0 (assumed)")
        n_ok, n_fail, n_nomatch = 0, 0, 0
        for i in df[assumed].index:
            r = df.loc[i]
            if not r["chr_grch37"] or not r["pos_grch37"]: n_nomatch += 1; continue
            pos38 = h38(r["chr_grch37"], int(float(r["pos_grch37"])))
            if not pos38: n_nomatch += 1; continue
            k = (r["chr_grch37"], pos38, r["ref"], r["alt"])
            m, u = masked.get(k), unmasked.get(k)
            if m is None and u is None: n_nomatch += 1; continue
            if m is not None: df.at[i, "SpliceAI_masked"] = f"{m:.4f}"
            df.at[i, "SpliceAI_unmasked"] = f"{u:.4f}" if u is not None else ""
            gate = u if u is not None else m   # spliceai.org = unmasked
            cls = r["classification"]
            if gate is not None and gate > SPLICE_MAX:
                # real splice signal -> not benign (agreed W2 gate, real data)
                df.at[i, "classification"] = (f"Synonymous, SpliceAI={gate:.2f} >0.1 "
                                              f"(splice signal - review, not benign)")
                n_fail += 1
            else:
                # real value <= 0.1 -> drop the '(assumed 0)' wording, keep call
                df.at[i, "classification"] = cls.replace("SpliceAI<=0.1 (assumed 0)",
                                                         "SpliceAI<=0.1 (measured)")
                n_ok += 1
        df.to_csv(p, sep="\t", index=False)
        print(f"  {p.name}: measured-ok {n_ok:,} | splice-fail relabeled {n_fail:,} | "
              f"no-match {n_nomatch:,}")
        # per-gene refresh
        pg = p.parent / pgname; pg.mkdir(exist_ok=True)
        for old in pg.glob("*.tsv"): old.unlink()
        for g in sorted(df["gene"].unique()):
            (pg/f"{g}__{tag}.tsv").write_text(df[df["gene"]==g].to_csv(sep="\t", index=False))

    # regenerate v3_5 from v3_4_roi
    src = "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv"
    v5dir = Path("cp_new_bundle/outputs/new_genes_68_v3_5")
    d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    d.to_csv(v5dir/"UPLOAD_v3_5.tsv", sep="\t", index=False)
    allg = sorted(set(open("cp_new_bundle/outputs/_shared/new_genes.txt").read().split()))
    hdr = "\t".join(d.columns)
    (v5dir/"per_gene").mkdir(exist_ok=True)
    for f in (v5dir/"per_gene").glob("*.tsv"): f.unlink()
    for g in allg:
        sub = d[d["gene"]==g]
        (v5dir/"per_gene"/f"{g}__v3_5.tsv").write_text(sub.to_csv(sep="\t",index=False) if len(sub) else hdr+"\n")
    print(f"\nv3_5 regenerated: {len(d):,} rows")

    # HOXB13 c.600A>G verification
    h = d[(d.gene=="HOXB13") & (d.hgvs_c=="c.600A>G")]
    if len(h):
        r = h.iloc[0]
        print(f"\nHOXB13 c.600A>G: masked={r['SpliceAI_masked']} unmasked={r['SpliceAI_unmasked']}")
        print(f"  -> {r['classification']}")


if __name__ == "__main__":
    main()

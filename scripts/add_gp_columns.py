#!/usr/bin/env python3
"""Add g. (hg19 + hg38), gnomAD variant ID, and p. columns to the v3.3/v3.4
uploads (full + ROI).

We are an hg19 shop, so the native g. is hg19 (our coordinates, no liftover).
For locating variants in gnomAD v4.1 (hg38) we ALSO provide an hg38 g. + gnomAD
ID produced by liftover -- uniformly for ALL genes (one consistent method, not
gnomAD's own VEP g. for some genes and liftover for others).

Columns added after `alt`:
  g_hg19        e.g. chr5:g.79950724G>C           (native, our build)
  gnomad_id     e.g. 5-80654905-G-C               (hg38, paste into gnomAD v4.1)
  g_hg38        e.g. chr5:g.80654905G>C           (hg38, via liftover)
  p_hgvs        e.g. p.Ala60Pro / p.Gln664=       (from gnomAD annotation; blank if none)
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

GP = Path("data/exports/cp_new/gnomad_hgvsc/_gp_map.tsv")
RESOLVED = Path("cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_resolved.tsv")
TABLES = [
    "cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_3/UPLOAD_v3_3_roi.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4.tsv",
    "cp_new_bundle/outputs/new_genes_68_v3_4/UPLOAD_v3_4_roi.tsv",
]


def nm_base(nm): return nm.split(".")[0] if nm else nm


def g_string(chrom, pos, ref, alt):
    """gnomAD-style g. Simple SNV -> chr:g.POSref>alt; indel shown with del/ins."""
    if not chrom or not pos:
        return ""
    if len(ref) == 1 and len(alt) == 1:
        return f"chr{chrom}:g.{pos}{ref}>{alt}"
    # deletion
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        s = int(pos) + 1; e = int(pos) + len(ref) - 1
        return f"chr{chrom}:g.{s}_{e}del" if e > s else f"chr{chrom}:g.{s}del"
    # insertion
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return f"chr{chrom}:g.{pos}_{int(pos)+1}ins{alt[1:]}"
    # delins / mnv
    e = int(pos) + len(ref) - 1
    return f"chr{chrom}:g.{pos}_{e}delins{alt}"


def main():
    # p. map (build-independent); pull may still be finishing -- ok, blanks fill
    lut = {}
    if GP.exists():
        gp = pd.read_csv(GP, sep="\t", dtype=str).fillna("")
        lut = {(r["gene"], nm_base(r["nm"]), r["cdot"]): r["p_hgvs"] for _, r in gp.iterrows()}
    print(f"p. map: {len(lut):,} keys")

    from pyliftover import LiftOver
    lo = LiftOver("hg19", "hg38"); _c = {}
    def lift(chrom, pos1):
        k = (chrom, pos1)
        if k in _c: return _c[k]
        res = lo.convert_coordinate(f"chr{chrom}", pos1 - 1)
        out = (res[0][1] + 1) if res else None
        _c[k] = out; return out
    lo_rev = LiftOver("hg38", "hg19"); _cr = {}
    def lift_rev(chrom, pos1):
        k = (chrom, pos1)
        if k in _cr: return _cr[k]
        res = lo_rev.convert_coordinate(f"chr{chrom}", pos1 - 1)
        out = (res[0][1] + 1) if res else None
        _cr[k] = out; return out

    # gnomad_common rows carry hg38 coords in the resolved file (empty hg19 in upload)
    hg38_map = {}       # (gene, transcript, hgvs_c) -> (chr38, pos38)
    # Known-coordinate overrides for gnomad_common rows not in the resolved file
    # (transcript-mismatch edge cases). hg38 verified via MyVariant/gnomAD.
    hg38_map_gc = {     # (gene, hgvs_c) -> (chr38, pos38)
        ("DPYD", "c.85T>C"): ("1", "97883329"),  # rs, DPYD*9A; A>G minus strand
    }
    if RESOLVED.exists():
        rv = pd.read_csv(RESOLVED, sep="\t", dtype=str).fillna("")
        for _, r in rv.iterrows():
            hg38_map[(r["gene"], r["transcript"], r["hgvs_c"])] = (r["chr_hg38"], r["pos_hg38"])
            hg38_map_gc.setdefault((r["gene"], r["hgvs_c"]), (r["chr_hg38"], r["pos_hg38"]))

    for path in TABLES:
        p = Path(path)
        if not p.exists():
            print(f"  SKIP (missing): {path}"); continue
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        # strip any prior g./p. columns so this is idempotent
        for c in ("gnomad_id", "g_hg38", "g_hg19", "p_hgvs"):
            if c in df.columns: df = df.drop(columns=[c])
        g19s, g38s, gids, phs = [], [], [], []
        n_lift = n_nolift = 0
        # UNIFORM METHOD (hg19 shop): every row is anchored on hg19, and the
        # hg38 columns are ALWAYS forward-liftover(hg19) -- one method, all genes.
        # hg19-native rows use their coords directly; gnomad_common rows (which
        # arrived as hg38) are first reverse-lifted to hg19 so they, too, are
        # anchored on hg19 and then forward-lifted like everyone else.
        for _, r in df.iterrows():
            chrom, pos19, ref, alt = r["chr_grch37"], r["pos_grch37"], r["ref"], r["alt"]
            phs.append(lut.get((r["gene"], nm_base(r["transcript"]), r["hgvs_c"]), ""))

            # Establish the hg19 anchor coordinate for this row
            if chrom and pos19:
                hg19_chr, hg19_pos = chrom, int(float(pos19))
            else:
                # gnomad_common: reverse-lift its gnomAD hg38 coord to hg19
                h38 = (hg38_map.get((r["gene"], r["transcript"], r["hgvs_c"]))
                       or hg38_map_gc.get((r["gene"], r["hgvs_c"])))
                hg19_chr = hg19_pos = None
                if h38 and h38[0] and h38[1]:
                    try:
                        h19 = lift_rev(h38[0], int(float(h38[1])))
                    except Exception:
                        h19 = None
                    if h19:
                        hg19_chr, hg19_pos = h38[0], h19

            if hg19_chr is None or hg19_pos is None:
                g19s.append(""); g38s.append(""); gids.append(""); n_nolift += 1
                continue

            g19s.append(g_string(hg19_chr, str(hg19_pos), ref, alt))
            # hg38 ALWAYS via forward liftover of the hg19 anchor
            try:
                hp = lift(hg19_chr, hg19_pos)
            except Exception:
                hp = None
            if hp:
                g38s.append(g_string(hg19_chr, str(hp), ref, alt))
                gids.append(f"{hg19_chr}-{hp}-{ref}-{alt}")
                n_lift += 1
            else:
                g38s.append(""); gids.append(""); n_nolift += 1
        ins = list(df.columns).index("alt") + 1
        df.insert(ins, "g_hg19", g19s)
        df.insert(ins + 1, "gnomad_id", gids)
        df.insert(ins + 2, "g_hg38", g38s)
        df.insert(ins + 3, "p_hgvs", phs)
        df.to_csv(p, sep="\t", index=False)
        pcov = sum(1 for x in phs if x)
        print(f"  {p.name}: {len(df):,} rows | lifted {n_lift:,}, no-lift {n_nolift:,} | p. filled {pcov:,}")

    for folder, stem, tag in [
        ("cp_new_bundle/outputs/new_genes_68_v3_3","UPLOAD_v3_3","v3_3"),
        ("cp_new_bundle/outputs/new_genes_68_v3_4","UPLOAD_v3_4","v3_4"),
    ]:
        for src, pgdir, suffix in [(f"{folder}/{stem}.tsv", f"{folder}/per_gene", tag),
                                   (f"{folder}/{stem}_roi.tsv", f"{folder}/per_gene_roi", f"{tag}_roi")]:
            if not Path(src).exists(): continue
            d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
            pg = Path(pgdir); pg.mkdir(exist_ok=True)
            for old in pg.glob("*.tsv"): old.unlink()
            for g in sorted(d["gene"].unique()):
                (pg/f"{g}__{suffix}.tsv").write_text(d[d["gene"]==g].to_csv(sep="\t", index=False))
    print("per-gene folders refreshed")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compare FAF from gnomAD v2.1.1 (hg19 native) against the FAF we pulled
from gnomAD v4.1 (hg38 -> dbNSFP hg19-lift). Measures how much liftover
variance the dbNSFP-bundled hg19 coords introduce vs querying gnomAD
natively in hg19.

For each cp_new gene:
  1. Get the hg19 region from per-gene TSV (min/max hg19_pos).
  2. tabix gnomAD v2.1.1 exomes sites VCF for that region (hg19, no 'chr' prefix).
  3. Extract AF_popmax + per-population faf95_{afr,amr,eas,nfe,sas}, take max.
  4. Join against per-gene TSV on (hg19_chr, hg19_pos, ref, alt).
  5. Compare against gnomad_v41_faf95_grpmax + gnomad_v41_af_joint.

Writes:
  data/exports/cp_new/cp_new_gnomad_v2_compare.tsv      -- per-variant comparison
  data/exports/cp_new/cp_new_gnomad_v2_summary.txt      -- aggregate stats

Caveats:
  - v2.1.1 (~125k exomes) is genuinely smaller than v4.1 (~807k joint), so
    a fraction of variants present in v4 won't be in v2 at all. That's
    biology + sampling, not liftover.
  - Compare percentage agreement on the variants present in BOTH releases.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

DATA = Path("data/exports/cp_new")
GNOMAD_V2_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"

POPS = ["afr", "amr", "eas", "nfe", "sas"]  # exclude fin/asj/oth per ACMG convention


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def safe_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def fetch_v2_region(chrom: str, start: int, end: int) -> pd.DataFrame:
    """Pull gnomAD v2.1.1 exomes sites VCF for a region, hg19 coords."""
    url = GNOMAD_V2_URL.format(chrom=chrom)
    try:
        r = subprocess.run(
            ["tabix", url, f"{chrom}:{start}-{end}"],
            capture_output=True, text=True, timeout=600, check=False,
        )
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {chrom}:{start}-{end}", file=sys.stderr)
        return pd.DataFrame()
    if r.returncode != 0:
        return pd.DataFrame()
    rows = []
    for line in r.stdout.splitlines():
        if line.startswith("#") or not line:
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, _id, ref, alt, _q, filt, info = parts[:8]
        ifields = parse_info(info)
        af_popmax = safe_float(ifields.get("AF_popmax"))
        # Per-population FAF95, take max excluding fin/asj/oth
        faf_per_pop = [safe_float(ifields.get(f"faf95_{p}")) for p in POPS]
        faf_per_pop = [v for v in faf_per_pop if v is not None]
        faf_max = max(faf_per_pop) if faf_per_pop else None
        rows.append({
            "v2_chr": chrom, "v2_pos": int(pos), "v2_ref": ref, "v2_alt": alt,
            "v2_af_popmax": af_popmax,
            "v2_faf95_popmax": faf_max,
            "v2_filter": filt,
        })
    return pd.DataFrame(rows)


def main() -> int:
    tsvs = sorted(t for t in DATA.glob("*.tsv")
                  if t.stem not in ("cp_new_seqnext_minimal", "cp_new_seqnext_strict"))
    print(f"Comparing {len(tsvs)} cp_new genes against gnomAD v2.1.1 hg19 native\n")

    all_chunks = []
    for tsv in tqdm(tsvs, desc="genes"):
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_values=[".",""], low_memory=False)
        if df.empty or "hg19_chr" not in df.columns:
            continue
        df["hg19_pos"] = pd.to_numeric(df["hg19_pos(1-based)"], errors="coerce").astype("Int64")
        df = df.dropna(subset=["hg19_pos","hg19_chr"])
        if df.empty:
            continue
        chrom = str(df["hg19_chr"].iloc[0])
        pos_min = int(df["hg19_pos"].min())
        pos_max = int(df["hg19_pos"].max())
        v2 = fetch_v2_region(chrom, pos_min, pos_max)
        if v2.empty:
            continue
        df["hg19_chr"] = df["hg19_chr"].astype(str)
        merged = df.merge(
            v2,
            left_on=["hg19_chr","hg19_pos","ref","alt"],
            right_on=["v2_chr","v2_pos","v2_ref","v2_alt"],
            how="left",
        )
        # Only keep rows that joined OR had gnomAD v4 data
        merged["v4_faf"] = pd.to_numeric(merged.get("gnomad_v41_faf95_grpmax"), errors="coerce")
        merged["v4_af"] = pd.to_numeric(merged.get("gnomad_v41_af_joint"), errors="coerce")
        # Compute disagreement category
        v2_present = merged["v2_af_popmax"].notna()
        v4_present = merged["v4_af"].notna()
        merged["status"] = "no_data"
        merged.loc[v2_present & ~v4_present, "status"] = "v2_only"
        merged.loc[v4_present & ~v2_present, "status"] = "v4_only"
        merged.loc[v2_present & v4_present, "status"] = "both"

        keep = ["genename","hg19_chr","hg19_pos","ref","alt",
                "v4_faf","v4_af","v2_af_popmax","v2_faf95_popmax","v2_filter","status"]
        all_chunks.append(merged[[c for c in keep if c in merged.columns]])

    out = pd.concat(all_chunks, ignore_index=True)
    out_path = DATA / "cp_new_gnomad_v2_compare.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {len(out):,} rows -> {out_path}")

    # === aggregate stats ===
    summary_path = DATA / "cp_new_gnomad_v2_summary.txt"
    lines = []
    def L(s=""):
        lines.append(s); print(s)
    L("=" * 70)
    L("gnomAD v4.1 hg38 (via dbNSFP-lift) vs v2.1.1 hg19 native")
    L("=" * 70)
    L(f"Total variants in cp_new genes:    {len(out):,}")
    L(f"  In both releases:                {(out['status']=='both').sum():,}")
    L(f"  v4.1 only (newer / hg38 native): {(out['status']=='v4_only').sum():,}")
    L(f"  v2.1.1 only (rare):              {(out['status']=='v2_only').sum():,}")
    L(f"  In neither (rare/novel):         {(out['status']=='no_data').sum():,}")
    L()

    both = out[out["status"] == "both"].copy()
    # Compare AF: v2 popmax vs v4 popmax-equivalent (af_joint is total, not popmax;
    # closest direct comparison is v2 AF_popmax vs v4 grpmax FAF -- not identical
    # metrics but the closest available across the two)
    both["v2_af_p"] = pd.to_numeric(both["v2_af_popmax"], errors="coerce")
    both["v4_faf_p"] = pd.to_numeric(both["v4_faf"], errors="coerce")
    pair = both.dropna(subset=["v2_af_p","v4_faf_p"]).copy()
    pair["delta"] = (pair["v2_af_p"] - pair["v4_faf_p"]).abs()
    pair["max_af"] = pair[["v2_af_p","v4_faf_p"]].max(axis=1)
    pair["rel_delta"] = pair["delta"] / pair["max_af"].clip(lower=1e-6)

    L(f"Variants with FAF/AF in BOTH releases: {len(pair):,}")
    L(f"  Exact agreement (delta <= 1e-7):     {(pair['delta'] <= 1e-7).sum():,}")
    L(f"  Within 10% relative:                 {(pair['rel_delta'] <= 0.10).sum():,}")
    L(f"  Within 50% relative:                 {(pair['rel_delta'] <= 0.50).sum():,}")
    L(f"  >2x apart:                           {(pair['rel_delta'] > 1.0).sum():,}")
    L()

    # W1/W3 rule fires: compare which release would have fired the rule
    L("Rule-fire comparison (W1 = FAF > 5%, W3 needs FAF > 0.1%):")
    L(f"  Both >5%:                            {((pair['v4_faf_p'] > 0.05) & (pair['v2_af_p'] > 0.05)).sum():,}")
    L(f"  v4 >5% but v2 not:                   {((pair['v4_faf_p'] > 0.05) & (pair['v2_af_p'] <= 0.05)).sum():,}")
    L(f"  v2 >5% but v4 not:                   {((pair['v2_af_p'] > 0.05) & (pair['v4_faf_p'] <= 0.05)).sum():,}")
    L(f"  Both >0.1%:                          {((pair['v4_faf_p'] > 0.001) & (pair['v2_af_p'] > 0.001)).sum():,}")
    L(f"  v4 >0.1% but v2 not:                 {((pair['v4_faf_p'] > 0.001) & (pair['v2_af_p'] <= 0.001)).sum():,}")
    L(f"  v2 >0.1% but v4 not:                 {((pair['v2_af_p'] > 0.001) & (pair['v4_faf_p'] <= 0.001)).sum():,}")
    L()

    L("=== Worst disagreements (top 15 by absolute FAF delta) ===")
    worst = pair.nlargest(15, "delta")[
        ["genename","hg19_chr","hg19_pos","ref","alt","v4_faf_p","v2_af_p","delta"]]
    L(worst.to_string(index=False))

    summary_path.write_text("\n".join(lines) + "\n")
    print(f"\nSummary -> {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

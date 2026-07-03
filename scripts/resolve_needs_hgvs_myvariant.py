#!/usr/bin/env python3
"""Resolve HGVSc for the 17,244 gnomad_common variants parked in needs_hgvs.tsv
using MyVariant.info batch queries (faster than VariantValidator).

For each variant, pull snpeff/VEP annotation for the project's RefSeq NM_
transcript (per the 68-gene transcript map). Output:
  - cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_resolved.tsv
    All rows with FAF >5% that we resolved -- ready to add to upload.
  - cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_unresolved.tsv
    Failures (rare -- multi-gene overlap, etc.)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

NEEDS = Path("cp_new_bundle/outputs/new_genes_68_v2/gnomad_common_needs_hgvs.tsv")
V3 = Path("cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv")
OUT_OK = Path("cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_resolved.tsv")
OUT_FAIL = Path("cp_new_bundle/outputs/new_genes_68_v3/gnomad_high_faf_unresolved.tsv")
CACHE = Path("data/cache/myvariant_resolve.json")
CACHE.parent.mkdir(parents=True, exist_ok=True)

MYVARIANT_BATCH = "https://myvariant.info/v1/variant"
BATCH_SIZE = 200
FAF_THRESHOLD = 0.05  # only resolve rows >5% per lab's request


def load_cache() -> dict:
    if CACHE.exists():
        try: return json.loads(CACHE.read_text())
        except Exception: return {}
    return {}


def save_cache(c: dict):
    CACHE.write_text(json.dumps(c))


def build_hgvs_g(chrom, pos, ref, alt) -> str:
    """Build hg38 g. notation. SNV: chr5:g.80678827T>C. Ins: chr5:g.80678827_80678828insXX. Del: chr5:g.80678827delT."""
    pos = int(pos)
    if len(ref) == 1 and len(alt) == 1:
        return f"chr{chrom}:g.{pos}{ref}>{alt}"
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        # insertion
        return f"chr{chrom}:g.{pos}_{pos+1}ins{alt[1:]}"
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        # deletion
        return f"chr{chrom}:g.{pos+1}_{pos+len(ref)-1}del"
    # complex / MNV
    return f"chr{chrom}:g.{pos}_{pos+len(ref)-1}delins{alt}"


def query_batch(ids: list[str], cache: dict) -> dict:
    """POST batch query to myvariant. Returns id -> result."""
    fresh = [i for i in ids if i not in cache]
    if fresh:
        for attempt in range(3):
            try:
                r = requests.post(MYVARIANT_BATCH, data={
                    "ids": ",".join(fresh),
                    "fields": "snpeff,_id",
                    "assembly": "hg38",
                }, timeout=60)
                r.raise_for_status()
                data = r.json()
                for item in data:
                    qid = item.get("query") or item.get("_id")
                    if qid:
                        cache[qid] = item
                break
            except Exception as e:
                print(f"  batch error attempt {attempt+1}: {e}", file=sys.stderr)
                time.sleep(2 ** attempt)
        else:
            for i in fresh: cache[i] = {"notfound": True, "error": "batch_failed"}
    return {i: cache.get(i, {"notfound": True}) for i in ids}


def extract_hgvsc(result: dict, target_nm: str) -> tuple[str, str, str] | None:
    """Walk snpeff.ann[], return (hgvs_c, effect, matched_nm) for the target NM_.
    Falls back to any NM_ for the same gene if exact-version match fails."""
    if not result or result.get("notfound"): return None
    snpeff = result.get("snpeff", {})
    if not snpeff: return None
    ann = snpeff.get("ann", [])
    if isinstance(ann, dict): ann = [ann]
    target_base = target_nm.split(".")[0]  # NM_002439 from NM_002439.5
    # Pass 1: exact NM_ match
    for a in ann:
        fid = a.get("feature_id", "")
        if fid == target_nm and a.get("hgvs_c"):
            return a["hgvs_c"], a.get("effect", ""), fid
    # Pass 2: same NM_ base (different version)
    for a in ann:
        fid = a.get("feature_id", "")
        if fid.split(".")[0] == target_base and a.get("hgvs_c"):
            return a["hgvs_c"], a.get("effect", ""), fid
    return None


def main() -> int:
    if not NEEDS.exists(): sys.exit(f"missing {NEEDS}")
    if not V3.exists(): sys.exit(f"missing {V3}")

    # Project NM_ per gene (from v3 upload's transcript column)
    v3 = pd.read_csv(V3, sep="\t", dtype=str).fillna("")
    project_nm = {g: tx for g, tx in v3[v3["transcript"]!=""][["gene","transcript"]].drop_duplicates().values}
    print(f"Project NM_ map: {len(project_nm)} genes")

    needs = pd.read_csv(NEEDS, sep="\t", dtype=str).fillna("")
    needs["_faf"] = pd.to_numeric(needs["FAF95_grpmax"], errors="coerce")
    high = needs[needs["_faf"] >= FAF_THRESHOLD].copy()
    print(f"needs_hgvs total: {len(needs):,}")
    print(f"needs_hgvs at FAF>={FAF_THRESHOLD}: {len(high):,}")

    high["_hgvs_g"] = high.apply(lambda r: build_hgvs_g(
        r["chr_hg38"], r["pos_hg38"], r["ref"], r["alt"]), axis=1)

    cache = load_cache()
    print(f"Cache: {len(cache):,} entries pre-loaded\n")

    # Batch
    all_ids = high["_hgvs_g"].tolist()
    print(f"Querying MyVariant in batches of {BATCH_SIZE}...")
    results = {}
    t0 = time.time()
    for i in range(0, len(all_ids), BATCH_SIZE):
        batch = all_ids[i:i+BATCH_SIZE]
        res = query_batch(batch, cache)
        results.update(res)
        if (i // BATCH_SIZE) % 5 == 0:
            print(f"  batch {i//BATCH_SIZE + 1}/{(len(all_ids)+BATCH_SIZE-1)//BATCH_SIZE}  "
                  f"({i+len(batch)}/{len(all_ids)}, {time.time()-t0:.0f}s)")
        if i % 1000 == 0: save_cache(cache)
    save_cache(cache)
    print(f"All batches done in {time.time()-t0:.0f}s")

    # Extract HGVSc
    rows_ok, rows_fail = [], []
    for _, r in high.iterrows():
        gene = r["gene"]
        target_nm = project_nm.get(gene)
        res = results.get(r["_hgvs_g"], {})
        if not target_nm:
            rows_fail.append({**r.to_dict(), "_fail": "no_project_NM_for_gene"})
            continue
        hit = extract_hgvsc(res, target_nm)
        if not hit:
            rows_fail.append({**r.to_dict(), "_fail": "no_hgvsc_match"})
            continue
        hgvs_c, effect, matched_nm = hit
        rows_ok.append({
            "gene": gene,
            "transcript": matched_nm,
            "hgvs_c": hgvs_c,
            "classification": "Benign FAF >5%",
            "PhastCons100way": "",
            "PhyloP100way": "",
            "REVEL": "",
            "SpliceAI_masked": "",
            "SpliceAI_unmasked": "",
            "FAF95_grpmax": r["FAF95_grpmax"],
            "CADD_phred": "",
            "AlphaMissense_pred": "",
            "ClinVar_sig": "",
            "chr_grch37": "",
            "pos_grch37": "",
            "ref": r["ref"],
            "alt": r["alt"],
            "vv_status": "myvariant_resolved",
            "source": "gnomad_common",
            "vep_effect": effect,
            "chr_hg38": r["chr_hg38"],
            "pos_hg38": r["pos_hg38"],
        })

    print(f"\n=== Resolved: {len(rows_ok):,}  | Failed: {len(rows_fail):,} ===")
    if rows_ok:
        pd.DataFrame(rows_ok).to_csv(OUT_OK, sep="\t", index=False)
        print(f"  -> {OUT_OK}")
    if rows_fail:
        pd.DataFrame(rows_fail).to_csv(OUT_FAIL, sep="\t", index=False)
        print(f"  -> {OUT_FAIL}")

    if rows_ok:
        ok_df = pd.DataFrame(rows_ok)
        print(f"\n=== Resolved by effect (top 10) ===")
        print(ok_df["vep_effect"].value_counts().head(10).to_string())
        print(f"\n=== Resolved by gene (top 10) ===")
        print(ok_df["gene"].value_counts().head(10).to_string())
        print(f"\n=== Resolved indels (len(ref) != len(alt)) ===")
        ok_df["_is_indel"] = ok_df.apply(lambda r: len(r["ref"]) != len(r["alt"]), axis=1)
        print(f"  indels: {ok_df['_is_indel'].sum():,}  / {len(ok_df):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

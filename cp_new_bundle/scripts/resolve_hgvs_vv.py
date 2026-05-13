#!/usr/bin/env python3
"""Resolve HGVS c. for cp_new SeqNext rows via VariantValidator REST.

For each row in cp_new_seqnext.tsv (and per-gene *_seqnext.tsv files),
query rest.variantvalidator.org with the GRCh37 coordinates and the BED's
RefSeq NM_ transcript. If the variant resolves on that transcript, replace
the existing hgvs_c (which came from whichever dbNSFP/snpEff transcript
happened to be in the row) with the authoritative one.

If the variant doesn't resolve on the requested transcript (e.g. non-coding
position on that isoform), mark the row with the existing c. plus a
'(VV_FLAGGED:reason)' suffix so manual review can catch it.

Results cached to disk at data/exports/cp_new/.vv_cache.json so re-runs are
incremental.

Usage:
  uv run python scripts/resolve_hgvs_vv.py [--in-place]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import requests
import pandas as pd
from tqdm import tqdm

VV_BASE = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
CACHE_PATH = Path("data/exports/cp_new/.vv_cache.json")
SEQNEXT_DIR = Path("data/exports/cp_new/seqnext")
GENOME_BUILD = "GRCh37"


def cache_load() -> dict:
    if CACHE_PATH.exists():
        with open(CACHE_PATH) as f:
            return json.load(f)
    return {}


def cache_save(cache: dict) -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_PATH, "w") as f:
        json.dump(cache, f)


def resolve_one(chrom: str, pos: int, ref: str, alt: str, tx: str,
                cache: dict, session: requests.Session) -> tuple[str, str, bool]:
    """Return (hgvs_c, status, hit_api) where status is 'ok' | 'flagged:<reason>'."""
    key = f"{chrom}-{pos}-{ref}-{alt}|{tx}"
    if key in cache:
        c, status = cache[key]
        return c, status, False

    url = f"{VV_BASE}/{GENOME_BUILD}/{chrom}-{pos}-{ref}-{alt}/{tx}"
    try:
        r = session.get(url, timeout=30)
        if r.status_code != 200:
            cache[key] = ("", f"flagged:http_{r.status_code}")
            return cache[key][0], cache[key][1], True
        data = r.json()
    except (requests.RequestException, ValueError) as e:
        return "", f"flagged:network_{type(e).__name__}", True

    # Successful response: top-level keys are like "NM_007294.4:c.5486A>T"
    # plus metadata keys ("flag", "metadata", "validation_warning_1" etc.)
    pattern = re.compile(rf"^{re.escape(tx)}:c\.")
    hits = [k for k in data.keys() if pattern.match(k)]
    if hits:
        match = hits[0]
        c_part = match.split(":", 1)[1]   # "c.5486A>T"
        cache[key] = (c_part, "ok")
        return cache[key][0], cache[key][1], True

    any_nm = [k for k in data.keys() if k.startswith("NM_") and ":c." in k]
    if any_nm:
        c_part = any_nm[0].split(":", 1)[1]
        cache[key] = (c_part, f"flagged:tx_mismatch:{any_nm[0].split(':',1)[0]}")
        return cache[key][0], cache[key][1], True

    flag = data.get("flag") or "no_match"
    cache[key] = ("", f"flagged:{flag}")
    return cache[key][0], cache[key][1], True


def annotate_file(path: Path, cache: dict, session: requests.Session, in_place: bool) -> tuple[int, int, int]:
    """Returns (ok, flagged, errors)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty or "transcript" not in df.columns:
        return 0, 0, 0

    new_c = []
    statuses = []
    for _, row in df.iterrows():
        chrom = str(row["chr_grch37"]).removeprefix("chr")
        pos = int(row["pos_grch37"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        tx = str(row["transcript"])
        c, status, hit_api = resolve_one(chrom, pos, ref, alt, tx, cache, session)
        if status == "ok":
            new_c.append(c)
        elif status.startswith("flagged:"):
            orig = str(row.get("hgvs_c") or "")
            # Strip any prior VV_FLAGGED suffix so re-runs don't accumulate them
            orig = re.sub(r"\s*\(VV_FLAGGED:[^)]*\)\s*$", "", orig)
            new_c.append(f"{orig} (VV_FLAGGED:{status[len('flagged:'):]})")
        else:
            new_c.append(row.get("hgvs_c") or "")
        statuses.append(status)
        # Rate-limit only when we actually hit the API. 1.2 s/call keeps
        # us well under VV's per-IP limit; cache hits skip the sleep.
        if hit_api:
            time.sleep(1.2)

    df["hgvs_c"] = new_c
    df["vv_status"] = statuses

    out_path = path if in_place else path.with_suffix(".vv.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    ok = sum(1 for s in statuses if s == "ok")
    flagged = sum(1 for s in statuses if s.startswith("flagged:"))
    return ok, flagged, len(statuses) - ok - flagged


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-place", action="store_true",
                   help="Overwrite the input *_seqnext.tsv files instead of writing *.vv.tsv")
    p.add_argument("--only-combined", action="store_true",
                   help="Process only cp_new_seqnext.tsv (skip per-gene files)")
    args = p.parse_args()

    if not SEQNEXT_DIR.is_dir():
        sys.exit(f"SeqNext directory not found: {SEQNEXT_DIR}")

    cache = cache_load()
    print(f"VV cache entries on disk: {len(cache):,}")

    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    # Process the combined file first; per-gene files get hits from cache after that
    combined = SEQNEXT_DIR / "cp_new_seqnext.tsv"
    files = [combined] if args.only_combined else [combined] + sorted(
        f for f in SEQNEXT_DIR.glob("*_seqnext.tsv") if f.name != combined.name
    )

    total_ok = total_flag = total_err = 0
    for f in tqdm(files, desc="files"):
        ok, flag, err = annotate_file(f, cache, session, args.in_place)
        total_ok += ok
        total_flag += flag
        total_err += err
        cache_save(cache)  # checkpoint after each file

    print(f"\nDone. ok={total_ok:,}  flagged={total_flag:,}  errors={total_err}")
    print(f"Cache saved with {len(cache):,} entries -> {CACHE_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

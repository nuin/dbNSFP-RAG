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

Backends (--backend):
  rest   (default) rest.variantvalidator.org, throttled to ~1 call / 1.2 s.
  local  biocommons `hgvs` (the engine VariantValidator is built on) against
         a UTA database + SeqRepo; no REST calls, no rate limit.
         Install with `uv pip install -e ".[hgvs]"`. `hgvs` reads UTA_DB_URL
         (default: public uta.biocommons.org; point it at a local UTA
         Postgres for bulk runs) and HGVS_SEQREPO_DIR (local SeqRepo
         snapshot; otherwise sequences are fetched remotely).
         Cached separately in data/exports/cp_new/.hgvs_local_cache.json.

Usage:
  uv run python scripts/resolve_hgvs_vv.py [--in-place] [--backend rest|local]
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
LOCAL_CACHE_PATH = Path("data/exports/cp_new/.hgvs_local_cache.json")
SEQNEXT_DIR = Path("data/exports/cp_new/seqnext")
GENOME_BUILD = "GRCh37"


def cache_load(path: Path = CACHE_PATH) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def cache_save(cache: dict, path: Path = CACHE_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
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


def vcf_to_hgvs_g(ac: str, pos: int, ref: str, alt: str) -> str:
    """Convert a VCF-style allele to an HGVS g. string on RefSeq accession `ac`.

    Trims shared trailing/leading bases so padded indels become del/ins/delins;
    the hgvs normalizer 3'-shifts the projected c. afterwards.
    """
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    if len(ref) == 1 and len(alt) == 1:
        return f"{ac}:g.{pos}{ref}>{alt}"
    if not ref:  # pure insertion between pos-1 and pos
        return f"{ac}:g.{pos - 1}_{pos}ins{alt}"
    end = pos + len(ref) - 1
    span = f"{pos}" if end == pos else f"{pos}_{end}"
    if not alt:
        return f"{ac}:g.{span}del"
    return f"{ac}:g.{span}delins{alt}"


class LocalResolver:
    """Project GRCh37 g. onto an NM_ with biocommons hgvs (UTA + SeqRepo)."""

    def __init__(self, assembly: str = GENOME_BUILD):
        import hgvs.assemblymapper
        import hgvs.dataproviders.uta
        import hgvs.exceptions
        import hgvs.parser
        from bioutils.assemblies import make_name_ac_map

        self.exc = hgvs.exceptions
        self.hdp = hgvs.dataproviders.uta.connect()
        self.parser = hgvs.parser.Parser()
        self.am = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, assembly_name=assembly, alt_aln_method="splign",
            replace_reference=True,
        )
        self.name_to_ac = make_name_ac_map("GRCh37.p13" if assembly == "GRCh37" else assembly)

    def g_to_c(self, chrom: str, pos: int, ref: str, alt: str, tx: str) -> str:
        """Return 'NM_...:c.XXX' or raise ValueError(<flag reason>)."""
        ac = self.name_to_ac.get(chrom)
        if ac is None:
            raise ValueError(f"unknown_chrom_{chrom}")
        try:
            var_g = self.parser.parse_hgvs_variant(vcf_to_hgvs_g(ac, pos, ref, alt))
            var_c = self.am.g_to_c(var_g, tx)
        except self.exc.HGVSDataNotAvailableError as e:
            raise ValueError("tx_not_available") from e
        except self.exc.HGVSInvalidIntervalError as e:
            raise ValueError("intergenic") from e
        except self.exc.HGVSError as e:
            raise ValueError(type(e).__name__) from e
        return str(var_c)


def resolve_one_local(chrom: str, pos: int, ref: str, alt: str, tx: str,
                      cache: dict, resolver) -> tuple[str, str, bool]:
    """Local-backend twin of resolve_one(): same (hgvs_c, status, hit_api) contract."""
    key = f"{chrom}-{pos}-{ref}-{alt}|{tx}"
    if key in cache:
        c, status = cache[key]
        return c, status, False
    try:
        full = resolver.g_to_c(chrom, pos, ref, alt, tx)
    except ValueError as e:
        cache[key] = ("", f"flagged:{e}")
        return cache[key][0], cache[key][1], True
    cache[key] = (full.split(":", 1)[1], "ok")
    return cache[key][0], cache[key][1], True


def annotate_file(path: Path, cache: dict, resolve, in_place: bool,
                  rate_limit: float = 1.2) -> tuple[int, int, int]:
    """Returns (ok, flagged, errors).

    `resolve(chrom, pos, ref, alt, tx, cache)` -> (hgvs_c, status, hit_api).
    """
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
        c, status, hit_api = resolve(chrom, pos, ref, alt, tx, cache)
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
        if hit_api and rate_limit:
            time.sleep(rate_limit)

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
    p.add_argument("--backend", choices=("rest", "local"), default="rest",
                   help="rest = rest.variantvalidator.org (default); "
                        "local = biocommons hgvs + UTA/SeqRepo (needs the [hgvs] extra)")
    args = p.parse_args()

    if not SEQNEXT_DIR.is_dir():
        sys.exit(f"SeqNext directory not found: {SEQNEXT_DIR}")

    cache_path = CACHE_PATH if args.backend == "rest" else LOCAL_CACHE_PATH
    cache = cache_load(cache_path)
    print(f"{args.backend} cache entries on disk: {len(cache):,}")

    if args.backend == "rest":
        session = requests.Session()
        session.headers.update({"Accept": "application/json"})

        def resolve(chrom, pos, ref, alt, tx, cache):
            return resolve_one(chrom, pos, ref, alt, tx, cache, session)
        rate_limit = 1.2
    else:
        try:
            resolver = LocalResolver()
        except ImportError as e:
            sys.exit(f'--backend local needs the hgvs extra (uv pip install -e ".[hgvs]"): {e}')

        def resolve(chrom, pos, ref, alt, tx, cache):
            return resolve_one_local(chrom, pos, ref, alt, tx, cache, resolver)
        rate_limit = 0.0

    # Process the combined file first; per-gene files get hits from cache after that
    combined = SEQNEXT_DIR / "cp_new_seqnext.tsv"
    files = [combined] if args.only_combined else [combined] + sorted(
        f for f in SEQNEXT_DIR.glob("*_seqnext.tsv") if f.name != combined.name
    )

    total_ok = total_flag = total_err = 0
    for f in tqdm(files, desc="files"):
        ok, flag, err = annotate_file(f, cache, resolve, args.in_place, rate_limit)
        total_ok += ok
        total_flag += flag
        total_err += err
        cache_save(cache, cache_path)  # checkpoint after each file

    print(f"\nDone. ok={total_ok:,}  flagged={total_flag:,}  errors={total_err}")
    print(f"Cache saved with {len(cache):,} entries -> {cache_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

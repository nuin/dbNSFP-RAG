#!/usr/bin/env python3
"""Validate our SpliceAI scores against the Broad SpliceAI Lookup API.

RUN THIS ON A MACHINE WITH INTERNET ACCESS to spliceailookup-api.broadinstitute.org
(it is egress-blocked from the build sandbox, so it must be run locally).

For a sample of exported variants it queries the API in BOTH masked (mask=1) and
unmasked (mask=0) modes at distance=50 (the same settings our local run used) and
compares to our SpliceAI_masked / SpliceAI_unmasked columns.

Usage:
  python scripts/validate_spliceai_api.py            # 50-variant sample
  python scripts/validate_spliceai_api.py 200        # 200-variant sample
  python scripts/validate_spliceai_api.py 200 KMT2D  # restrict to a gene
"""

import json
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

UPLOAD = "cp_new_bundle/outputs/new_genes_68_v3_7/UPLOAD_v3_7.tsv"
API = (
    "https://spliceailookup-api.broadinstitute.org/spliceai/"
    "?hg=38&distance=50&mask={mask}&variant={variant}"
)
TOL = 0.02  # allow 0.02 rounding tolerance


class InvalidSpliceAIResponse(TypeError):
    """SpliceAI returned JSON that does not match its score response schema."""


def response_ds_max(data):
    if not isinstance(data, dict):
        raise InvalidSpliceAIResponse("SpliceAI response must be an object")
    scores = data.get("scores")
    if scores is None:
        return 0.0
    if not isinstance(scores, list):
        raise InvalidSpliceAIResponse("SpliceAI scores must be a list")
    best = 0.0
    for score in scores:
        if not isinstance(score, dict):
            raise InvalidSpliceAIResponse("SpliceAI score entries must be objects")
        for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
            value = score.get(key)
            if value is None:
                continue
            try:
                parsed = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise InvalidSpliceAIResponse(f"Invalid SpliceAI {key} score: {value!r}") from exc
            best = max(best, parsed)
    return best


def api_ds_max(gnomad_id, mask):
    url = API.format(mask=mask, variant=gnomad_id)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                d = json.load(r)
            return response_ds_max(d)
        except (
            urllib.error.URLError,
            TimeoutError,
            OSError,
            ValueError,
            InvalidSpliceAIResponse,
        ) as exc:
            print(
                f"SpliceAI API attempt {attempt + 1} failed for {gnomad_id}: {exc}", file=sys.stderr
            )
            if attempt == 2:
                return None
            time.sleep(3)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    gene = sys.argv[2] if len(sys.argv) > 2 else None
    d = pd.read_csv(UPLOAD, sep="\t", dtype=str).fillna("")
    if gene:
        d = d[d["gene"] == gene]
    # only rows with a measured SpliceAI value and a gnomad_id
    d = d[(d["SpliceAI_masked"] != "") & (d["gnomad_id"].str.count("-") == 3)]
    step = max(1, len(d) // n)
    sample = d.iloc[::step].head(n).reset_index(drop=True)
    print(f"Validating {len(sample)} variants against Broad SpliceAI Lookup API\n")

    ok_m = ok_u = bad = err = 0
    print(f"{'variant':<26}{'ours_m':>8}{'api_m':>8}{'ours_u':>8}{'api_u':>8}  status")
    print("-" * 70)
    for _, r in sample.iterrows():
        am = api_ds_max(r["gnomad_id"], 1)
        au = api_ds_max(r["gnomad_id"], 0)
        om, ou = num(r["SpliceAI_masked"]), num(r["SpliceAI_unmasked"])
        if am is None or au is None:
            err += 1
            status = "API error"
        else:
            mmatch = om is not None and abs(om - am) <= TOL
            umatch = ou is not None and abs(ou - au) <= TOL
            ok_m += mmatch
            ok_u += umatch
            status = "OK" if (mmatch and umatch) else "DIFF"
            if status == "DIFF":
                bad += 1
        print(
            f"{r['gene'] + ' ' + r['hgvs_c']:<26}"
            f"{('' if om is None else f'{om:.2f}'):>8}{('' if am is None else f'{am:.2f}'):>8}"
            f"{('' if ou is None else f'{ou:.2f}'):>8}{('' if au is None else f'{au:.2f}'):>8}  {status}"
        )
        time.sleep(1)  # be polite to the API

    print(f"\n=== {len(sample)} variants ===")
    print(f"  masked   within {TOL}: {ok_m}")
    print(f"  unmasked within {TOL}: {ok_u}")
    print(f"  differing: {bad} | API errors: {err}")


if __name__ == "__main__":
    main()

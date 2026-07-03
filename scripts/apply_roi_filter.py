#!/usr/bin/env python3
"""Apply the intronic ROI [-20,+10] filter to a full upload table -> _roi.tsv.
Keeps coding/UTR-exonic (no offset) + intronic within [-20,+10]; drops deep
intronic. Usage: apply_roi_filter.py <upload.tsv>
"""
import re, sys
from pathlib import Path
import pandas as pd

OFF = re.compile(r"c\.[\*\-]?\d+([+\-])(\d+)")
LO, HI = -20, 10

def offset(h):
    m = OFF.search(h or "")
    if not m: return None
    return -int(m.group(2)) if m.group(1)=="-" else int(m.group(2))

def main():
    src = Path(sys.argv[1])
    df = pd.read_csv(src, sep="\t", dtype=str).fillna("")
    off = df["hgvs_c"].apply(offset)
    keep = off.isna() | off.between(LO, HI)
    kept, dropped = df[keep], df[~keep]
    out = src.with_name(src.stem + "_roi.tsv")
    kept.to_csv(out, sep="\t", index=False)
    dropped.to_csv(src.with_name(src.stem + "_roi_dropped.tsv"), sep="\t", index=False)
    print(f"{src.name}: {len(df):,} -> {len(kept):,} ROI [-20,+10] ({len(dropped):,} deep-intronic dropped)")

if __name__ == "__main__":
    main()

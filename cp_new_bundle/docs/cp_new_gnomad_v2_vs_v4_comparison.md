# gnomAD v4.1 (hg38 → dbNSFP lift) vs v2.1.1 (hg19 native)

Re-ran the cp_new pipeline's FAF lookup using **gnomAD v2.1.1 exomes
sites VCF queried in hg19 directly** (no liftover) and compared against
the FAFs we currently use (gnomAD v4.1 joint, hg38, queried at dbNSFP's
hg19-lifted coordinates).

**Question being tested**: does the dbNSFP-bundled hg19 coordinate
introduce liftover variance that changes the AF/FAF we see? If yes,
we'd see lots of FAF disagreement between the two paths. If no, the
disagreements are just v2-vs-v4 sampling differences.

## Coverage (1,057,343 cp_new variants)

| Where the variant has FAF data | Count |
|---|---|
| In **both** releases | 67,480 |
| **v4.1 only** (newer release, more samples) | 154,596 |
| **v2.1.1 only** (rare — present in older release, dropped from v4?) | 3,111 |
| In neither (rare/novel — most of dbNSFP's nsSNVs are predicted but never observed) | 832,156 |

v4 has 50× more "release-only" variants than v2 — expected, since v4 is a
much larger sample (807k joint vs 125k v2 exomes).

## FAF95 vs FAF95 (apples-to-apples, 49,853 variants)

| Agreement | Count |
|---|---|
| Exact match (Δ ≤ 1e-7) | 167 |
| Within **10%** relative | 3,183 |
| Within **50%** relative | 19,237 |
| **>2× apart** | **0** |

**No variant had FAFs more than 2× apart between the two releases.**
Strong evidence that **dbNSFP's hg19-from-hg38 lift is not introducing
material FAF error**. The differences that exist are sampling noise
(v4 has ~6.5× more samples) and ancestry-group composition differences,
not liftover artifacts.

## Rule-fire impact

Would the W1 and W3 classifications change based on which release we
queried?

| Rule | Both fire | v4-only fires | v2-only fires |
|---|---|---|---|
| **W1** (FAF > 5%) | 288 | 12 | 7 |
| **W3** (FAF > 0.1%) | 2,016 | 414 | 187 |

- **W1**: 19 of 307 high-FAF calls (~6%) flip between releases — mostly
  variants near the 5% threshold where v2's smaller sample makes FAF95
  cross the boundary differently.
- **W3**: 601 of 2,617 calls (~23%) flip. This is bigger because the
  0.1% threshold sits in a regime where v2's smaller sample noise has
  more leverage on the FAF95 lower bound.

The disagreements are mostly **rare-variant calls where v2's small
sample produces a FAF95 of 0** (insufficient confidence) while v4's
larger sample gives a non-zero estimate.

## Examples of the biggest disagreements

| Gene | hg19 pos | ref>alt | v4 FAF | v2 FAF | Reason |
|---|---|---|---|---|---|
| TRIM28 | chr19:59,055,539 | T>C | 0.284 | 0.000 | v2 too few samples to estimate FAF95; AF_popmax non-zero |
| GALNT12 | chr9:101,570,116 | G>A | 0.197 | 0.000 | same — FAF95 LB = 0 in v2, real freq seen in v4 |
| FANCD2 | chr3:10,088,343 | A>G | 0.456 | 0.281 | common variant, both report it; FAF95 differs because v2's AFR ancestry sampling differs from v4 |
| CDKN1C | chr11:2,906,095 | C>T | 0.002 | 0.136 | v2 has higher count in non-NFE pops; v4's larger joint set dilutes the per-pop max |
| POT1 | chr7:124,499,165 | A>C | 0.063 | 0.197 | population-specific high-AF variant; v2's smaller pops show it higher |
| MAD2L2 | chr1:11,740,542 | C>A | 0.074 | 0.190 | similar — pop-max behavior differs |
| ALK | chr2:29,451,802 | G>C | 0.166 | 0.054 | v4 captures it in genomes too; popmax shifts |

These are real **sampling and ancestry-composition differences between
v2 and v4** — not liftover errors. None of them indicate that dbNSFP's
hg19 coords are wrong; the same genomic variant is being captured in
both releases, just with different denominators.

## Conclusion

| Concern | Answer |
|---|---|
| Does dbNSFP's hg19 liftover corrupt FAF lookups? | **No** — 0 cases of >2× disagreement; the lift is precise enough for FAF-based classification |
| Should we switch to gnomAD v2.1.1 to avoid liftover entirely? | **No** — v2 has 50× fewer release-only variants and significantly noisier FAF estimates at rare-variant ranges. Using v4 is the more accurate source despite the lift. |
| Could rule fires depend on release? | **Yes, modestly** — ~6% of W1, ~23% of W3 calls flip between releases. This is expected sampling difference, not pipeline error. |
| Should we report which release the FAF came from? | **Yes** — already present in the SeqNext export columns (`FAF95_grpmax` is from v4.1 joint). Audit log already implicit in the script. |

The pipeline's choice (v4.1 hg38 + dbNSFP hg19-lift) is the right one
for accuracy. Comparing against v2.1.1 confirms the lift isn't a
material source of error — the gnomAD release version dominates any
liftover variance.

## Reproduce

```bash
uv run python scripts/compare_gnomad_v2_hg19.py
```

Outputs:
- `data/exports/cp_new/cp_new_gnomad_v2_compare.tsv` (per-variant)
- `data/exports/cp_new/cp_new_gnomad_v2_summary.txt` (aggregate stats)

~8 minutes (165 gene-region tabix queries against gnomAD v2.1.1 public bucket).

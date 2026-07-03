# v3.3 — gnomAD-anchored filtering (lab feedback 2026-06-19)

Lab: "Still need more filtering. Nonsense variants in the list, pipeline fails
need removing, the W2 increase is the issue, anchor on gnomad."

## Filters applied to v3.2

| Filter | Rows dropped | Rationale |
|---|---:|---|
| **gnomAD anchor** | 67,992 | only keep variants observed in gnomAD v4.1 |
| **rule_fails removed** | 3,539 | pipeline fails not imported to SeqNext |
| **nonsense removed** | 10 | pipeline codon-degeneracy mis-tags |

## v3.2 → v3.3

| | v3.2 | v3.3 |
|---|---:|---:|
| Total | 113,055 | **41,514** |
| Benign W1 (FAF>5%) | 11,714 | 11,711 |
| Benign W2 (syn) | 96,109 | 29,316 |
| Likely_benign W3 | 495 | 487 |
| rule_fails | 3,539 | 0 |
| gene coverage | 67 | 67 |

The W2 collapse (96k → 29k) is the gnomAD anchor removing ~67k synonymous
variants from the local enumeration catalog that were never observed in any
real individual. Every remaining row is a variant that actually exists in
gnomAD v4.1.

(TERC remains the only zero-coverage gene — a non-coding RNA where no W1/W2/W3
rule applies.)

## The gnomAD anchor

Pulled ALL gnomAD v4.1 observed variants (any frequency) for the 68 gene
regions: **1,401,013 observed variants** across all 67 coding genes
(`data/exports/cp_new/gnomad_observed/`). A row stays only if
`(chr_hg38, pos_hg38, ref, alt)` is in that set. gnomad_common-source rows are
observed by definition (they came from gnomAD).

NOTE: 4 large gene regions (KIF1B, KMT2D, LZTR1, MAD2L2) silently failed their
first tabix pull (empty result, no error) and had to be re-pulled — otherwise
all their variants would have been falsely dropped as "not in gnomAD". Worth a
non-empty-result assertion in the pull script for next time.

## Two bugs found & fixed during this build

### Bug A — anchor falsely dropped real synonymous variants
First attempt translated the upload's hg19 coords to hg38 using dbNSFP's
position map. But dbNSFP excludes synonymous variants by design, so synonymous
positions weren't in the map and couldn't be translated — they were all
falsely dropped (only 16/24 MSH3 lab variants survived).
**Fix**: use pyliftover (UCSC hg19→hg38 chain) which handles all positions.
Recovered 7 of the 8 lost MSH3 syn.

### Bug B — nonsense filter conflated transcripts
The nonsense filter keyed on `(gene, hgvs_c)`. The same c. string can name two
different genomic variants on different transcripts — e.g. MSH3 c.1992G>A is
synonymous on the canonical NM_002439.5 (hg38 80768028) but ALSO appears as
p.Trp664* on an alternate transcript for a DIFFERENT genomic variant (hg38
80775744, c.2304G>A canonical). The loose key dropped the synonymous one.
**Fix**: exempt synonymous_catalog rows from the nonsense filter — they're
BioPython-verified synonymous on the project NM_ and cannot be nonsense.
The nonsense filter now drops only 10 true pipeline mis-tags.

## MSH3 verification

All 24 of the lab's expected MSH3 variants are present in v3.3:
- 5 high-FAF dels/dups + 1 intronic (gnomad_common, Benign FAF >5%)
- 11 missense (pipeline, Likely_benign W3)
- 8 synonymous (synonymous_catalog, Benign W2)

MSH3 total: 1,547 rows (all gnomAD-observed).

## Files
```
UPLOAD_v3.tsv                  35,962 rows -- the SeqNext upload (v3.3)
UPLOAD_v3_3_dropped.tsv        77,093 rows -- audit with drop reason per row
per_gene/{GENE}__v3.tsv        66 files
```

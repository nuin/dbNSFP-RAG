**Subject:** cp_new v3.3 — gnomAD-anchored, all filters applied

All four of your filtering requests are in. v3.3 at
`cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv` — **41,514 rows**.

**1. Anchored on gnomAD.** Pulled every gnomAD v4.1 variant (any frequency) in
the 68 gene regions — 1.4M observed variants — and kept only upload rows that
actually appear in gnomAD. This removed ~68,000 synonymous variants from my
local enumeration that were never seen in any real individual. That was the
"W2 increase" you flagged: most of it was theoretical.

**2. Pipeline fails (rule_fails) removed** — all 3,539 dropped.

**3. Nonsense removed** — but only 10, all true pipeline mis-tags. (My first
pass dropped a real synonymous variant, MSH3 c.1992G>A, because it shares the
"c.1992G>A" string with a nonsense variant on a different transcript at a
different genomic position. Fixed.)

**v3.2 → v3.3**

| | v3.2 | v3.3 |
|---|---:|---:|
| Total | 113,055 | 41,514 |
| Benign W1 (FAF>5%) | 11,714 | 11,711 |
| Benign W2 (synonymous) | 96,109 | 29,316 |
| Likely_benign W3 | 495 | 487 |
| rule_fails | 3,539 | 0 |

**MSH3 — all 24 of your listed variants are present** (5 dels/dups + 1 intronic
as Benign FAF>5%, 11 missense as W3 LB, 8 synonymous as W2 Benign). MSH3 total
1,547 rows, all gnomAD-observed.

**On synonymous_catalog** — to answer your earlier question, it's a local
enumeration of every possible synonymous SNV in each gene's CDS (built because
dbNSFP excludes synonymous and ClinVar only has submitted ones). The gnomAD
anchor now restricts it to synonymous variants that have actually been observed,
which is what you wanted — no more theoretical calls that aren't in gnomAD.

Two coordinate bugs found and fixed during this build (both documented in
`docs/V3_3_FILTERS.md`):
- the gnomAD anchor was using dbNSFP's coordinate map, which excludes
  synonymous positions, so it was falsely dropping real syn variants — switched
  to proper UCSC liftover
- 4 large gene regions silently failed their gnomAD pull and had to be
  re-pulled, else their variants would have been falsely dropped

Every row in v3.3 is a variant that exists in gnomAD, classified purely by
W1/W2/W3, no ClinVar evidence used. Please spot-check MSH3 and any other gene.

Paulo

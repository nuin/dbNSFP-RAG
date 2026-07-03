# v3.4 — liftover-free gnomAD anchor (c.-based)

Built in response to: "the c. on gnomAD 4 is enough" — the lab's point that
genome-build liftover isn't needed because HGVS c. is transcript-relative.

They're right, with one caveat: gnomAD's *displayed* c. defaults to MANE Select,
which isn't always the project NM_. So v3.4 extracts gnomAD's HGVSc for the
specific MANE_SELECT NM_ from gnomAD v4.1's own VEP annotation and matches on
`(gene, NM_, c.)`. This is:
- **liftover-free** (c. is build-independent)
- **transcript-correct** (matched to the project NM_, not gnomAD's default)
- **collision-proof** (each canonical c. is unique to its genomic variant —
  fixes the c.1992G>A cross-transcript bug that broke position-matching)

## How the anchor works

1. `pull_gnomad_hgvsc.py` pulls gnomAD v4.1 exomes+genomes VEP for each gene
   region and extracts `(SYMBOL, MANE_SELECT_NM, c.)` for every observed variant
   → 2,747,960 keys (`data/exports/cp_new/gnomad_hgvsc/_observed_cdot.tsv`).
2. `build_upload_v3_4.py` keeps an upload row only if its
   `(gene, NM_, c.)` is in that set.

### Transcript-mismatch fallback (4 genes)

For 4 genes the project NM_ differs from gnomAD's MANE_SELECT, so c. numbering
can't be reconciled:

| gene | project NM_ | gnomAD MANE |
|---|---|---|
| DPYD | ENST00000876341 | NM_000110 |
| GPR161 | NM_001349632.1 | NM_001375883 |
| MITF | NM_000248.3 | NM_001354604 |
| SMARCA4 | NM_001128849.2 | NM_003072 |

For just these 4, v3.4 falls back to position-matching (pyliftover) against the
gnomAD observed-coordinate set. The other 63 genes are fully liftover-free.

## v3.4 numbers

| | rows |
|---|---:|
| **Total** | **41,450** |
| Benign W1 (FAF>5%) | 11,709 |
| Benign W2 (synonymous) | 29,316 |
| Likely_benign W3 | 425 |

Source: synonymous_catalog 29,303 / gnomad_common 11,481 / pipeline 345 /
clinvar 321. Gene coverage 67/68 (TERC is non-coding RNA, no rule applies).

## vs v3.3 (position-anchored)

| | v3.3 (position+liftover) | v3.4 (c.-based) |
|---|---:|---:|
| Total | 41,514 | 41,450 |
| Liftover used | yes (all genes) | only 4 mismatch genes |
| Anchor key | (chr,pos,ref,alt) hg38 | (gene, NM_, c.) |

Nearly identical totals — confirms both anchors agree. v3.4 is the cleaner
design: liftover-free for 63/67 coding genes and immune to cross-transcript
c.-collisions.

## MSH3 verification

24/25 lab variants present by exact c.; the 25th (c.181_189del) is present as
its sibling c.181_189dup — gnomAD carries both del and dup at that 9-bp repeat,
and the upstream gnomAD pull captured the dup allele. Same locus, FAF ~0.19.

## Files
```
UPLOAD_v3_4.tsv              41,450 rows -- the SeqNext upload
UPLOAD_v3_4_dropped.tsv      71,605 rows -- audit with drop reason
per_gene/{GENE}__v3_4.tsv    67 files
docs/README.md               this file
```

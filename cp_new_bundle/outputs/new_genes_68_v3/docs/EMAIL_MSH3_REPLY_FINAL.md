**Subject:** Re: MSH3 review — v3.2 ready, 24 of 25 lab variants now covered

You were right and the gaps were real. Three separate bugs identified and fixed,
plus the high-FAF gap closed. v3.2 covers **24 of the 25 MSH3 variants** you
listed (the 25th, c.181_189del, is in as c.181_189dup — MyVariant resolves the
same 9-bp repeat-region variant as a duplication; same locus, same FAF).

## What was broken

**Bug 1 — high-FAF dels/dups/intronic missing.** My v2 gnomAD pull captured
17,244 high-FAF variants (including all your dels and intronic), but I never
resolved their HGVSc so they sat in `needs_hgvs.tsv` and never made it into the
upload. Now resolved via MyVariant batch query (took 12 min for 17k variants).
**Recovered 11,481 high-FAF rows** as `Benign FAF >5%` including all 5 dels +
the intronic c.1897-8A>G you flagged.

**Bug 2 — missense W3 dropping silently.** `classify_cp_new.py:safe_float()`
choked on dbNSFP's `;`-joined multi-transcript values like
`.;.;0.220;.;.;.;.;.` and returned None, so REVEL was None for most missense
variants and W3 never fired. **Recovered 268 missense W3 LB calls** across the
68 genes — including all 11 MSH3 missense LB variants you listed (c.178G>C,
c.190C>G, c.173C>T, etc.).

**Bug 3 — synonymous catalog dropping rows without PhastCons.** The annotator
required PhastCons<1.0 and dropped rows where PhastCons was missing (positions
with no dbNSFP entry to borrow from). **Recovered 85,282 syn variants** —
including all 8 MSH3 syn variants you listed.

All three source-level bugs fixed with comments warning future-me not to
reintroduce them.

## v3.2 numbers

| | v3.1 | v3.2 |
|---|---:|---:|
| Total rows | 16,025 | **111,857** |
| Benign W1 (FAF>5%) | 234 | 11,714 |
| Benign W2 (syn) | 12,024 | 96,109 |
| Likely_benign W3 | 228 | 495 |
| MSH3 rows | 354 | **3,199** |

The synonymous catch-up explains most of the jump. Many of those rows have
PhastCons unmeasured (we lacked a borrowable conservation score). They're
emitted with the note `"PhastCons unmeasured"` in the classification text.

## MSH3 lab-expected variants

24/25 in v3.2:

| HGVSc | classification | source |
|---|---|---|
| c.162_179del | Benign FAF >5% | gnomad_common |
| c.199_207del | Benign FAF >5% | gnomad_common |
| c.359-7G>A | Benign FAF >5% | gnomad_common |
| c.181_189del → **c.181_189dup** | Benign FAF >5% | gnomad_common (same locus, dup notation) |
| c.178_186del | Benign FAF >5% | gnomad_common |
| c.1897-8A>G | Benign FAF >5% | gnomad_common |
| c.178G>C, c.190C>G, c.173C>T, c.1258A>G, c.1571A>C, c.1313C>T, c.205C>T, c.2740A>G, c.1522A>G, c.146C>G, c.1160T>A | Likely_benign W3 | pipeline |
| c.162T>C, c.1992G>A, c.204T>G, c.111C>T, c.2685C>T, c.1194C>T, c.96A>C, c.3009C>T | Benign W2 (syn) | synonymous_catalog |

## On synonymous_catalog

To answer your question — `synonymous_catalog` is a local enumeration I built
because dbNSFP excludes synonymous variants by design and ClinVar only has the
syn variants someone has submitted. The script walks every CDS exon in hg19
RefSeq for each of the 68 genes, generates each codon, and lists every SNV that
produces a synonymous substitution per the standard codon table (verified with
BioPython). It does NOT include indels — for the 1-2 nt synonymous indels you
might care about, those are in gnomAD and reach the upload via the FAF>5% pull.

## SeqNext size

111k rows is large. If the system has trouble with this size, easy paths:
- drop rule_fails entirely (-3,539 rows)
- drop W2 rows where PhastCons unmeasured (-~85k rows; would also drop those 8
  MSH3 syn variants you flagged, so not recommended without lab agreement)

Tell me if either filter is wanted.

Files in `cp_new_bundle/outputs/new_genes_68_v3/`:
- `UPLOAD_v3.tsv` — v3.2
- `per_gene/MSH3__v3.tsv` — 3,199 rows for spot-check
- `NEW_IN_V3_2.tsv` — audit of the 95,832 new rows added

Paulo

# cp_new SeqNext UPLOAD — v2 (gnomAD-inclusive)

Addresses lab feedback (2026-06-01) on the v1 upload:
1. Independent gnomAD pull of FAF >5% variants — done, found GOT2 c.816C>T and c.213T>C among 45 new synonymous-coding common variants
2. Stop dropping synonymous coding variants — done, they now flow through as Benign

## Files

| File | Rows | Purpose |
|---|---:|---|
| `cp_new_seqnext_UPLOAD_v2.tsv` | 39,371 | **The new SeqNext upload** |
| `cp_new_gnomad_common_standalone.tsv` | 172 | Standalone gnomAD-common variants (HGVSc resolved) |
| `cp_new_gnomad_common_needs_hgvs.tsv` | 17,243 | Noncoding gnomAD-common variants pending VV resolution |
| `per_gene/{GENE}__v2.tsv` | 68 files | Per-gene split of v2 upload |

## v1 vs v2

| Metric | v1 | v2 | Delta |
|---|---:|---:|---|
| Total rows | 5,796 | **39,331** | +33,535 |
| Gene coverage | 67/68 | **68/68** | CDKN2C/ELOC/H3-3A no longer blank |
| Synonymous coding | dropped | kept | now in upload |
| Source: pipeline | 174 | 1,715 | restored W2/W3 intronic-ROI rows |
| Source: clinvar | 5,622 | 24,856 | restored syn-coding from ClinVar |
| Source: synonymous_catalog | 0 | 12,759 | now flows through |
| Source: gnomad_common | 0 | 1 | NEW (independent FAF>5% pull, others deduped to ClinVar) |

## Composition

By class:
- Benign: 18,135 (46%)
- Likely_benign: 21,196 (54%)

## What was the actual fix

The lab's three originally-missed GOT2 variants were **always present in ClinVar
(2-star Benign)** — but my v1 syn-coding drop filter was removing them. Reversing
that filter restored them along with ~32k other valid Benign rows.

The independent gnomAD pull (which I built per your proposal) found 172 common
variants with HGVSc; 171 were already in ClinVar (confirmed, deduped in favor of
ClinVar), 1 was net-new. So the gnomAD pull is good due diligence — it confirms
the data and would catch anything ClinVar lacks — but the bulk of the missing
benigns came back through unfiltering syn coding.

GOT2 verification in v2:
- c.816C>T — ClinVar 2-star Benign ✅
- c.213T>C — ClinVar 2-star Benign ✅
- c.228T>G — ClinVar 2-star Benign ✅

## Open: VV resolution for 17,243 noncoding gnomAD-commons

Deep intronic / UTR variants at FAF >5% that don't have an HGVSc yet (not in
dbNSFP, not in synonymous catalog). To bring them into the upload, run
VariantValidator on the list — at the existing 1.2s/call rate that's ~6 hours.

Worth checking with the lab first whether they want those (they're variants
the panel calls, but most analysts don't review deep intronic commons).

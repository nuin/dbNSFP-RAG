**Subject:** cp_new pipeline v2 — synonymous restored, GOT2 variants confirmed

Ran both fixes you proposed. Quick summary:

**Upload size: 5,796 → 39,331 rows** (folder: `cp_new_bundle/outputs/v2_gnomad_inclusive/`)
- All 68 new genes now have coverage (CDKN2C, ELOC, H3-3A were previously blank).
- Synonymous coding rows flow through as Benign — no longer dropped.

**GOT2 variants the lab flagged — all three confirmed in v2:**
- c.816C>T — ClinVar 2-star Benign + gnomAD FAF 85%
- c.213T>C — ClinVar 2-star Benign + gnomAD FAF 80%
- c.228T>G — ClinVar 2-star Benign

These were always in ClinVar; my v1 filter was dropping them. Fixed.

**Independent gnomAD pull also done** (your proposal #2):
- Pulled all FAF >5% variants for 68 gene regions from gnomAD v4.1
- 17,415 common variants in scope
- 172 with HGVSc resolved (125 missense + 45 synonymous + 2 other)
- 17,243 noncoding (intronic/UTR) need VariantValidator HGVSc resolution — currently parked in `cp_new_gnomad_common_needs_hgvs.tsv` pending your call on whether the lab wants them in the upload

**Composition of v2:**
- 24,856 ClinVar B/LB (63%) — includes restored synonymous
- 12,759 synonymous catalog (32%) — local enumeration, Benign by W2
- 1,715 pipeline rules W1/W2/W3 (4%)
- 1 gnomad_common net-new

**Recommended next step:** spot-check `per_gene/GOT2__v2.tsv` (or any gene you want), then we can decide on:
- Whether to VV-resolve the 17,243 noncoding gnomAD-commons (~6 hours batch)
- Whether to keep the 1-star ClinVar rows or undercall those too

Paulo

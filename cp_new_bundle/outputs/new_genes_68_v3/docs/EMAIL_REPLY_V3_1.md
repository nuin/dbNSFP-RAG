**Subject:** Re: HOXB13 review — confirmed bug, v3.1 fixed

You're right and the catch is important. Thank you both for the spot-check.

**HOXB13 c.108C>A, c.144T>A etc. were definitely mis-tagged.**

Root cause: my `classify_cp_new.py` had a shortcut treating any variant at a
dbNSFP "2-fold degenerate" codon position as synonymous. That's wrong — at a
2-fold site only 1 of 3 alts is silent, the other 2 are missense. The bug
mis-tagged 1,083 missense + 134 nonsense variants across all 68 genes.

**Fixes in v3.1** (now at `UPLOAD_v3.tsv`):

1. Re-verified every "Benign synonymous" call against dbNSFP's per-row
   aaref/aaalt comparison. Dropped 1,217 mis-tagged missense/nonsense rows.
   Three of those were borderline missense and got reclassified to
   `Likely_benign FAF >0.1%, REVEL<0.29, SpliceAI<=0.1` (W3 passes).
2. Removed all 20,991 `rule_fails:no_rule_applies` rows from the upload as
   requested. Other `rule_fails:` rows (W3_faf_too_low, W2_spliceai etc.)
   are still in the file as a review queue, but if you want all rule_fails
   dropped just say the word.
3. Fixed the bug at source in `scripts/classify_cp_new.py` with a comment
   warning so it doesn't recur.

**v3 → v3.1**: 38,314 → 16,025 rows.

HOXB13 spot-check confirmation:
- c.108C>A, c.108C>G, c.124C>A, c.124C>G, c.144T>A, c.144T>G — all now DROPPED
- c.832G>T (R278L) — Likely_benign LB via W3 (REVEL 0.272, SpliceAI 0)
- c.513T>C, c.366C>T — Benign FAF >5% (ClinVar source, rule-validated)

Please re-spot-check HOXB13 in v3.1 and any other gene you want me to verify.
If you find more discrepancies the same dbNSFP re-check can be run on them.

Full v3.1 change details in `docs/V3_1_FIX.md`.

Paulo

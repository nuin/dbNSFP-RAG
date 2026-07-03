**Subject:** Re: cp_new — ClinVar removed + SpliceAI re-running unmasked

Both concerns acknowledged and acted on.

**1. ClinVar removed entirely.**
v3 upload is at `cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv` —
**13,871 rows**, down from v2's 39,331. Only the three project rules apply:
- W1: FAF >5% → Benign
- W2: synonymous + SpliceAI ≤0.1 + PhastCons <1.0 → Benign
- W3: missense + FAF >0.1% + REVEL <0.29 + SpliceAI ≤0.1 → Likely_benign

Sources kept: synonymous_catalog (W2, 12,253), pipeline (W1/W2/W3, 1,617),
gnomad_common (W1, 1). No ClinVar evidence used anywhere.

As a side effect of the SpliceAI backfill, 603 W2 calls that v2 marked Benign
based on a "0 (assumed)" SpliceAI placeholder turned out to have measured
SpliceAI >0.1 — those are dropped in v3.

Audit trail of the 25,460 rows removed from v2 → v3 is at
`UPLOAD_v3_dropped.tsv` with the reason per row.

**2. SpliceAI mismatch — re-running unmasked.**
Our scores are from a local SpliceAI v1.3.1 run with `-M 1` (Illumina **masked**
mode — zeros out signals at canonical splice sites to highlight altered
activity). The spliceai.org lookup tool uses **unmasked** by default — so it
shows raw scores at all positions including canonical splice. That's the
discrepancy.

Unmasked re-run (`-M 0`) is in progress on the same 24,101 variant sites,
~4 hours. When it completes I'll re-validate W2/W3 thresholds against the
unmasked scores (which is the right value to compare against spliceai.org).

Will send the updated v3 as soon as unmasked SpliceAI finishes. If you have
2–3 spot-check variants you've already looked up on spliceai.org, send them
over and I'll confirm the numbers match before pushing v3 to SeqNext.

Paulo

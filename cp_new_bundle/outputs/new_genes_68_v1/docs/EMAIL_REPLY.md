**Subject:** Re: cp_new pipeline — fixing the synonymous FAF>5% gap

You're right on both counts.

**On the missed >5% synonymous variants** — confirmed root cause. The W1 (FAF>5%) rule today runs *after* a dbNSFP join, and dbNSFP 5.3 only covers non-synonymous variants. So any synonymous common variant has no path through W1, regardless of how high its gnomAD MAF is. The Option A synonymous catalog rows confirm it — every one shows an empty FAF column. The data is missing from our annotation, not from gnomAD.

**On dropping synonymous coding** — that was my filter and I had it backwards. You wanted them in the upload as auto-Benign, not excluded. I'll reverse it.

**Proposal (matches yours):**

1. **Independent gnomAD common-variant pull**, exactly as you suggested. Direct tabix on the 68-gene BED regions against gnomAD v4.1 joint sites, filter to grpmax FAF95 > 5%, mark all Benign. No dbNSFP gate, no transcript gate, no consequence gate. Adds a new `source = gnomad_common` to the merge.
2. **Stop filtering synonymous coding** in the upload builder. Synonymous calls from any source (W2, ClinVar, the new gnomAD pull, the Option A catalog) flow into the upload as Benign.
3. **Backfill FAF onto the Option A catalog** — same gnomAD region pull, joined by coords — so W2 calls also show FAF.

Expected outcome: the upload set grows, but every added row is a population-common variant flagged Benign with explicit FAF evidence — the safest possible auto-call. Should resolve the GOT2-style misses your team flagged earlier.

Will run it and send a fresh upload. If you want to spot-check a specific gene first (GOT2? BRCA1?) before I do the full 68, say the word.

Paulo

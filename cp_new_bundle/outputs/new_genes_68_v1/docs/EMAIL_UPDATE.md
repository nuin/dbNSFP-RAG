**Subject:** cp_new pipeline — conservative upload + per-gene splits

Quick update on the cp_new (68 new genes, hg19) classification work.

**New deliverable: conservative upload set**
- `cp_new_seqnext_UPLOAD_conservative.tsv` — **1,335 rows** (885 Benign + 450 Likely_benign)
- Down from 5,796 in the prior upload. Removes ClinVar 1-star (single-submitter), "Conflicting interpretations", and W1 FAF>5% calls where AlphaMissense=P/LP or REVEL>0.5 or CADD>25 disagree with the benign call.
- Audit trail at `cp_new_seqnext_UPLOAD_conservative_dropped.tsv` — 4,461 rows with `drop_reason`.

**Per-gene files** — now consolidated in `per_gene_all/` with tagged filenames:
- `{GENE}__full.tsv` — all sources, pre-filter (46,344 rows total)
- `{GENE}__pipeline.tsv` — W1/W2/W3 rule hits only (2,189 rows)
- `{GENE}__with_clinvar.tsv` — original upload (5,796 rows)
- `{GENE}__conservative.tsv` — new conservative cut (1,335 rows)

**Coverage gaps to flag** — three genes have zero conservative-cut coverage and need manual review of all variants: **CDKN2C, ELOC, H3-3A**. ELOC has 36 candidates in the full set, all dropped by intronic/coord filters — worth re-investigating the BED interval.

Full QA notes in `UPLOAD_ISSUES.md` (10 items, priority-ordered).

Paulo

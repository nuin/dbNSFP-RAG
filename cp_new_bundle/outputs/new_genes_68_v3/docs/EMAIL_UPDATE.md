**Subject:** cp_new v3 ready — unmasked SpliceAI + ClinVar removed

v3 upload is at `cp_new_bundle/outputs/new_genes_68_v3/UPLOAD_v3.tsv` —
**13,490 rows**, both lab concerns addressed.

**1. SpliceAI matched to spliceai.org**
Re-ran SpliceAI in unmasked mode (`-M 0`) on 24,101 sites. Both columns kept in
the upload so you can audit:
- `SpliceAI_masked` — original `-M 1` Illumina masked
- `SpliceAI_unmasked` — `-M 0`, matches spliceai.org default

The threshold check (W2/W3 require SpliceAI ≤0.1) is now applied against the
**unmasked** value. This dropped 381 rows that had been wrongly called Benign
because masked mode was hiding the canonical splice signal. Examples of what
got correctly removed:

| Gene | HGVSc | masked | unmasked |
|---|---|---:|---:|
| ANKRD26 | c.1816A>C | 0.07 | **0.98** |
| MBD4 | c.252G>A | 0.00 | **0.88** |
| GPR161 | c.-139+1G>C | 0.02 | **0.76** (canonical +1 donor) |

**2. ClinVar removed entirely**
All 24,856 ClinVar-source rows dropped. Only the three project rules apply:

| Rule | Criteria | Call | rows in v3 |
|---|---|---|---:|
| W1 | FAF >5% | Benign | (subset of pipeline + 1 gnomad_common) |
| W2 | synonymous + SpliceAI ≤0.1 + PhastCons <1.0 | Benign | 11,883 |
| W3 | missense + FAF >0.1% + REVEL <0.29 + SpliceAI ≤0.1 | LB | 22 |

Final v3:
- 13,468 Benign + 22 Likely_benign
- 67/68 genes (TERC blank — non-coding RNA, no rules apply)
- Audit trail of dropped rows at `UPLOAD_v3_dropped.tsv`

**Open**
- 11% of the upload has unmasked SpliceAI populated only from the new run; rows
  from the original (May) masked-only VCF don't have an unmasked value yet.
  Can re-run unmasked over the older VCF too if you want 100% coverage on
  both columns.
- If you can share 2–3 spot-check variants you looked up on spliceai.org, I'll
  confirm our `SpliceAI_unmasked` matches before pushing v3 to SeqNext.

Full notes at `docs/SPLICEAI_NOTES.md` and `docs/README.md`.

Paulo

# SpliceAI handling in v3

## Two columns kept

| Column | Mode | Reference |
|---|---|---|
| `SpliceAI_masked` | `-M 1` (Illumina masked) | Zeros out signals at canonical splice sites |
| `SpliceAI_unmasked` | `-M 0` (Illumina raw) | Matches spliceai.org default |

## Why both

The lab compares against spliceai.org, which uses **unmasked** by default. Our
original v2 run used masked mode, which is why scores didn't match.

Both columns are now in the upload so analysts can see both values. The
**unmasked column is what spliceai.org shows**.

## Project-criterion threshold

The W2/W3 rules require SpliceAI ≤0.1. In v3 we apply this threshold against the
**unmasked** value (matching spliceai.org). 381 rows that v3 (initial) marked
Benign using only the masked check were dropped because they have unmasked > 0.1.

## Why masked-vs-unmasked diverge at canonical splice sites

Masked mode zeros the score at known canonical splice sites so it shows only
*altered* splice activity. Unmasked shows the raw model output, including the
canonical signal. So a variant directly at the +1 donor position will read 0 in
masked mode (the canonical signal is masked) but high in unmasked.

Example: `GPR161 c.-139+1G>C` — masked 0.02, **unmasked 0.76** (at +1 donor).

## Coverage gap

`SpliceAI_unmasked` is 88.7% populated vs `SpliceAI_masked` 99.6%. The 11% gap
is positions in the original (May 2026) masked SpliceAI VCF that weren't part of
the 24,101-site set re-scored unmasked. To close the gap, re-run unmasked over
the older VCFs too.

## Re-run command (for the record)

```bash
~/data/spliceai/venv-spliceai/bin/spliceai \
  -I data/exports/cp_new/cp_new_v2_needed.vcf \
  -O data/exports/cp_new/cp_new_v2_needed.spliceai_unmasked.vcf \
  -R ~/data/hg38/hg38.fa \
  -A grch38 \
  -M 0
```

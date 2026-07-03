# v3.1 — bug fixes from lab review (2026-06-16)

Lab reported two issues during HOXB13 spot-check:

## 1. Missense and nonsense variants wrongly tagged "Benign synonymous"

**Examples**:
- HOXB13 c.108C>A — REVEL 0.396, CADD 22.1, AlphaMissense B → **missense, not synonymous**
- HOXB13 c.124C>A — REVEL 0.302, CADD 23.6, AlphaMissense B → **missense**
- HOXB13 c.144T>A — CADD 34.0 → **nonsense**

**Root cause**: `scripts/classify_cp_new.py:is_synonymous()` had this shortcut:

```python
# codon_degeneracy: 2 = synonymous in dbNSFP convention
cd = str(row.get("codon_degeneracy", "")).strip()
if cd in ("2", "4"):
    return True
```

dbNSFP's `codon_degeneracy` describes the **codon position's** degeneracy, not the
specific alt's effect. At a 2-fold degenerate site only 1 of 3 alts is silent;
the other 2 are missense. The buggy code treated all 3 alts as synonymous.

**Fix**: removed the codon_degeneracy shortcut in
`scripts/classify_cp_new.py`. Only `HGVSp ends with '='` or explicit
`aaref == aaalt` count as synonymous now.

**v3.1 patch** (`scripts/build_upload_v3_1.py`) re-checks every pipeline row
labeled "Benign synonymous" against dbNSFP's per-row aaref/aaalt and:
- Truly synonymous → kept as-is
- Missense passing W3 → reclassified `Likely_benign FAF >0.1%, REVEL<0.29, SpliceAI<=0.1`
- Missense failing W3, nonsense (aaalt=X) → dropped from upload

Pass-1 results on the 1,377 pipeline syn-labeled rows in v3:
- 0 truly synonymous
- 3 reclassified to W3 LB
- 1,217 dropped (1,083 missense + 134 nonsense)
- 76 unverified (no dbNSFP AA lookup, no missense predictors → kept)
- 81 dropped via predictor heuristic (had REVEL/AlphaMissense → must be missense)

The 1,377 number reflects how broken the original W2 pipeline was on
non-dbNSFP-AA-tagged variants. The 76 remaining unverified rows are
intronic-ROI / UTR (`c.*XXXX`) rows where there's no protein to compare and
the W2 conservation/splice argument still applies.

## 2. `rule_fails:no_rule_applies` rows removed from upload

Lab: *"If we are not going to import them into SeqPilot, I don't think we
should include them now so please have them filtered out."*

**Fix**: 20,991 `rule_fails:no_rule_applies` rows dropped from v3.1. Other
`rule_fails:X` rows (e.g. `rule_fails:W3_faf_0.001`) are still kept because
they carry useful evidence — the lab can choose to filter them or use them as
a review queue.

## v3 → v3.1

| | v3 | v3.1 |
|---|---:|---:|
| Total rows | 38,314 | **16,025** |
| Rule-based Benign (W1/W2) | 13,557 | 12,258 |
| Rule-based Likely_benign (W3) | 225 | 228 |
| `rule_fails:no_rule_applies` | 20,991 | **0** |
| Other `rule_fails:X` | 3,539 | 3,539 |
| Mis-tagged missense/nonsense | 1,217 | **0** |
| Sources | clinvar+pipeline+syn+gnomad | clinvar+pipeline+syn+gnomad |

## HOXB13 spot-check — fixed examples

| HGVSc | v3 (broken) | v3.1 (correct) |
|---|---|---|
| c.108C>A | Benign synonymous | DROPPED (missense) |
| c.108C>G | Benign synonymous | DROPPED (missense) |
| c.124C>A | Benign synonymous | DROPPED (missense) |
| c.124C>G | Benign synonymous | DROPPED (missense) |
| c.144T>A | Benign synonymous | DROPPED (nonsense) |
| c.144T>G | Benign synonymous | DROPPED (nonsense) |

## Audit trail

All 23,587 dropped rows (mis-tagged + no_rule_applies) appended to
`UPLOAD_v3_dropped.tsv` with reason.

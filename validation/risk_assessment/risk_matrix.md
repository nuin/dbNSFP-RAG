# Risk Classification Matrix

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Version** | 1.0 |
| **Date** | 2026-01-11 |
| **Standard** | ISO 14971:2019 |

---

## 1. Risk Matrix Visualization

### Severity vs Probability Matrix

```
                           PROBABILITY
              ┌─────────┬─────────┬─────────┬─────────┬─────────┐
              │Improbable│ Remote  │Occasional│Probable │Frequent │
              │   (1)   │   (2)   │   (3)   │   (4)   │   (5)   │
    ┌─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
    │Catastro-│         │         │         │         │         │
S   │phic (5) │   LOW   │ MEDIUM  │  HIGH   │  HIGH   │  HIGH   │
E   ├─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
V   │Major    │         │         │         │         │         │
E   │  (4)    │   LOW   │   LOW   │ MEDIUM  │  HIGH   │  HIGH   │
R   ├─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
I   │Moderate │         │         │         │         │         │
T   │  (3)    │   LOW   │   LOW   │   LOW   │ MEDIUM  │ MEDIUM  │
Y   ├─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
    │Minor    │         │         │         │         │         │
    │  (2)    │   LOW   │   LOW   │   LOW   │   LOW   │ MEDIUM  │
    ├─────────┼─────────┼─────────┼─────────┼─────────┼─────────┤
    │Negligib-│         │         │         │         │         │
    │le (1)   │   LOW   │   LOW   │   LOW   │   LOW   │   LOW   │
    └─────────┴─────────┴─────────┴─────────┴─────────┴─────────┘
```

### Risk Level Definitions

| Level | Color | Action Required |
|-------|-------|-----------------|
| **HIGH** | Red | Unacceptable - Must implement risk controls before release |
| **MEDIUM** | Yellow | Conditionally acceptable - Requires documented justification |
| **LOW** | Green | Acceptable - No additional controls required |

---

## 2. Failure Mode Distribution

### Before Mitigation

```
HIGH RISK (RPN >= 50)
├── FM-06: Reference allele mismatch      [RPN: 75] ████████████████████
├── FM-01: Invalid variant accepted       [RPN: 64] ████████████████
└── FM-18: GRCh37/38 confusion           [RPN: 60] ███████████████

MEDIUM RISK (RPN 20-49)
├── FM-05: Multi-allelic variant         [RPN: 48] ████████████
├── FM-12: Wrong confidence score        [RPN: 48] ████████████
├── FM-14: Stale ClinVar data           [RPN: 48] ████████████
├── FM-07: Pathogenic -> Benign         [RPN: 45] ███████████
├── FM-08: Benign -> Pathogenic         [RPN: 36] █████████
├── FM-09: VUS over-called              [RPN: 36] █████████
├── FM-24: Embedding mismatch           [RPN: 32] ████████
├── FM-13: Database corruption          [RPN: 30] ███████
├── FM-03: Negative/zero position       [RPN: 27] ██████
├── FM-17: Score parsing error          [RPN: 27] ██████
├── FM-26: Data exposure                [RPN: 24] ██████
└── FM-20: Position off-by-one          [RPN: 24] ██████

LOW RISK (RPN < 20)
├── FM-10: Invalid ACMG criteria        [RPN: 18]
├── FM-15: Missing variants             [RPN: 18]
├── FM-19: Liftover failure             [RPN: 18]
├── FM-21: Model timeout                [RPN: 18]
├── FM-25: Prompt injection             [RPN: 18]
├── FM-27: Unauthorized access          [RPN: 18]
├── FM-02: Malformed chromosome         [RPN: 16]
├── FM-04: Non-ACGT allele             [RPN: 12]
├── FM-11: Missing criteria             [RPN: 12]
├── FM-16: Duplicate variants           [RPN: 8]
├── FM-22: Model crash                  [RPN: 6]
└── FM-23: API unavailable              [RPN: 6]
```

### After Mitigation (Target State)

```
HIGH RISK (RPN >= 50)
└── (None - all mitigated)

MEDIUM RISK (RPN 20-49)
├── FM-07: Pathogenic -> Benign         [RPN: 30] (Accepted)
├── FM-08: Benign -> Pathogenic         [RPN: 24] (Accepted)
├── FM-09: VUS over-called              [RPN: 24] (Accepted)
├── FM-14: Stale ClinVar data           [RPN: 24] (Mitigated)
├── FM-18: GRCh37/38 confusion          [RPN: 20] (Mitigated)
├── FM-06: Reference allele mismatch    [RPN: 20] (Mitigated)
└── ... (others similar or lower)

LOW RISK (RPN < 20)
└── (Remaining 20+ failure modes)
```

---

## 3. Risk Acceptability Criteria

### Acceptability Matrix

| Severity | Acceptable Probability |
|----------|----------------------|
| Catastrophic (5) | Remote (2) or lower after controls |
| Major (4) | Occasional (3) or lower after controls |
| Moderate (3) | Probable (4) or lower |
| Minor (2) | Any |
| Negligible (1) | Any |

### Risk/Benefit Analysis

For medium-risk items that cannot be fully mitigated:

| Failure Mode | Residual Risk | Benefit | Accept? |
|--------------|---------------|---------|---------|
| FM-07: P->B misclass | <10% of pathogenic variants | Provides guidance for 90%+ | Yes |
| FM-08: B->P misclass | <10% of benign variants | Flags variants for review | Yes |
| FM-09: VUS over-call | ~30% of classifications | Conservative = safe | Yes |

---

## 4. Control Effectiveness Requirements

| Control Type | Required Effectiveness |
|--------------|----------------------|
| Prevention | Must reduce P by >= 2 levels |
| Detection | Must reduce D by >= 2 levels |
| Alert/Warning | Must reduce D by >= 1 level |
| Guidance/Education | Supplementary only |

---

## 5. Traceability

| Failure Mode | FMEA Section | Test Case | Validation Protocol |
|--------------|--------------|-----------|-------------------|
| FM-01 | 3.1 | test_input_validation.py | OQ-001 |
| FM-06 | 3.1 | test_input_validation.py | OQ-002 |
| FM-07 | 3.2 | test_classification.py | PQ-001 |
| FM-18 | 3.4 | test_api.py | OQ-003 |

---

## 6. Document History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-01-11 | Initial matrix |

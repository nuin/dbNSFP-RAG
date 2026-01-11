# Validation Report

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Version** | [VERSION] |
| **Report Date** | [DATE] |
| **Prepared By** | [NAME] |
| **Reviewed By** | [NAME] |

---

## 1. Executive Summary

### 1.1 Purpose
This report documents the validation activities performed on the ACMG Variant Classification System to demonstrate that it meets its intended use requirements.

### 1.2 Conclusion
[PASS/FAIL with summary]

### 1.3 Key Findings

| Category | Result |
|----------|--------|
| Installation Qualification | [PASS/FAIL] |
| Operational Qualification | [PASS/FAIL] |
| Performance Qualification | [PASS/FAIL] |
| Deviations | [COUNT] |

---

## 2. Test Environment

### 2.1 Hardware

| Component | Specification |
|-----------|--------------|
| System | [MODEL] |
| CPU | [CPU] |
| RAM | [RAM] |
| Storage | [STORAGE] |
| GPU | [GPU or N/A] |

### 2.2 Software

| Component | Version |
|-----------|---------|
| Operating System | [OS] |
| Python | [VERSION] |
| FAISS | [VERSION] |
| sentence-transformers | [VERSION] |
| Model | [VERSION/DATE] |
| Database | [BUILD DATE] |

---

## 3. Installation Qualification Results

| Test ID | Description | Expected | Actual | Result |
|---------|-------------|----------|--------|--------|
| IQ-001 | Python version | >= 3.10 | | |
| IQ-002 | Dependencies | All installed | | |
| IQ-003 | Model files | Present | | |
| IQ-004 | Database files | Present | | |
| IQ-005 | Configuration | Loads | | |
| IQ-006 | FAISS import | Success | | |
| IQ-007 | Embeddings | Load | | |
| IQ-008 | API server | Responds | | |

---

## 4. Operational Qualification Results

### 4.1 Input Validation (OQ-001)

| Test | Input | Expected | Actual | Result |
|------|-------|----------|--------|--------|
| Valid SNV | chr=17, pos=41197801 | Accepted | | |
| Invalid chr | chr=99 | Rejected | | |
| Negative pos | pos=-100 | Rejected | | |
| Invalid allele | ref=X | Rejected | | |
| Multi-allelic | alt=G,T | Rejected | | |

**Overall Result**: [PASS/FAIL]

### 4.2 API Endpoints (OQ-002)

| Endpoint | Method | Expected | Actual | Result |
|----------|--------|----------|--------|--------|
| /health | GET | 200 | | |
| /lookup | GET | 200/404 | | |
| /classify | GET | Response | | |
| /classify | POST | Response | | |
| /gene/{symbol} | GET | List | | |

**Overall Result**: [PASS/FAIL]

### 4.3 Database Operations (OQ-003)

| Operation | Expected | Actual | Result |
|-----------|----------|--------|--------|
| get_by_id() | Returns data | | |
| get_by_gene() | Returns list | | |
| search() | Returns ranked | | |

**Overall Result**: [PASS/FAIL]

---

## 5. Performance Qualification Results

### 5.1 Classification Accuracy (PQ-001)

| Metric | Threshold | Actual | Result |
|--------|-----------|--------|--------|
| Pathogenic sensitivity | >= 90% | | |
| Benign specificity | >= 90% | | |
| Catastrophic errors | 0% | | |

**Gold Standard Dataset**: [N] variants from ClinVar

**Detailed Results**:

| Classification | Total | Correct | Accuracy |
|---------------|-------|---------|----------|
| Pathogenic | | | |
| Likely Pathogenic | | | |
| Benign | | | |
| Likely Benign | | | |
| VUS | | | |

**Overall Result**: [PASS/FAIL]

### 5.2 Performance Metrics (PQ-002)

| Metric | Threshold | Actual | Result |
|--------|-----------|--------|--------|
| Classification latency (p95) | < 2000ms | | |
| Lookup latency (p95) | < 500ms | | |
| Health check latency | < 100ms | | |
| Model load time | < 60s | | |

**Overall Result**: [PASS/FAIL]

### 5.3 Consistency (PQ-003)

| Test | Expected | Actual | Result |
|------|----------|--------|--------|
| Repeated classification | Identical | | |
| Input normalization | Identical | | |

**Overall Result**: [PASS/FAIL]

---

## 6. Deviations

### Deviation 1: [TITLE]

| Field | Value |
|-------|-------|
| Test ID | |
| Description | |
| Expected | |
| Actual | |
| Root Cause | |
| Impact Assessment | |
| Resolution | |
| Retest Result | |

---

## 7. Risk Assessment Update

Based on validation results, the following risk assessments are updated:

| FM ID | Original RPN | Post-Validation RPN | Status |
|-------|-------------|---------------------|--------|
| FM-01 | 64 | | |
| FM-06 | 75 | | |
| FM-07 | 45 | | |
| FM-18 | 60 | | |

---

## 8. Conclusion

### 8.1 Summary

[Summary of validation activities and results]

### 8.2 Recommendation

[Recommendation: Release / Conditional Release / Do Not Release]

### 8.3 Conditions (if applicable)

[List any conditions for release]

---

## 9. Approval

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Test Executor | | | |
| Quality Assurance | | | |
| Software Lead | | | |
| Clinical Advisor | | | |
| Regulatory Affairs | | | |

---

## Appendix A: Test Output Logs

[Attach pytest output]

```
[PASTE TEST OUTPUT HERE]
```

---

## Appendix B: Coverage Report

[Attach coverage report summary]

---

## Appendix C: Reference Documents

| Document | Version | Location |
|----------|---------|----------|
| FMEA | 1.0 | validation/risk_assessment/fmea.md |
| Risk Matrix | 1.0 | validation/risk_assessment/risk_matrix.md |
| Validation Protocol | 1.0 | validation/protocols/validation_protocol.md |
| IEC 62304 Checklist | 1.0 | validation/compliance/iec62304_checklist.md |
| ISO 14971 Checklist | 1.0 | validation/compliance/iso14971_checklist.md |

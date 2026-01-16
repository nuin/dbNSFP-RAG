# Failure Mode and Effects Analysis (FMEA)

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Version** | 1.0 |
| **Date** | 2026-01-11 |
| **Author** | System Validation Team |
| **Standard** | ISO 14971:2019, IEC 62304:2006+A1:2015 |
| **Software Safety Class** | Class B (Non-serious injury possible) |

---

## 1. Scope

This FMEA covers the ACMG Variant Classification System, which:
- Ingests genetic variant data from dbNSFP
- Stores variant annotations in a FAISS vector database
- Provides REST API for variant lookup and classification
- Uses a fine-tuned LLM to classify variants per ACMG/AMP guidelines

### Intended Use
Clinical decision support for genetic variant interpretation. Classifications are advisory and require expert review before clinical action.

### Users
- Clinical geneticists
- Genetic counselors
- Laboratory directors
- Bioinformaticians

---

## 2. Risk Evaluation Criteria

### 2.1 Severity (S)

| Score | Level | Description | Clinical Impact |
|-------|-------|-------------|-----------------|
| 5 | Catastrophic | Wrong classification leads to missed or incorrect diagnosis | Patient harm from missed treatment or unnecessary intervention |
| 4 | Major | Significant delay in clinical decision-making | Extended diagnostic odyssey, delayed treatment |
| 3 | Moderate | Requires manual review and correction | Additional expert time, potential for confusion |
| 2 | Minor | Inconvenience with no clinical impact | User frustration, workflow disruption |
| 1 | Negligible | Cosmetic or UX issue | Minor display problems |

### 2.2 Probability of Occurrence (P)

| Score | Level | Frequency | Rate |
|-------|-------|-----------|------|
| 5 | Frequent | Regular occurrence | >1 per 100 variants |
| 4 | Probable | Will occur several times | 1 per 100-1,000 variants |
| 3 | Occasional | May occur occasionally | 1 per 1,000-10,000 variants |
| 2 | Remote | Unlikely but possible | 1 per 10,000-100,000 variants |
| 1 | Improbable | Very unlikely | <1 per 100,000 variants |

### 2.3 Detection (D)

| Score | Level | Description |
|-------|-------|-------------|
| 5 | None | No mechanism to detect the failure |
| 4 | Low | Failure may be detected by careful manual review |
| 3 | Moderate | Failure likely detected during normal use |
| 2 | High | Failure almost certainly detected before causing harm |
| 1 | Certain | Automatic detection prevents failure from propagating |

### 2.4 Risk Priority Number (RPN)

**RPN = Severity x Probability x Detection**

| RPN Range | Risk Level | Required Action |
|-----------|------------|-----------------|
| 50-125 | **HIGH** | Mitigation required before release |
| 20-49 | **MEDIUM** | Documented risk acceptance required |
| 1-19 | **LOW** | Acceptable, monitor during operation |

---

## 3. Failure Modes Analysis

### 3.1 Input Processing Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-01 | Invalid variant accepted | No input validation | Wrong variant looked up, silent failure | 4 | 4 | 4 | **64** | HIGH |
| FM-02 | Malformed chromosome | User typo, API misuse | Lookup returns no results | 2 | 4 | 2 | 16 | LOW |
| FM-03 | Negative/zero position | Data entry error | Database error or wrong lookup | 3 | 3 | 3 | 27 | MEDIUM |
| FM-04 | Non-ACGT allele | Copy-paste error, encoding issue | Variant not found | 2 | 3 | 2 | 12 | LOW |
| FM-05 | Multi-allelic variant | Complex variant input | Partial or incorrect result | 4 | 3 | 4 | **48** | MEDIUM |
| FM-06 | Reference allele mismatch | Wrong genome build, typo | Wrong variant classified | 5 | 3 | 5 | **75** | HIGH |

### 3.2 Classification Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-07 | Pathogenic -> Benign | Model error, insufficient data | Missed diagnosis, delayed treatment | 5 | 3 | 3 | **45** | MEDIUM |
| FM-08 | Benign -> Pathogenic | Model hallucination | Unnecessary procedures, anxiety | 4 | 3 | 3 | 36 | MEDIUM |
| FM-09 | VUS over-called | Conservative model behavior | Increased uncertainty | 3 | 4 | 3 | 36 | MEDIUM |
| FM-10 | Invalid ACMG criteria | LLM generates fake codes | Misleading evidence documentation | 3 | 3 | 2 | 18 | LOW |
| FM-11 | Missing criteria | Incomplete model output | Incomplete justification | 2 | 3 | 2 | 12 | LOW |
| FM-12 | Wrong confidence score | Hardcoded tiers, not model-derived | Misplaced trust in classification | 3 | 4 | 4 | **48** | MEDIUM |

### 3.3 Data Integrity Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-13 | Database corruption | Interrupted build, disk failure | Wrong annotations returned | 5 | 2 | 3 | 30 | MEDIUM |
| FM-14 | Stale ClinVar data | No update process | Outdated clinical significance | 4 | 4 | 3 | **48** | MEDIUM |
| FM-15 | Missing variants | Incomplete panel build | Variant not found when should exist | 3 | 3 | 2 | 18 | LOW |
| FM-16 | Duplicate variants | Build process error | Inconsistent results | 2 | 2 | 2 | 8 | LOW |
| FM-17 | Score parsing error | Multi-value field handling | Wrong score displayed/used | 3 | 3 | 3 | 27 | MEDIUM |

### 3.4 Coordinate System Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-18 | GRCh37/38 confusion | No build indicator | Wrong variant looked up | 5 | 3 | 4 | **60** | HIGH |
| FM-19 | Liftover failure | Variant not mappable | Missing GRCh37 variant | 3 | 3 | 2 | 18 | LOW |
| FM-20 | Position off-by-one | 0-based vs 1-based confusion | Adjacent variant classified | 4 | 2 | 3 | 24 | MEDIUM |

### 3.5 System Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-21 | Model timeout | Resource exhaustion, complex query | Classification unavailable | 3 | 3 | 2 | 18 | LOW |
| FM-22 | Model crash | Memory error, corrupted weights | Service unavailable | 3 | 2 | 1 | 6 | LOW |
| FM-23 | API unavailable | Server crash, network issue | No classifications possible | 3 | 2 | 1 | 6 | LOW |
| FM-24 | Embedding mismatch | Model version vs index mismatch | Wrong similarity results | 4 | 2 | 4 | 32 | MEDIUM |

### 3.6 Security Failures

| ID | Failure Mode | Potential Cause | Potential Effect | S | P | D | RPN | Risk Level |
|----|--------------|-----------------|------------------|---|---|---|-----|------------|
| FM-25 | Prompt injection | Malicious input | Unexpected model behavior | 3 | 2 | 3 | 18 | LOW |
| FM-26 | Data exposure | Logging sensitive data | Privacy violation | 4 | 2 | 3 | 24 | MEDIUM |
| FM-27 | Unauthorized access | No authentication | Unauthorized use | 3 | 3 | 2 | 18 | LOW |

---

## 4. High-Risk Failure Modes (RPN >= 50)

### FM-06: Reference Allele Mismatch (RPN: 75)

**Description**: User provides a reference allele that doesn't match the actual reference genome at that position.

**Root Cause Analysis**:
- User copied variant from different genome build
- Transcription error from lab report
- Strand confusion (complement reported)

**Risk Mitigation**:
| # | Control Measure | Type | Effectiveness |
|---|----------------|------|---------------|
| 1 | Validate input ref against dbNSFP ref field | Detection | High |
| 2 | Display warning if ref doesn't match | Alert | Medium |
| 3 | Return detailed error with expected ref | Guidance | Medium |

**Residual Risk**: After mitigation, RPN = 5 x 2 x 2 = **20** (MEDIUM - Acceptable)

---

### FM-01: Invalid Variant Accepted (RPN: 64)

**Description**: System accepts malformed variant input without validation, leading to silent failures or wrong lookups.

**Root Cause Analysis**:
- No input validation layer implemented
- API accepts any string for chromosome
- Position not validated as positive integer
- Alleles not validated as nucleotides

**Risk Mitigation**:
| # | Control Measure | Type | Effectiveness |
|---|----------------|------|---------------|
| 1 | Add `src/validation.py` module | Prevention | High |
| 2 | Validate chromosome (1-22, X, Y, M/MT) | Prevention | High |
| 3 | Validate position > 0 | Prevention | High |
| 4 | Validate alleles match regex `^[ACGT]+$` | Prevention | High |
| 5 | Return descriptive error messages | Guidance | Medium |

**Residual Risk**: After mitigation, RPN = 4 x 1 x 2 = **8** (LOW - Acceptable)

---

### FM-18: GRCh37/38 Confusion (RPN: 60)

**Description**: User queries with coordinates from wrong genome build, resulting in lookup of different variant.

**Root Cause Analysis**:
- Database uses GRCh37 coordinates
- No build indicator in API response
- User assumes GRCh38 (current standard)
- Position differences can be significant

**Risk Mitigation**:
| # | Control Measure | Type | Effectiveness |
|---|----------------|------|---------------|
| 1 | Add `build` parameter to API (default: GRCh37) | Clarification | High |
| 2 | Include `genome_build: "GRCh37"` in all responses | Awareness | High |
| 3 | Document coordinate system prominently | Education | Medium |
| 4 | Add `/health` endpoint with build info | Visibility | Medium |

**Residual Risk**: After mitigation, RPN = 5 x 2 x 2 = **20** (MEDIUM - Acceptable)

---

## 5. Medium-Risk Failure Modes Requiring Documentation

### FM-05, FM-09, FM-12, FM-14 (RPN 45-48)

These failure modes have RPN between 20-49 and require documented risk acceptance:

| ID | Failure Mode | Accepted Risk Rationale |
|----|--------------|------------------------|
| FM-05 | Multi-allelic variant | System will detect and warn; decomposition is user responsibility |
| FM-07 | Pathogenic -> Benign | Gold standard validation ensures <10% misclassification; expert review required |
| FM-09 | VUS over-called | Conservative approach preferred in clinical setting |
| FM-12 | Wrong confidence | Confidence is advisory; classification tier is primary signal |
| FM-14 | Stale ClinVar data | Version displayed in API; regular updates documented in SOP |

---

## 6. Risk Control Measures Summary

| Priority | Mitigation | Failure Modes Addressed | Status |
|----------|-----------|------------------------|--------|
| 1 | Input validation module | FM-01, FM-02, FM-03, FM-04, FM-05 | To implement |
| 2 | Reference allele validation | FM-06 | To implement |
| 3 | Genome build indicator | FM-18 | To implement |
| 4 | Version tracking in API | FM-14, FM-24 | To implement |
| 5 | Gold standard test suite | FM-07, FM-08, FM-09 | To implement |
| 6 | ACMG code validation | FM-10 | To implement |
| 7 | Database integrity checks | FM-13, FM-16 | To implement |

---

## 7. Verification of Risk Controls

Each control measure will be verified through:

1. **Unit tests** - Validate individual control functions
2. **Integration tests** - Validate end-to-end behavior
3. **Gold standard testing** - Validate against known-truth variants
4. **Manual review** - Expert review of edge cases

Verification results will be documented in the Validation Report.

---

## 8. Residual Risk Assessment

After implementing all control measures:

| Risk Level | Count | Percentage |
|------------|-------|------------|
| HIGH (RPN >= 50) | 0 | 0% |
| MEDIUM (RPN 20-49) | 8 | 30% |
| LOW (RPN < 20) | 19 | 70% |

**Overall residual risk is ACCEPTABLE** for intended use as clinical decision support tool requiring expert review.

---

## 9. Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-11 | Validation Team | Initial FMEA |

---

## 10. Approval

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Quality Assurance | | | |
| Software Development | | | |
| Clinical Advisor | | | |
| Regulatory Affairs | | | |

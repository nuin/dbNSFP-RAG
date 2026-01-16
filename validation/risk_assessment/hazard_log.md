# Hazard Log

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Version** | 1.0 |
| **Date** | 2026-01-11 |
| **Standard** | ISO 14971:2019 |

---

## Hazard Tracking

| Hazard ID | Description | Source | Status | Owner | Due Date |
|-----------|-------------|--------|--------|-------|----------|
| HAZ-001 | Invalid variant input leads to wrong classification | FMEA FM-01 | Open | Dev Team | - |
| HAZ-002 | Reference allele mismatch causes silent wrong lookup | FMEA FM-06 | Open | Dev Team | - |
| HAZ-003 | Genome build confusion (GRCh37 vs GRCh38) | FMEA FM-18 | Open | Dev Team | - |
| HAZ-004 | Pathogenic variant misclassified as Benign | FMEA FM-07 | Open | ML Team | - |
| HAZ-005 | Stale ClinVar data leads to outdated classifications | FMEA FM-14 | Open | Data Team | - |
| HAZ-006 | Multi-allelic variants not properly handled | FMEA FM-05 | Open | Dev Team | - |

---

## Hazard Details

### HAZ-001: Invalid Variant Input

**Description**: System accepts malformed variant coordinates without validation

**Harm**: Wrong variant looked up, leading to incorrect classification

**Sequence of Events**:
1. User submits variant with invalid chromosome (e.g., "chr99")
2. System does not validate input
3. Lookup fails silently or returns wrong variant
4. User receives incorrect or no classification

**Risk Control Measures**:
- [ ] Implement input validation module (`src/validation.py`)
- [ ] Validate chromosome against allowed values
- [ ] Validate position is positive integer
- [ ] Validate alleles are valid nucleotides
- [ ] Return descriptive error messages

**Verification**: test_input_validation.py

---

### HAZ-002: Reference Allele Mismatch

**Description**: User provides reference allele that doesn't match genome

**Harm**: Different variant classified than intended

**Sequence of Events**:
1. User copies variant from report using different genome build
2. Reference allele at GRCh37 position differs from input
3. System looks up variant with mismatched ref
4. Wrong variant or no variant returned

**Risk Control Measures**:
- [ ] Compare input ref with database ref field
- [ ] Warn if ref doesn't match
- [ ] Include expected ref in error response

**Verification**: test_input_validation.py::test_ref_allele_validation

---

### HAZ-003: Genome Build Confusion

**Description**: User queries with coordinates from wrong genome build

**Harm**: Completely different genomic position queried

**Sequence of Events**:
1. Database built with GRCh37 coordinates
2. User assumes GRCh38 (current standard)
3. Position difference may be hundreds to thousands of bases
4. Different or no variant returned

**Risk Control Measures**:
- [ ] Add `build` parameter to API endpoints
- [ ] Include `genome_build: "GRCh37"` in all responses
- [ ] Document coordinate system in API docs
- [ ] Add build info to `/health` endpoint

**Verification**: test_api.py::test_genome_build_indicator

---

### HAZ-004: Pathogenic Misclassification

**Description**: Pathogenic variant incorrectly classified as Benign

**Harm**: Missed diagnosis, delayed or absent treatment

**Sequence of Events**:
1. Pathogenic variant submitted for classification
2. Model lacks sufficient training data or evidence
3. Classification returned as Benign
4. Clinician may not pursue further testing

**Risk Control Measures**:
- [ ] Gold standard validation with known variants
- [ ] Require >= 90% sensitivity on pathogenic variants
- [ ] Include confidence score and criteria
- [ ] Document that expert review is required

**Verification**: test_classification.py::test_pathogenic_sensitivity

---

### HAZ-005: Stale ClinVar Data

**Description**: Database contains outdated ClinVar annotations

**Harm**: Classifications based on superseded clinical evidence

**Sequence of Events**:
1. dbNSFP version includes ClinVar from specific date
2. ClinVar updates classifications over time
3. System returns outdated significance
4. Classification may not reflect current evidence

**Risk Control Measures**:
- [ ] Track dbNSFP version in database metadata
- [ ] Display version in `/health` endpoint
- [ ] Display version in classification response
- [ ] Document update schedule in SOP

**Verification**: test_api.py::test_version_tracking

---

### HAZ-006: Multi-allelic Variants

**Description**: System cannot properly handle multi-allelic sites

**Harm**: Partial or incorrect variant classification

**Sequence of Events**:
1. User submits variant with multiple alt alleles (e.g., "A" -> "G,T")
2. System not designed for multi-allelic representation
3. May classify only first allele or fail
4. User receives incomplete result

**Risk Control Measures**:
- [ ] Detect multi-allelic input (comma in alt)
- [ ] Return error with guidance to decompose
- [ ] Document multi-allelic handling in API docs

**Verification**: test_input_validation.py::test_multi_allelic_detection

---

## Status Definitions

| Status | Description |
|--------|-------------|
| Open | Hazard identified, controls not implemented |
| In Progress | Controls being implemented |
| Verification | Controls implemented, awaiting verification |
| Closed | Controls verified effective |
| Accepted | Risk accepted without additional controls |

---

## Document History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-01-11 | Initial hazard log |

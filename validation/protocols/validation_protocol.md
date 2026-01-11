# Validation Protocol

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Version** | 1.0 |
| **Date** | 2026-01-11 |
| **Standard** | IEC 62304:2006+A1:2015, ISO 13485:2016 |

---

## 1. Purpose

This protocol defines the validation approach for the ACMG Variant Classification System to ensure it meets its intended use requirements and performs safely and effectively.

---

## 2. Scope

### In Scope
- Software functionality validation
- Classification accuracy validation
- API interface validation
- Data integrity validation
- Performance validation

### Out of Scope
- Hardware validation
- Network infrastructure
- User training effectiveness
- Clinical utility studies

---

## 3. Validation Phases

### 3.1 Installation Qualification (IQ)

**Objective**: Verify the system is installed correctly with all dependencies.

| Test ID | Description | Method | Acceptance Criteria |
|---------|-------------|--------|---------------------|
| IQ-001 | Python version | `python --version` | Python 3.10+ |
| IQ-002 | Dependencies installed | `pip list` | All requirements.txt packages present |
| IQ-003 | Model files present | File check | model/ directory contains weights |
| IQ-004 | Database files present | File check | vectordb/ contains index files |
| IQ-005 | Configuration valid | Import test | config.py loads without error |
| IQ-006 | FAISS library | Import test | `import faiss` succeeds |
| IQ-007 | Sentence transformers | Import test | Embedding model loads |
| IQ-008 | API server starts | `uvicorn` start | Server responds to /health |

**Execution**:
```bash
# IQ-001: Python version
python --version

# IQ-002: Dependencies
pip check

# IQ-003-004: File presence
ls -la models/acmg-classifier/model/
ls -la data/vectordb/grch37-ngsgenes/

# IQ-005-007: Import tests
python -c "from src import config, vectorstore, chunker"
python -c "import faiss"
python -c "from sentence_transformers import SentenceTransformer"

# IQ-008: Server start
uvicorn api.server:app --host 0.0.0.0 --port 8000 &
curl http://localhost:8000/health
```

---

### 3.2 Operational Qualification (OQ)

**Objective**: Verify all features function correctly under normal operating conditions.

#### OQ-001: Input Validation

| Test | Input | Expected Result |
|------|-------|-----------------|
| Valid SNV | chr=17, pos=41197801, ref=T, alt=A | Accepted |
| Invalid chromosome | chr=99 | Rejected with error |
| Negative position | pos=-100 | Rejected with error |
| Invalid allele | ref=X | Rejected with error |
| Multi-allelic | alt=G,T | Rejected with guidance |

**Execution**:
```bash
pytest validation/test_suite/test_input_validation.py -v
```

#### OQ-002: API Endpoints

| Endpoint | Method | Test | Expected |
|----------|--------|------|----------|
| /health | GET | Basic request | 200 OK |
| /lookup | GET | Valid variant | 200 or 404 |
| /classify | GET | Valid variant | Classification response |
| /classify | POST | JSON body | Classification response |
| /gene/{symbol} | GET | Valid gene | Variant list |

**Execution**:
```bash
pytest validation/test_suite/test_api.py -v
```

#### OQ-003: Database Operations

| Test | Operation | Expected |
|------|-----------|----------|
| Variant lookup | get_by_id() | Returns variant or None |
| Gene query | get_by_gene() | Returns variant list |
| Semantic search | search() | Returns ranked results |

**Execution**:
```bash
pytest validation/test_suite/test_vectorstore.py -v
```

#### OQ-004: Genome Build Handling

| Test | Scenario | Expected |
|------|----------|----------|
| GRCh37 coordinates | Default database | Correct lookup |
| Build indicator | API response | genome_build field present |
| Version tracking | /health endpoint | Version info returned |

---

### 3.3 Performance Qualification (PQ)

**Objective**: Verify the system performs accurately on representative clinical data.

#### PQ-001: Classification Accuracy

**Gold Standard Dataset**: 35 ClinVar variants with expert-reviewed classifications

| Metric | Threshold | Method |
|--------|-----------|--------|
| Pathogenic sensitivity | >= 90% | P/LP variants classified as P or LP |
| Benign specificity | >= 90% | B/LB variants classified as B or LB |
| Catastrophic errors | 0% | P never classified as B, vice versa |

**Execution**:
```bash
pytest validation/test_suite/test_classification.py -v -m "not slow"
```

#### PQ-002: Performance Metrics

| Metric | Threshold |
|--------|-----------|
| Single classification latency (p95) | < 2000ms |
| Lookup latency (p95) | < 500ms |
| Health check latency | < 100ms |
| Model load time | < 60s |

**Execution**:
```bash
pytest validation/test_suite/ -v -m slow
```

#### PQ-003: Consistency

| Test | Method | Threshold |
|------|--------|-----------|
| Same variant consistency | 3 repeated classifications | 100% identical |
| Input normalization | chr17 vs 17 | Same result |

---

## 4. Test Environment

### 4.1 Hardware Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 20 GB | 50+ GB |
| GPU | None (CPU inference) | Apple M1+ or NVIDIA |

### 4.2 Software Requirements

| Component | Version |
|-----------|---------|
| Python | 3.10+ |
| Operating System | macOS 12+, Ubuntu 20.04+, Windows 10+ |
| Database | FAISS with grch37-ngsgenes index |
| Model | acmg-classifier LoRA weights |

---

## 5. Acceptance Criteria Summary

| Phase | Criteria | Pass/Fail |
|-------|----------|-----------|
| IQ | All components installed | |
| OQ-Input | 100% invalid inputs rejected | |
| OQ-API | All endpoints respond correctly | |
| OQ-DB | All operations work | |
| PQ-Accuracy | Pathogenic sensitivity >= 90% | |
| PQ-Accuracy | Benign specificity >= 90% | |
| PQ-Accuracy | 0 catastrophic errors | |
| PQ-Performance | Latency within thresholds | |
| PQ-Consistency | 100% reproducible | |

---

## 6. Deviation Handling

Any test failure or deviation from expected results shall be:

1. **Documented** - Record exact failure details
2. **Investigated** - Determine root cause
3. **Assessed** - Evaluate impact on safety/effectiveness
4. **Resolved** - Fix issue or document risk acceptance
5. **Retested** - Verify resolution

---

## 7. Documentation Requirements

### 7.1 Test Records

Each test execution shall record:
- Test ID and description
- Date and time of execution
- Tester name/ID
- Software version tested
- Test environment details
- Actual results
- Pass/Fail determination
- Any deviations or observations

### 7.2 Validation Report

The validation report shall include:
- Executive summary
- Test environment description
- IQ/OQ/PQ results summary
- Deviation summary
- Conclusion statement
- Approval signatures

---

## 8. Revalidation Triggers

Revalidation is required when:

| Change | Revalidation Scope |
|--------|-------------------|
| Model retrained | PQ-001, PQ-003 |
| Database rebuilt | OQ-003, PQ-001 |
| API changes | OQ-002, OQ-004 |
| Major code changes | Full OQ + PQ |
| Dependencies updated | IQ + regression tests |
| dbNSFP version update | PQ-001 (accuracy check) |

---

## 9. Approval

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Quality Assurance | | | |
| Software Lead | | | |
| Clinical Advisor | | | |
| Regulatory Affairs | | | |

---

## Appendix A: Test Commands

```bash
# Full validation suite
pytest validation/test_suite/ -v --tb=short

# IQ only
pytest validation/test_suite/ -v -m "not requires_model and not requires_database"

# OQ only
pytest validation/test_suite/test_input_validation.py validation/test_suite/test_api.py -v

# PQ only
pytest validation/test_suite/test_classification.py -v

# With coverage
pytest validation/test_suite/ --cov=src --cov=api --cov-report=html

# Generate JUnit XML for CI
pytest validation/test_suite/ --junitxml=validation/reports/results.xml
```

---

## Appendix B: Traceability Matrix

| Requirement | Test ID | FMEA Reference |
|-------------|---------|----------------|
| Input validation | OQ-001 | FM-01, FM-02, FM-03, FM-04 |
| Reference validation | OQ-001 | FM-06 |
| Multi-allelic detection | OQ-001 | FM-05 |
| Build indication | OQ-004 | FM-18 |
| Classification accuracy | PQ-001 | FM-07, FM-08 |
| Version tracking | OQ-004 | FM-14 |
| Performance | PQ-002 | FM-21 |

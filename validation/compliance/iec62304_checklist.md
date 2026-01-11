# IEC 62304 Compliance Checklist

## Document Information

| Field | Value |
|-------|-------|
| **Product** | ACMG Variant Classification System |
| **Standard** | IEC 62304:2006+A1:2015 |
| **Safety Class** | Class B |
| **Date** | 2026-01-11 |

---

## Software Safety Classification

### Classification Rationale

| Factor | Assessment |
|--------|------------|
| **Can the software cause or contribute to a hazardous situation?** | Yes - incorrect classification could affect clinical decisions |
| **Severity of potential harm** | Non-serious injury (delayed/incorrect diagnosis) |
| **Independent risk control measures** | Yes - expert review required before clinical action |

**Conclusion**: Class B (software that could contribute to hazardous situation resulting in non-serious injury)

---

## Clause 5: Software Development Process

### 5.1 Software Development Planning

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.1.1 Software development plan established | [ ] | CLAUDE.md, this document |
| 5.1.2 Keep plan updated | [ ] | Version history |
| 5.1.3 Reference or include: | | |
| - Development life cycle model | [ ] | Agile with validation gates |
| - Plans for development standards | [ ] | Coding standards (ruff) |
| - Plans for integration | [ ] | CI/CD pipeline |
| - Plans for verification | [ ] | Test suite, code review |
| - Plans for risk management | [ ] | FMEA document |
| - Plans for documentation | [ ] | This checklist, protocols |
| - Plans for configuration management | [ ] | Git, version tags |
| - Plans for problem resolution | [ ] | GitHub issues |

### 5.2 Software Requirements Analysis

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.2.1 Define software requirements | [ ] | Functional specs |
| 5.2.2 Include: | | |
| - Functional requirements | [ ] | API endpoints, classification |
| - Input/output requirements | [ ] | Variant formats, JSON responses |
| - Interfaces to other systems | [ ] | REST API specification |
| - Security requirements | [ ] | Input validation |
| - Usability requirements | [ ] | Error messages |
| - Data definitions | [ ] | Variant, Classification schemas |
| - Installation requirements | [ ] | README, CLAUDE.md |
| - Operating requirements | [ ] | Hardware requirements |
| - Maintenance requirements | [ ] | Update procedures |
| 5.2.3 Re-evaluate medical device risk analysis | [ ] | FMEA updates |
| 5.2.4 Update requirements as needed | [ ] | Change control |
| 5.2.5 Verify requirements | [ ] | Requirements review |
| 5.2.6 Evaluate requirements for risk | [ ] | FMEA traceability |

### 5.3 Software Architectural Design

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.3.1 Document software architecture | [ ] | CLAUDE.md architecture section |
| 5.3.2 Develop architecture implementing requirements | [ ] | Design documents |
| 5.3.3 Document interfaces between software items | [ ] | API documentation |
| 5.3.4 Specify functional/performance requirements of SOUP | [ ] | Dependencies list |
| 5.3.5 Identify segregation for risk control | [ ] | Module separation |
| 5.3.6 Verify architecture | [ ] | Architecture review |

### 5.4 Software Detailed Design

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.4.1 Subdivide into units | [ ] | src/, api/ modules |
| 5.4.2 Develop detailed design for units | [ ] | Code documentation |
| 5.4.3 Develop detailed design for interfaces | [ ] | API specs, docstrings |
| 5.4.4 Verify detailed design | [ ] | Code review |

### 5.5 Software Unit Implementation

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.5.1 Implement software units | [ ] | Source code |
| 5.5.2 Establish unit verification process | [ ] | Unit tests |
| 5.5.3 Unit acceptance criteria | [ ] | Test pass, coverage |
| 5.5.4 Additional acceptance criteria | [ ] | Code review checklist |
| 5.5.5 Verify software unit | [ ] | pytest results |

### 5.6 Software Integration and Integration Testing

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.6.1 Integrate software units | [ ] | Build process |
| 5.6.2 Verify software integration | [ ] | Integration tests |
| 5.6.3 Integration test content | [ ] | test_api.py |
| 5.6.4 Integration regression tests | [ ] | CI pipeline |
| 5.6.5 Verify integration test results | [ ] | Test reports |
| 5.6.6 Integration test documentation | [ ] | Test protocols |
| 5.6.7 Evaluate integration for risk | [ ] | FMEA updates |

### 5.7 Software System Testing

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.7.1 Establish tests for requirements | [ ] | Validation protocol |
| 5.7.2 Test for proper event handling | [ ] | Error handling tests |
| 5.7.3 System test from device risk | [ ] | Risk-based testing |
| 5.7.4 Verify test results | [ ] | Validation report |
| 5.7.5 Test documentation | [ ] | Test records |

### 5.8 Software Release

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 5.8.1 Ensure software verified | [ ] | Validation complete |
| 5.8.2 Document known anomalies | [ ] | Release notes |
| 5.8.3 Evaluate known anomalies | [ ] | Risk assessment |
| 5.8.4 Document released versions | [ ] | Version tags |
| 5.8.5 Document how software delivered | [ ] | Deployment guide |
| 5.8.6 Archive software release | [ ] | Git tags, artifacts |
| 5.8.7 Ensure repeatability | [ ] | Build scripts |
| 5.8.8 Release identification | [ ] | Semantic versioning |

---

## Clause 6: Software Maintenance Process

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 6.1 Establish maintenance plan | [ ] | Maintenance SOP |
| 6.2 Problem and modification analysis | [ ] | Issue tracking |
| 6.3 Implement modifications | [ ] | Change control |

---

## Clause 7: Software Risk Management Process

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 7.1 Identify hazardous situations | [ ] | FMEA |
| 7.2 Document risk control measures | [ ] | FMEA mitigations |
| 7.3 Verify risk control measures | [ ] | Test traceability |
| 7.4 Risk management of changes | [ ] | Change impact analysis |

---

## Clause 8: Software Configuration Management Process

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 8.1 Configuration identification | [ ] | Version control |
| 8.2 Change control | [ ] | Git branches, PRs |
| 8.3 Configuration status accounting | [ ] | Release notes |

---

## Clause 9: Software Problem Resolution Process

| Requirement | Status | Evidence |
|-------------|--------|----------|
| 9.1 Prepare problem reports | [ ] | Issue templates |
| 9.2 Investigate problems | [ ] | Investigation logs |
| 9.3 Advise relevant parties | [ ] | Notification procedure |
| 9.4 Use change control | [ ] | PR process |
| 9.5 Maintain records | [ ] | Issue history |
| 9.6 Analyze problems for trends | [ ] | Periodic review |
| 9.7 Verify problem resolution | [ ] | Regression tests |
| 9.8 Document test content | [ ] | Test records |

---

## SOUP (Software of Unknown Provenance) Management

### Identified SOUP Components

| SOUP | Version | Risk | Verification |
|------|---------|------|--------------|
| Python | 3.10+ | Low | Version check |
| FastAPI | 0.100+ | Low | API tests |
| FAISS | Latest | Medium | Index tests |
| sentence-transformers | Latest | Medium | Embedding tests |
| mlx-lm | Latest | High | Model tests |
| pandas | 2.0+ | Low | Data tests |

### SOUP Risk Assessment

| SOUP | Failure Mode | Risk Control |
|------|--------------|--------------|
| FAISS | Index corruption | Integrity checks, rebuild capability |
| sentence-transformers | Wrong embeddings | Version pinning, validation tests |
| mlx-lm | Model errors | Output validation, fallback options |

---

## Document History

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2026-01-11 | Initial checklist | Validation Team |

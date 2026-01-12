# ACMG Variant Classification System - Limitations

**FOR RESEARCH USE ONLY**

This document describes the known limitations of the ACMG Variant Classification System. Users should review these limitations before using the system for any research purposes.

---

## 1. Classification Accuracy

### Current Performance

| Metric | Current | Target |
|--------|---------|--------|
| Pathogenic Sensitivity | 0% | ≥90% |
| Benign Specificity | 0% | ≥90% |
| Likely Pathogenic Accuracy | 0% | ≥85% |
| Likely Benign Accuracy | 0% | ≥85% |

**Note**: The classification model requires fine-tuning on ACMG-labeled training data to meet performance targets. Currently, the system defaults to "Uncertain_significance" for most variants when using the fallback model.

### VUS Rate

Without a fine-tuned model, the system reports a high VUS (Variant of Uncertain Significance) rate. This is a conservative approach that avoids incorrect pathogenicity calls but provides limited clinical utility.

### Catastrophic Misclassification

The validation suite confirms 0% catastrophic misclassification (Pathogenic→Benign or Benign→Pathogenic), meeting the critical safety requirement.

---

## 2. Database Coverage

### Gene Panel

- **Panel**: NGSgenes (314 genes)
- **Coverage**: Cardiac + hereditary cancer genes
- **Variants Indexed**: ~269,000 variants

### Variants NOT Covered

- Variants outside the NGSgenes panel
- Variants without GRCh37 coordinates (liftover failures)
- Structural variants (CNVs, large deletions)
- Variants not present in dbNSFP 5.3.1a

### Coordinate System

- **Current**: GRCh37 (hg19) only
- GRCh38 queries will return incorrect results
- Users must ensure coordinates match the genome build

---

## 3. Input Limitations

### Supported Variant Types

| Type | Supported | Notes |
|------|-----------|-------|
| SNVs | Yes | Full support |
| Small indels | Partial | Limited annotation coverage |
| MNVs | No | Not supported |
| Multi-allelic | No | Must be decomposed first |
| Structural variants | No | Not supported |

### Input Validation

The system validates:
- Chromosome: 1-22, X, Y, M/MT
- Position: Positive integers only
- Alleles: A, C, G, T only (IUPAC codes not supported)

---

## 4. Annotation Data Limitations

### dbNSFP Version

- **Version**: 5.3.1a
- **Build Date**: See database metadata
- Annotations reflect the state of computational predictors at time of release

### Predictor Scores

Not all variants have scores for all predictors. Missing scores indicate:
- The predictor doesn't cover that variant type
- The variant is in a region not assessed by that predictor
- Data was not available at dbNSFP compilation time

### ClinVar Data

- ClinVar annotations may be outdated
- Star ratings and review status are not currently displayed
- Conflicting interpretations are not reconciled

---

## 5. Model Limitations

### Classification Algorithm

- Uses LLM-based classification (requires fine-tuning)
- Falls back to Ollama llama3.2:3b when fine-tuned model unavailable
- Fallback model not trained for ACMG classification

### Confidence Scores

- Confidence values are estimates, not calibrated probabilities
- High confidence does not guarantee correctness
- Always requires expert review

### ACMG Criteria Extraction

- Criteria codes are parsed from model output
- May include hallucinated or invalid criteria codes
- Extracted criteria require manual verification

---

## 6. Performance Limitations

### Response Times

| Operation | P95 Latency | Notes |
|-----------|-------------|-------|
| Health check | <100ms | |
| Variant lookup | <500ms | Database only |
| Classification | <2000ms | With model inference |

### Concurrent Requests

- No explicit rate limiting configured
- Heavy load may increase latency
- Model inference is sequential (not parallelized)

---

## 7. Security Considerations

### Data Privacy

- Variants are stored locally in FAISS database
- No data is sent to external services (except Ollama if configured)
- API has no authentication by default

### Input Sanitization

- Basic input validation implemented
- SQL injection N/A (no SQL database)
- XSS protection via JSON-only API responses

---

## 8. Regulatory Status

### Not for Clinical Use

This system is:
- **NOT** FDA cleared or approved
- **NOT** CE marked
- **NOT** validated for diagnostic use
- **FOR RESEARCH USE ONLY**

### Required Reviews

All classifications from this system require:
1. Independent expert review
2. Verification against primary literature
3. Confirmation with orthogonal methods
4. Clinical correlation by qualified geneticist

---

## 9. Known Issues

### Test Suite Status

From validation testing (2024):
- 117 tests passing
- 4 tests failing (classification accuracy thresholds)
- All input validation tests pass
- All API endpoint tests pass

### Future Enhancements Needed

1. Fine-tune classification model on ACMG training data
2. Add GRCh38 coordinate support
3. Implement reference allele validation
4. Add API input validation layer
5. Display ClinVar review status

---

## 10. References

- [ACMG/AMP Guidelines (Richards et al., 2015)](https://doi.org/10.1038/gim.2015.30)
- [dbNSFP Documentation](https://sites.google.com/site/jpaborern/dbnsfp)
- [ClinGen Sequence Variant Interpretation](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/)

---

## Contact

For questions about these limitations, please open an issue in the project repository.

**Last Updated**: January 2026

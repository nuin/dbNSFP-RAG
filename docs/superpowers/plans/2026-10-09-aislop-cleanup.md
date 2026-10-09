# Aislop cleanup implementation plan

**Goal:** Resolve the reported 311 findings while preserving variant scoring and command behavior.

**Architecture:** Refactor by responsibility, keep CLI entry points and API routes compatible, and declare optional runtime dependencies. Preserve existing uncommitted changes and bundled script independence. Do not weaken scanner rules.

**Tech stack:** Python, FastAPI, pandas, SQLite, FAISS, pytest, Ruff, aislop.

**Spec:** User-provided findings and repository AGENTS.md.

## Work groups

- [x] API and scoring: split request schemas, model backends, route helpers, and criterion evaluators; fix error handling and unused values. Verify API and scoring tests.
- [x] Panel/inspection scripts and bundle: extract reusable parsing, reporting, and processing helpers; preserve standalone bundle operation and existing classification differences. Verify representative CLI and synthetic fixtures.
- [x] Upload/backfill scripts: narrow parsing failures, close files, extract coherent workflow stages, remove dead values. Verify synthetic data transformations.
- [x] Core ingestion/storage/training: split CLI build workflows, retain checkpoints and lazy optional dependencies, replace unsafe metadata deserialization, declare dependency extras, and fix logging and remaining lint findings. Verify storage round trips, malformed metadata, interrupted builds, and CLI smoke tests.
- [x] Integration: run repository tests and fresh full aislop scan; inspect changes and resolve regressions until no actionable findings remain.

## Review focus

- Existing metadata files must remain readable without executing pickle globals.
- Interrupted builds must save resumable state and return a non-success exit status.
- Missing or malformed external annotations must produce observable errors without fabricating evidence.
- Bundled scripts must work with their deployed directory structure.
- Medical classification thresholds and rules must remain unchanged by refactoring.

## Verification outcome

All 311 original findings were addressed without changing scanner rules. The final
scan reports 98/100 with one new false positive: `pickletools` is a Python standard
library module (`pickletools in sys.stdlib_module_names` is true), but aislop 0.18.1
reports it as an undeclared import. It remains imported normally because its opcode
parser prevents executable legacy pickle payloads, including cached extension
opcodes, from reaching the restricted reader.

Verification covers the complete pytest suite, original Ruff rule set on source
files, formatting, whitespace, API/OpenAPI and scoring comparisons against the
initial working-tree code, standalone bundle output comparisons, SQLite build
integration for ZIP/directory inputs, interrupted and limited build resumes,
cross-chromosome liftover, and malicious legacy metadata payloads. Full model
training and production database rebuilds were not run.

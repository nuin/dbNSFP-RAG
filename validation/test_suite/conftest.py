"""Pytest fixtures for ACMG Variant Classification System validation."""

import json
import os
import sys
from pathlib import Path
from typing import Generator

import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Return project root directory."""
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def gold_standard_path(project_root: Path) -> Path:
    """Return path to gold standard benchmark file."""
    return project_root / "validation" / "test_suite" / "gold_standard" / "clinvar_benchmark.json"


@pytest.fixture(scope="session")
def gold_standard(gold_standard_path: Path) -> list[dict]:
    """Load gold standard benchmark variants."""
    if not gold_standard_path.exists():
        pytest.skip(f"Gold standard file not found: {gold_standard_path}")
    with open(gold_standard_path) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def pathogenic_variants(gold_standard: list[dict]) -> list[dict]:
    """Filter gold standard to pathogenic variants only."""
    return [v for v in gold_standard if v["expected"] == "Pathogenic"]


@pytest.fixture(scope="session")
def likely_pathogenic_variants(gold_standard: list[dict]) -> list[dict]:
    """Filter gold standard to likely pathogenic variants only."""
    return [v for v in gold_standard if v["expected"] == "Likely_pathogenic"]


@pytest.fixture(scope="session")
def benign_variants(gold_standard: list[dict]) -> list[dict]:
    """Filter gold standard to benign variants only."""
    return [v for v in gold_standard if v["expected"] == "Benign"]


@pytest.fixture(scope="session")
def likely_benign_variants(gold_standard: list[dict]) -> list[dict]:
    """Filter gold standard to likely benign variants only."""
    return [v for v in gold_standard if v["expected"] == "Likely_benign"]


@pytest.fixture(scope="session")
def vus_variants(gold_standard: list[dict]) -> list[dict]:
    """Filter gold standard to VUS variants only."""
    return [v for v in gold_standard if v["expected"] == "Uncertain_significance"]


@pytest.fixture(scope="session")
def database_path(project_root: Path) -> Path:
    """Return path to variant database."""
    db_path = project_root / "data" / "vectordb" / "grch37-ngsgenes"
    if not db_path.exists():
        pytest.skip(f"Database not found: {db_path}")
    return db_path


@pytest.fixture(scope="session")
def vectorstore(database_path: Path):
    """Load vector store for testing."""
    from src.vectorstore import VariantVectorStore
    store = VariantVectorStore(str(database_path))
    store.load()
    return store


@pytest.fixture(scope="session")
def model_path(project_root: Path) -> Path:
    """Return path to classification model."""
    model_path = project_root / "models" / "acmg-classifier" / "model"
    return model_path


@pytest.fixture(scope="module")
def api_client():
    """Create test client for API."""
    from fastapi.testclient import TestClient
    from api.server import app
    return TestClient(app)


# Valid variant test data
@pytest.fixture
def valid_snv() -> dict:
    """Return a valid SNV for testing."""
    return {
        "chr": "17",
        "pos": 41197801,
        "ref": "T",
        "alt": "A"
    }


@pytest.fixture
def valid_brca1_variant() -> dict:
    """Return a well-known BRCA1 variant."""
    return {
        "chr": "17",
        "pos": 41244936,
        "ref": "G",
        "alt": "A",
        "gene": "BRCA1",
        "expected": "Pathogenic"
    }


# Invalid variant test data
@pytest.fixture
def invalid_chromosome_variants() -> list[dict]:
    """Return variants with invalid chromosomes."""
    return [
        {"chr": "99", "pos": 100, "ref": "A", "alt": "G"},
        {"chr": "0", "pos": 100, "ref": "A", "alt": "G"},
        {"chr": "chr99", "pos": 100, "ref": "A", "alt": "G"},
        {"chr": "", "pos": 100, "ref": "A", "alt": "G"},
        {"chr": None, "pos": 100, "ref": "A", "alt": "G"},
        {"chr": "ABC", "pos": 100, "ref": "A", "alt": "G"},
    ]


@pytest.fixture
def invalid_position_variants() -> list[dict]:
    """Return variants with invalid positions."""
    return [
        {"chr": "17", "pos": -100, "ref": "A", "alt": "G"},
        {"chr": "17", "pos": 0, "ref": "A", "alt": "G"},
        {"chr": "17", "pos": None, "ref": "A", "alt": "G"},
        {"chr": "17", "pos": "abc", "ref": "A", "alt": "G"},
    ]


@pytest.fixture
def invalid_allele_variants() -> list[dict]:
    """Return variants with invalid alleles."""
    return [
        {"chr": "17", "pos": 100, "ref": "X", "alt": "G"},
        {"chr": "17", "pos": 100, "ref": "A", "alt": "X"},
        {"chr": "17", "pos": 100, "ref": "", "alt": "G"},
        {"chr": "17", "pos": 100, "ref": "A", "alt": ""},
        {"chr": "17", "pos": 100, "ref": None, "alt": "G"},
        {"chr": "17", "pos": 100, "ref": "A", "alt": None},
        {"chr": "17", "pos": 100, "ref": "A1", "alt": "G"},
        {"chr": "17", "pos": 100, "ref": "A", "alt": "G2"},
    ]


@pytest.fixture
def multi_allelic_variants() -> list[dict]:
    """Return multi-allelic variants."""
    return [
        {"chr": "17", "pos": 100, "ref": "A", "alt": "G,T"},
        {"chr": "17", "pos": 100, "ref": "A", "alt": "G,T,C"},
    ]


# Performance testing fixtures
@pytest.fixture
def performance_threshold_ms() -> int:
    """Maximum acceptable response time in milliseconds."""
    return 2000


@pytest.fixture
def classification_accuracy_threshold() -> float:
    """Minimum acceptable classification accuracy."""
    return 0.90


# Environment markers
def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "requires_model: marks tests that require the ML model"
    )
    config.addinivalue_line(
        "markers", "requires_database: marks tests that require the vector database"
    )
    config.addinivalue_line(
        "markers", "integration: marks integration tests"
    )

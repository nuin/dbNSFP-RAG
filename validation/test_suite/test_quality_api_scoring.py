"""Regression coverage for public scoring and API boundaries."""

import gc
import weakref
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from api import server
from src.acmg_scoring import (
    combine_criteria,
    create_criterion,
    evaluate_all_criteria,
    evaluate_population_criteria,
)
from src.clinvar_api import ClinVarClient


@pytest.mark.parametrize(
    ("frequency", "expected"),
    [
        (None, ["PM2"]),
        (0.000099, ["PM2"]),
        (0.0001, []),
        (0.01, []),
        (0.010001, ["BS1"]),
        (0.05, ["BA1"]),
    ],
)
def test_population_threshold_boundaries(frequency, expected):
    assert [c.code for c in evaluate_population_criteria(frequency)] == expected


@pytest.mark.parametrize(
    ("codes", "classification", "confidence", "rule"),
    [
        (["BA1", "PVS1"], "Benign", 0.95, "BA1: Allele frequency ≥5%"),
        (["BS1", "BP4", "PVS1"], "Likely_benign", 0.8, "1 BS + 1 BP"),
        (["PVS1", "PS1"], "Pathogenic", 0.95, "PVS1 + ≥1 PS"),
        (["PVS1", "PM2"], "Likely_pathogenic", 0.8, "PVS1 + 1 PM"),
        (["PS1", "PM1", "PM2", "PM4"], "Pathogenic", 0.85, "1 PS + ≥3 PM"),
        (["PS1", "PP2", "PP3"], "Likely_pathogenic", 0.7, "1 PS + ≥2 PP"),
        (["PM1", "PM2", "PP2", "PP3"], "Likely_pathogenic", 0.7, "2 PM + ≥2 PP"),
        (["PM2", "BS1"], "Uncertain_significance", 0.5, ""),
    ],
)
def test_combination_priority_and_confidence(codes, classification, confidence, rule):
    result = combine_criteria([create_criterion(code, met=True) for code in codes])
    assert (result.classification.value, result.confidence, result.rule_applied) == (
        classification,
        confidence,
        rule,
    )
    assert [c.code for c in result.criteria_met] == codes


def test_all_criteria_preserves_display_order_and_splice_boundary():
    criteria = evaluate_all_criteria({"consequence": "synonymous_variant", "spliceai_ag": 0.2})
    assert [c.code for c in criteria] == [
        "PVS1",
        "PS2",
        "PS3",
        "PS4",
        "PM1",
        "PM2",
        "PM3",
        "PM4",
        "PS1",
        "PM5",
        "PM6",
        "PP1",
        "PP2",
        "PP3",
        "PP4",
        "PP5",
        "BA1",
        "BS1",
        "BS2",
        "BS3",
        "BS4",
        "BP1",
        "BP2",
        "BP3",
        "BP4",
        "BP5",
        "BP6",
        "BP7",
    ]
    bp7 = criteria[-1]
    assert bp7.status.value == "not_met"
    assert "SpliceAI max=0.20 ≥ 0.2" in bp7.evidence


def test_clinvar_cache_is_owned_by_client(monkeypatch):
    client = ClinVarClient("research@example.org")
    requests = []

    def make_request(endpoint, params):
        requests.append(endpoint)
        return {"esearchresult": {"idlist": []}}

    monkeypatch.setattr(client, "_make_request", make_request)
    assert client.get_pathogenic_in_gene("BRCA1") == []
    assert client.get_pathogenic_in_gene("BRCA1") == []
    assert requests == ["esearch.fcgi"]
    reference = weakref.ref(client)
    monkeypatch.undo()
    del client
    gc.collect()
    assert reference() is None


@pytest.fixture
def record_client(monkeypatch):
    metadata = {
        "chr": "1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": "TEST",
        "consequence": "synonymous_variant",
        "gnomad_af": 0.1,
    }
    record = {"id": "1_100_A_G", "metadata": metadata, "document": "Stored annotation"}
    store = SimpleNamespace(
        get_by_id=lambda variant_id: record if variant_id == record["id"] else None,
        count=lambda: 1,
        list_chromosomes=lambda: [{"chromosome": "1", "count": 1}],
    )
    monkeypatch.setattr(server, "get_database", lambda build: store)
    monkeypatch.setattr(server, "CLINVAR_CLIENT", None)
    with TestClient(server.app) as client:
        yield client


def test_classification_coordinates_and_hgvs_agree(record_client):
    coordinates = record_client.get("/classify?chr=1&pos=100&ref=A&alt=G")
    hgvs = record_client.get("/classify?hgvs=1:g.100A>G")
    assert coordinates.status_code == hgvs.status_code == 200
    assert coordinates.json() == hgvs.json()
    assert coordinates.json()["acmg_classification"] == "Benign"
    assert len(coordinates.json()["all_criteria"]) == 28


def test_lookup_and_browse_use_api_database_provider(record_client):
    lookup = record_client.get("/lookup?hgvs=1:g.100A>G")
    assert lookup.status_code == 200
    assert lookup.json()["full_annotation"] == "Stored annotation"
    assert lookup.json()["variant_id"] == "1_100_A_G"
    assert record_client.get("/chromosomes").json() == {
        "genome_build": "GRCh37",
        "chromosomes": [{"chromosome": "1", "count": 1}],
    }


def test_missing_variant_and_bad_input_are_distinct(record_client):
    missing = record_client.get("/classify?chr=1&pos=999&ref=A&alt=G")
    assert missing.status_code == 200
    assert missing.json()["acmg_classification"] == "Not_in_database"
    assert missing.json()["confidence"] == 0.0
    assert record_client.get("/lookup?chr=1&pos=999&ref=A&alt=G").status_code == 404
    assert record_client.get("/classify?hgvs=invalid").status_code == 400

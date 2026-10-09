"""Malformed upstream records must leave evidence unevaluated, not crash scoring."""

from types import SimpleNamespace

import pytest

from src import repeat_api
from src.acmg_scoring import evaluate_all_criteria
from src.clinvar_api import ClinVarClient


@pytest.mark.parametrize(
    "payload",
    [
        [],
        {"esearchresult": []},
        {"esearchresult": {"idlist": "123"}},
        {"esearchresult": {"idlist": [123]}},
    ],
)
def test_malformed_clinvar_search_is_not_evaluated(monkeypatch, payload):
    client = ClinVarClient("research@example.org")
    monkeypatch.setattr(client, "_make_request", lambda *args: payload)
    criteria = evaluate_all_criteria(
        {"gene": "TEST", "aa_pos": 12, "aa_ref": "Arg", "aa_alt": "Gly"}, client
    )
    affected = [c for c in criteria if c.code in ("PS1", "PM5")]
    assert [c.status.value for c in affected] == ["not_evaluated", "not_evaluated"]
    assert all("ClinVar lookup failed" in c.evidence for c in affected)


@pytest.mark.parametrize(
    "payload",
    [
        {"result": []},
        {"result": {"uids": ["123"], "123": None}},
        {"result": {"uids": ["123"], "123": {"clinical_significance": ["Pathogenic"]}}},
        {"result": {"uids": ["123"], "123": {"genes": [None]}}},
        {"result": {"uids": ["123"], "123": {"protein_change": ["p.Arg12Gly"]}}},
    ],
)
def test_malformed_clinvar_summary_is_not_evaluated(monkeypatch, payload):
    client = ClinVarClient("research@example.org")

    def response(endpoint, params):
        if endpoint == "esearch.fcgi":
            return {"esearchresult": {"idlist": ["123"]}}
        return payload

    monkeypatch.setattr(client, "_make_request", response)
    criteria = evaluate_all_criteria(
        {"gene": "TEST", "aa_pos": 12, "aa_ref": "Arg", "aa_alt": "Gly"}, client
    )
    affected = [c for c in criteria if c.code in ("PS1", "PM5")]
    assert [c.status.value for c in affected] == ["not_evaluated", "not_evaluated"]
    assert all("ClinVar lookup failed" in c.evidence for c in affected)


@pytest.mark.parametrize("payload", [[], {"rmsk": [None]}, {"rmsk": "invalid"}])
def test_malformed_repeat_response_is_not_evaluated(monkeypatch, payload):
    response = SimpleNamespace(raise_for_status=lambda: None, json=lambda: payload)
    monkeypatch.setattr(repeat_api.requests, "get", lambda *args, **kwargs: response)
    repeat_api.is_in_repeat_region.cache_clear()
    try:
        criteria = evaluate_all_criteria(
            {"chr": "1", "pos": 100, "consequence": "inframe_deletion"}
        )
        bp3 = next(c for c in criteria if c.code == "BP3")
        assert bp3.status.value == "not_evaluated"
        assert "Repeat region lookup failed" in bp3.evidence
    finally:
        repeat_api.is_in_repeat_region.cache_clear()


def test_clinvar_valid_summaries_retain_ps1_and_pm5_evidence(monkeypatch):
    client = ClinVarClient("research@example.org")

    def response(endpoint, params):
        if endpoint == "esearch.fcgi":
            return {"esearchresult": {"idlist": ["123", "456"]}}
        return {
            "result": {
                "uids": ["123", "456"],
                "123": {
                    "protein_change": "p.Arg12Gly",
                    "clinical_significance": {"description": "Pathogenic"},
                    "genes": [{"symbol": "TEST"}],
                },
                "456": {
                    "title": "TEST:c.1A>G (p.Arg12His)",
                    "clinical_significance": {"description": "Pathogenic"},
                },
            }
        }

    monkeypatch.setattr(client, "_make_request", response)
    criteria = evaluate_all_criteria(
        {"gene": "TEST", "aa_pos": 12, "aa_ref": "Arg", "aa_alt": "Gly"}, client
    )
    affected = [c for c in criteria if c.code in ("PS1", "PM5")]
    assert [c.status.value for c in affected] == ["met", "met"]
    assert all("p.Arg12Gly" in c.evidence and "p.Arg12His" in c.evidence for c in affected)


@pytest.mark.parametrize(
    "repeats, expected", [([], "not_met"), ([{"repClass": "LINE", "repFamily": "L1"}], "met")]
)
def test_valid_repeat_annotations_keep_status(monkeypatch, repeats, expected):
    response = SimpleNamespace(raise_for_status=lambda: None, json=lambda: {"rmsk": repeats})
    monkeypatch.setattr(repeat_api.requests, "get", lambda *args, **kwargs: response)
    repeat_api.is_in_repeat_region.cache_clear()
    try:
        criteria = evaluate_all_criteria(
            {"chr": "1", "pos": 100, "consequence": "inframe_deletion"}
        )
        assert next(c for c in criteria if c.code == "BP3").status.value == expected
    finally:
        repeat_api.is_in_repeat_region.cache_clear()


def test_unexpected_client_programming_errors_propagate():
    def broken_client(gene):
        raise AttributeError("client implementation bug")

    client = SimpleNamespace(get_pathogenic_in_gene=broken_client)
    with pytest.raises(AttributeError, match="implementation bug"):
        evaluate_all_criteria(
            {"gene": "TEST", "aa_pos": 12, "aa_ref": "Arg", "aa_alt": "Gly"}, client
        )

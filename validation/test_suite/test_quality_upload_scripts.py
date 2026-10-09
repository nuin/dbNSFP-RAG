"""Regression coverage for the upload and backfill processing stages."""

import json
from io import StringIO
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from scripts import (
    backfill_spliceai_dual,
    backfill_syn_spliceai,
    build_upload_v3_1,
    build_upload_v3_2,
    build_upload_v3_proper,
    finalize_spliceai_labels,
    reannotate_gnomad_common,
    resolve_needs_hgvs_myvariant,
    validate_spliceai_api,
)


@pytest.mark.parametrize(
    ("value", "expected"),
    [(None, None), (".", None), ("0 (assumed)", 0), (".;bad;0.29;0.8", 0.29), ("bad", None)],
)
def test_multitranscript_score_parsing(value, expected):
    assert build_upload_v3_2._f(value) == expected


@pytest.mark.parametrize(
    ("faf", "revel", "spliceai", "passes", "reason"),
    [
        ("0.01", "0.2", "0.1", True, "Likely_benign"),
        ("0.001", "0.2", "0.1", False, "faf_0.00100"),
        ("0.01", "0.29", "0.1", False, "revel_0.290"),
        ("0.01", "0.2", "0.11", False, "spliceai_0.110"),
        ("", "", "", False, "faf_none|revel_none"),
    ],
)
def test_missense_reclassification_thresholds(faf, revel, spliceai, passes, reason):
    row = {"FAF95_grpmax": faf, "REVEL": revel, "SpliceAI_masked": spliceai}
    classification, matched = build_upload_v3_1.classify_missense(row, "A", "G")
    assert matched is passes
    assert reason in classification


def test_clinvar_removal_precedes_source_deduplication():
    variant = {
        "gene": "GOT2",
        "transcript": "NM_002080.5",
        "hgvs_c": "c.816C>T",
        "chr_grch37": "16",
        "pos_grch37": "58742000",
        "ref": "C",
        "alt": "T",
    }
    final = pd.DataFrame(
        [
            {**variant, "source": "clinvar", "classification": "Benign ClinVar"},
            {**variant, "source": "synonymous_catalog", "classification": "Benign FAF >5%"},
        ]
    )
    gnomad = pd.DataFrame([{**variant, "classification": "Benign gnomAD"}])
    result = build_upload_v3_proper.combine_non_clinvar_sources(final, gnomad)
    assert result["source"].tolist() == ["synonymous_catalog"]
    assert result["classification"].tolist() == ["Benign FAF >5%"]
    assert not any(column.startswith("_") for column in result.columns)


def test_spliceai_backfill_preserves_measured_values():
    row = {
        "chr_grch37": "1",
        "pos_grch37": "100",
        "ref": "A",
        "alt": "G",
        "SpliceAI_masked": "0 (assumed)",
        "SpliceAI_unmasked": "",
    }
    frame = pd.DataFrame([row, {**row, "SpliceAI_masked": "0.25", "SpliceAI_unmasked": "0.4"}])
    backfill_spliceai_dual.fill_spliceai_scores(
        frame, {("1", 100): 200}, {("1", 200, "A", "G"): 0.1}, {("1", 200, "A", "G"): 0.2}
    )
    assert frame["SpliceAI_masked"].tolist() == ["0.1000", "0.25"]
    assert frame["SpliceAI_unmasked"].tolist() == ["0.2000", "0.4"]


@pytest.mark.parametrize(
    "reader", [backfill_syn_spliceai.parse_vcf, finalize_spliceai_labels.read_spliceai_scores]
)
def test_spliceai_reader_keeps_valid_entries_after_malformed_score(tmp_path, caplog, reader):
    path = tmp_path / "scores.vcf"
    path.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE|bad|0|0|0,G|GENE|0.2|0.3|0.1|0\n"
    )
    assert reader(path) == {("1", 100, "A", "G"): 0.3}
    assert "Invalid SpliceAI entry" in caplog.text


def test_reannotation_uses_synonymous_catalog_when_coding_lookup_misses():
    row = {
        "gene": "GENE",
        "chr_hg38": "1",
        "pos_hg38": "200",
        "ref": "A",
        "alt": "G",
        "FAF95_grpmax": "0.2",
        "gnomad_v41_af_joint": "0.3",
        "gnomad_v41_grpmax_anc": "nfe",
        "gnomad_v41_nhomalt": "2",
        "gnomad_v41_filter": "PASS",
    }
    annotation = reannotate_gnomad_common.lookup_annotation(
        row["gene"],
        ("1", 200, "A", "G"),
        lambda gene: {},
        {("1", 200): "100"},
        {("1", 100, "A", "G"): {"transcript": "NM_1.2", "hgvs_c": "c.3A>G"}},
    )
    assert annotation == ("NM_1.2", "c.3A>G", "1", "100", "synonymous", "syn_catalog")


def test_myvariant_resolution_prefers_exact_transcript_version():
    result = {
        "snpeff": {
            "ann": [
                {"feature_id": "NM_1.1", "hgvs_c": "c.1A>G", "effect": "old"},
                {"feature_id": "NM_1.2", "hgvs_c": "c.2A>G", "effect": "current"},
                {"feature_id": "NM_2.2", "hgvs_c": "c.3A>G", "effect": "other"},
            ]
        }
    }
    assert resolve_needs_hgvs_myvariant.extract_hgvsc(result, "NM_1.2") == (
        "c.2A>G",
        "current",
        "NM_1.2",
    )
    assert resolve_needs_hgvs_myvariant.extract_hgvsc(result, "NM_3.1") is None


def test_myvariant_cache_reports_corruption(tmp_path, monkeypatch, capsys):
    path = Path(tmp_path) / "cache.json"
    path.write_text("{broken")
    monkeypatch.setattr(resolve_needs_hgvs_myvariant, "CACHE", path)
    assert resolve_needs_hgvs_myvariant.load_cache() == {}
    assert "Unable to read MyVariant cache" in capsys.readouterr().err


@pytest.mark.parametrize("payload", [None, {"error": "upstream failure"}, [None], [{"query": []}]])
def test_myvariant_invalid_batch_retries_and_returns_failure(payload, monkeypatch, capsys):
    response = Mock()
    response.json.return_value = payload
    post = Mock(return_value=response)
    monkeypatch.setattr(resolve_needs_hgvs_myvariant.requests, "post", post)
    monkeypatch.setattr(resolve_needs_hgvs_myvariant.time, "sleep", lambda delay: None)
    assert resolve_needs_hgvs_myvariant.query_batch(["chr1:g.1A>G"], {}) == {
        "chr1:g.1A>G": {"notfound": True, "error": "batch_failed"}
    }
    assert post.call_count == 3
    assert capsys.readouterr().err.count("batch error attempt") == 3


def test_myvariant_recovers_after_invalid_batch(monkeypatch):
    response = Mock()
    hit = {"query": "chr1:g.1A>G", "snpeff": {"ann": []}}
    response.json.side_effect = [{"error": "upstream failure"}, [hit]]
    post = Mock(return_value=response)
    monkeypatch.setattr(resolve_needs_hgvs_myvariant.requests, "post", post)
    monkeypatch.setattr(resolve_needs_hgvs_myvariant.time, "sleep", lambda delay: None)
    assert resolve_needs_hgvs_myvariant.query_batch([hit["query"]], {}) == {hit["query"]: hit}
    assert post.call_count == 2


@pytest.mark.parametrize(
    "payload", [None, [], {"scores": {}}, {"scores": [None]}, {"scores": [{"DS_AG": []}]}]
)
def test_spliceai_invalid_response_retries_and_returns_none(payload, monkeypatch, capsys):
    request = Mock(side_effect=lambda *args, **kwargs: StringIO(json.dumps(payload)))
    monkeypatch.setattr(validate_spliceai_api.urllib.request, "urlopen", request)
    monkeypatch.setattr(validate_spliceai_api.time, "sleep", lambda delay: None)
    assert validate_spliceai_api.api_ds_max("1-1-A-G", 1) is None
    assert request.call_count == 3
    assert capsys.readouterr().err.count("SpliceAI API attempt") == 3


def test_spliceai_recovers_after_invalid_response(monkeypatch):
    responses = [StringIO("null"), StringIO('{"scores": [{"DS_AG": "0.1", "DS_AL": 0.3}]}')]
    request = Mock(side_effect=responses)
    monkeypatch.setattr(validate_spliceai_api.urllib.request, "urlopen", request)
    monkeypatch.setattr(validate_spliceai_api.time, "sleep", lambda delay: None)
    assert validate_spliceai_api.api_ds_max("1-1-A-G", 1) == 0.3
    assert request.call_count == 2


def test_myvariant_does_not_retry_unrelated_type_errors(monkeypatch):
    post = Mock(side_effect=TypeError("invalid request configuration"))
    monkeypatch.setattr(resolve_needs_hgvs_myvariant.requests, "post", post)
    with pytest.raises(TypeError, match="invalid request configuration"):
        resolve_needs_hgvs_myvariant.query_batch(["chr1:g.1A>G"], {})
    assert post.call_count == 1


def test_spliceai_does_not_retry_unrelated_type_errors(monkeypatch):
    request = Mock(side_effect=TypeError("invalid request configuration"))
    monkeypatch.setattr(validate_spliceai_api.urllib.request, "urlopen", request)
    with pytest.raises(TypeError, match="invalid request configuration"):
        validate_spliceai_api.api_ds_max("1-1-A-G", 1)
    assert request.call_count == 1

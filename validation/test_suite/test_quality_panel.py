"""Regression coverage for the standalone panel pipeline and its deployment bundle."""

from __future__ import annotations

import gzip
import importlib
import json
import sys
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(params=["scripts", "cp_new_bundle/scripts"])
def panel_module(request, monkeypatch):
    directory = ROOT / request.param
    monkeypatch.syspath_prepend(str(directory))
    names = (
        "cp_new",
        "cp_new_data",
        "cp_new_rules",
        "classify_cp_new",
        "inspect_variant",
        "inspector_sources",
        "inspector_report",
        "pull_clinvar_benign",
        "resolve_hgvs_vv",
    )
    for name in names:
        monkeypatch.delitem(sys.modules, name, raising=False)

    def load(name):
        return importlib.import_module(name)

    yield load, request.param.startswith("cp_new_bundle")
    for name in names:
        sys.modules.pop(name, None)


def test_classification_keeps_deployment_specific_rules(panel_module):
    load, bundled = panel_module
    rules = load("cp_new_rules")
    missense = {
        "aaref": "A",
        "aaalt": "G",
        "codon_degeneracy": "2",
        "spliceai_ds_max_masked": "0.01",
        "phastCons100way_vertebrate": "0.5",
    }
    assert rules.is_synonymous(missense) is bundled
    result = rules.classify(missense)
    assert (result is not None) is bundled
    assert rules.safe_float(".;0.2;.") == (None if bundled else 0.2)
    synonymous = {**missense, "aaref": "A", "aaalt": "A", "phyloP100way_vertebrate": "0.01"}
    assert rules.classify(synonymous, "c.100+6A>G") is not None
    assert rules.classify(synonymous, "c.100+7A>G") is None


def test_export_uses_resolution_cache_and_deduplicates(panel_module, tmp_path):
    load, _ = panel_module
    runner = load("classify_cp_new")
    row = {
        "#chr": "1",
        "genename": "OTHER;GENE",
        "hg19_chr": "1",
        "hg19_pos(1-based)": "50",
        "ref": "A",
        "alt": "G",
        "gnomad_v41_faf95_grpmax": "0.06",
        "HGVSc_snpEff": "c.10A>G;c.11A>G",
    }
    tsv = tmp_path / "GENE.tsv"
    pd.DataFrame([row, row]).to_csv(tsv, sep="\t", index=False)
    key = "1-50-A-G|NM_1.1"
    run = runner.ClassificationRun(
        {("1", "GENE"): [(0, 100, "NM_1.1")]},
        {"GENE": "NM_1.1"},
        {key: "c.20A>G"},
        {key: "ok"},
        runner.initial_counters(),
    )
    rows = runner.classify_file(tsv, run)
    assert len(rows) == 1
    assert rows[0]["gene"] == "GENE"
    assert rows[0]["hgvs_c"] == "c.20A>G"
    assert rows[0]["classification"] == "Benign FAF >5%"
    assert rows[0]["vv_status"] == "ok"
    assert run.counters["total_rows"] == 2
    assert run.counters["classified"] == 1


def test_annotation_joins_vcf_coordinates_and_writes_integer_positions(
    panel_module, tmp_path, monkeypatch
):
    load, _ = panel_module
    panel = load("cp_new")
    gene = pd.DataFrame(
        [
            {
                "#chr": "1",
                "pos(1-based)": 100.0,
                "hg19_pos(1-based)": 90.0,
                "ref": "A",
                "alt": "G",
                "_gene_upper": "GENE",
            }
        ]
    )
    hit = pd.DataFrame(
        [
            {
                "gnomad_chr": "1",
                "gnomad_pos": 100,
                "gnomad_ref": "A",
                "gnomad_alt": "G",
                "gnomad_v41_faf95_grpmax": "0.06",
            }
        ]
    )
    monkeypatch.setattr(panel, "fetch_gnomad_region", lambda *args: hit)
    options = panel.AnnotationOptions(tmp_path, None, None, False, True)
    assert panel.annotate_and_write_gene("GENE", gene, options) == 1
    result = pd.read_csv(tmp_path / "GENE.tsv", sep="\t", dtype=str)
    assert result.loc[0, "pos(1-based)"] == "100"
    assert result.loc[0, "hg19_pos(1-based)"] == "90"
    assert result.loc[0, "gnomad_v41_faf95_grpmax"] == "0.06"
    assert "_gene_upper" not in result.columns


def test_inspector_recovers_invalid_json_cache_with_diagnostic(
    panel_module, tmp_path, monkeypatch, caplog
):
    load, _ = panel_module
    sources = load("inspector_sources")
    monkeypatch.setattr(sources, "_HTTP_CACHE_DIR", tmp_path)
    cache = tmp_path / "variant.json"
    cache.write_text("{corrupt")
    response = Mock()
    response.json.return_value = {"vcf": {"position": 50}}
    session = Mock()
    session.get.return_value = response
    monkeypatch.setattr(sources, "SESSION", session)
    result = sources.http_get_json("https://example.test/variant", "variant")
    assert result == {"vcf": {"position": 50}}
    assert json.loads(cache.read_text()) == result
    assert "Cannot read API cache" in caplog.text
    session.get.assert_called_once()


def test_inspector_no_external_flag_skips_rsid_resolution(panel_module, monkeypatch, capsys):
    load, _ = panel_module
    inspector = load("inspect_variant")
    request = Mock()
    monkeypatch.setattr(inspector, "fetch_myvariant", request)
    monkeypatch.setattr(inspector, "find_in_classifications", lambda variant: [])
    monkeypatch.setattr(
        sys, "argv", ["inspect_variant.py", "lookup", "--rsid", "rs123", "--no-external", "--json"]
    )
    assert inspector.main() == 0
    assert json.loads(capsys.readouterr().out)["variant"]["rsid"] == "rs123"
    request.assert_not_called()


def test_clinvar_reader_preserves_genomic_alleles_and_deduplicates(panel_module, tmp_path):
    load, _ = panel_module
    clinvar = load("pull_clinvar_benign")
    record = {
        "GeneSymbol": "GENE",
        "ClinicalSignificance": "Likely benign",
        "ReviewStatus": "criteria provided, single submitter",
        "Assembly": "GRCh37",
        "Start": "50",
        "Name": "NM_2.1(GENE):c.10C>T (p.Ala4Val)",
        "ReferenceAllele": "na",
        "AlternateAllele": "na",
        "Chromosome": "1",
        "Type": "SNV",
    }
    path = tmp_path / "clinvar.gz"
    with gzip.open(path, "wt") as handle:
        pd.DataFrame([record, record, {**record, "Assembly": "GRCh38"}]).to_csv(
            handle, sep="\t", index=False
        )
    rows = clinvar.read_benign_rows(path, {"GENE"}, {"GENE": "NM_1.1"})
    assert len(rows) == 1
    assert rows[0]["transcript"] == "NM_1.1"
    assert rows[0]["hgvs_c"] == "c.10C>T"
    assert rows[0]["ref"] == rows[0]["alt"] == ""


def test_variantvalidator_cache_and_transcript_mismatch(panel_module):
    load, _ = panel_module
    resolver = load("resolve_hgvs_vv")
    variant = resolver.ResolutionRequest("1", 50, "A", "G", "NM_1.1")
    response = Mock(status_code=200)
    response.json.return_value = {"NM_2.1:c.10A>G": {}}
    session = Mock()
    session.get.return_value = response
    cache = {}
    assert resolver.resolve_one(variant, cache, session) == (
        "c.10A>G",
        "flagged:tx_mismatch:NM_2.1",
        True,
    )
    assert resolver.resolve_one(variant, cache, session) == (
        "c.10A>G",
        "flagged:tx_mismatch:NM_2.1",
        False,
    )
    session.get.assert_called_once()

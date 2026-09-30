"""Tests for scripts/resolve_hgvs_vv.py (REST + local hgvs backends)."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "resolve_hgvs_vv.py"
_spec = importlib.util.spec_from_file_location("resolve_hgvs_vv", SCRIPT)
vv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vv)

AC17 = "NC_000017.10"


class TestVcfToHgvsG:
    def test_snv(self):
        assert vv.vcf_to_hgvs_g(AC17, 41197801, "T", "A") == f"{AC17}:g.41197801T>A"

    def test_padded_deletion(self):
        assert vv.vcf_to_hgvs_g(AC17, 100, "ACT", "A") == f"{AC17}:g.101_102del"

    def test_single_base_deletion(self):
        assert vv.vcf_to_hgvs_g(AC17, 100, "AC", "A") == f"{AC17}:g.101del"

    def test_padded_insertion(self):
        assert vv.vcf_to_hgvs_g(AC17, 100, "A", "AGG") == f"{AC17}:g.100_101insGG"

    def test_delins(self):
        assert vv.vcf_to_hgvs_g(AC17, 100, "AT", "GC") == f"{AC17}:g.100_101delinsGC"

    def test_mnv_with_shared_suffix_trimmed(self):
        assert vv.vcf_to_hgvs_g(AC17, 100, "ATG", "CTG") == f"{AC17}:g.100A>C"

    def test_lowercase_input(self):
        assert vv.vcf_to_hgvs_g(AC17, 5, "c", "t") == f"{AC17}:g.5C>T"


class FakeResolver:
    def __init__(self, answers):
        self.answers = answers
        self.calls = 0

    def g_to_c(self, chrom, pos, ref, alt, tx):
        self.calls += 1
        ans = self.answers[(chrom, pos, ref, alt, tx)]
        if isinstance(ans, Exception):
            raise ans
        return ans


class TestResolveOneLocal:
    def test_ok_strips_transcript_and_caches(self):
        r = FakeResolver({("17", 41197801, "T", "A", "NM_007294.4"): "NM_007294.4:c.5486A>T"})
        cache = {}
        assert vv.resolve_one_local("17", 41197801, "T", "A", "NM_007294.4", cache, r) == (
            "c.5486A>T", "ok", True)
        # second call is a cache hit, no resolver call, hit_api False
        assert vv.resolve_one_local("17", 41197801, "T", "A", "NM_007294.4", cache, r) == (
            "c.5486A>T", "ok", False)
        assert r.calls == 1

    def test_deterministic_flag_is_cached(self):
        r = FakeResolver({("17", 1, "A", "G", "NM_X.1"): ValueError("intergenic")})
        cache = {}
        assert vv.resolve_one_local("17", 1, "A", "G", "NM_X.1", cache, r) == (
            "", "flagged:intergenic", True)
        assert cache["17-1-A-G|NM_X.1"] == ("", "flagged:intergenic")

    def test_data_not_available_not_cached(self):
        # may be a transient SeqRepo fetch failure, so it is retried next run
        r = FakeResolver({("17", 1, "A", "G", "NM_X.1"): ValueError("tx_not_available")})
        cache = {}
        assert vv.resolve_one_local("17", 1, "A", "G", "NM_X.1", cache, r) == (
            "", "flagged:tx_not_available", True)
        assert cache == {}

    def test_connection_error_flagged_not_cached(self):
        r = FakeResolver({("17", 1, "A", "G", "NM_X.1"): ConnectionError("uta down")})
        cache = {}
        assert vv.resolve_one_local("17", 1, "A", "G", "NM_X.1", cache, r) == (
            "", "flagged:network_ConnectionError", True)
        assert cache == {}


class FakeResponse:
    def __init__(self, status_code, payload):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, resp):
        self.resp = resp
        self.urls = []

    def get(self, url, timeout):
        self.urls.append(url)
        return self.resp


class TestResolveOneRest:
    def test_rest_ok_unchanged(self):
        s = FakeSession(FakeResponse(200, {"NM_007294.4:c.5486A>T": {}, "flag": "gene_variant"}))
        out = vv.resolve_one("17", 41197801, "T", "A", "NM_007294.4", {}, s)
        assert out == ("c.5486A>T", "ok", True)
        assert s.urls[0].endswith("/GRCh37/17-41197801-T-A/NM_007294.4")

    def test_rest_http_error_flagged(self):
        s = FakeSession(FakeResponse(429, {}))
        assert vv.resolve_one("17", 1, "A", "G", "NM_X.1", {}, s) == ("", "flagged:http_429", True)


def test_annotate_file_with_pluggable_resolver(tmp_path):
    f = tmp_path / "GENE_seqnext.tsv"
    pd.DataFrame([
        {"chr_grch37": "chr17", "pos_grch37": "41197801", "ref": "T", "alt": "A",
         "transcript": "NM_007294.4", "hgvs_c": "c.9999A>T"},
        {"chr_grch37": "17", "pos_grch37": "5", "ref": "A", "alt": "G",
         "transcript": "NM_X.1", "hgvs_c": "c.1A>G (VV_FLAGGED:old)"},
    ]).to_csv(f, sep="\t", index=False)

    r = FakeResolver({
        ("17", 41197801, "T", "A", "NM_007294.4"): "NM_007294.4:c.5486A>T",
        ("17", 5, "A", "G", "NM_X.1"): ValueError("intergenic"),
    })

    def resolve(chrom, pos, ref, alt, tx, cache):
        return vv.resolve_one_local(chrom, pos, ref, alt, tx, cache, r)

    ok, flagged, err = vv.annotate_file(f, {}, resolve, in_place=True, rate_limit=0)
    assert (ok, flagged, err) == (1, 1, 0)
    out = pd.read_csv(f, sep="\t", dtype=str)
    assert out["hgvs_c"].tolist() == ["c.5486A>T", "c.1A>G (VV_FLAGGED:intergenic)"]
    assert out["vv_status"].tolist() == ["ok", "flagged:intergenic"]


@pytest.mark.integration
@pytest.mark.slow
def test_local_resolver_brca1_against_uta():
    """Needs the [hgvs] extra and network access to UTA/SeqRepo."""
    pytest.importorskip("hgvs")
    try:
        resolver = vv.LocalResolver()
    except Exception as e:  # noqa: BLE001 -- no UTA connectivity
        pytest.skip(f"UTA unavailable: {e}")
    assert resolver.g_to_c("17", 41197801, "T", "A", "NM_007294.4") == "NM_007294.4:c.5486A>T"

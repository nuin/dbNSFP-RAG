"""Safe FAISS metadata persistence and legacy migration."""

import copyreg
import json
import pickle

import pytest

from src.vector_metadata import load_metadata, save_metadata


def test_metadata_round_trip(tmp_path):
    metadata = {"ids": ["1_12_A_G"], "documents": ["variant"], "metas": [{"pos": 12}]}
    save_metadata(tmp_path / "metadata.json", metadata)
    assert load_metadata(tmp_path, 1) == metadata


@pytest.mark.parametrize("protocol", [2, 4, 5])
def test_legacy_metadata_remains_readable(tmp_path, protocol):
    metadata = {"ids": ["1_12_A_G"], "documents": ["variant"], "metas": [{"score": None}]}
    (tmp_path / "metadata.pkl").write_bytes(pickle.dumps(metadata, protocol=protocol))
    assert load_metadata(tmp_path, 1) == metadata


def test_legacy_metadata_rejects_executable_globals(tmp_path):
    (tmp_path / "metadata.pkl").write_bytes(b"cos\nsystem\n.")
    with pytest.raises(pickle.UnpicklingError, match="global"):
        load_metadata(tmp_path, 0)


def test_metadata_rejects_index_length_mismatch(tmp_path):
    (tmp_path / "metadata.json").write_text(json.dumps({"ids": [], "documents": [], "metas": []}))
    with pytest.raises(ValueError, match="index"):
        load_metadata(tmp_path, 1)


def test_legacy_extension_cannot_execute_cached_callable(tmp_path):
    invoked = []
    extension_code = 250
    payload = b"\x80\x02\x82\xfa)R."
    copyreg._extension_cache[extension_code] = lambda: invoked.append(True)
    try:
        (tmp_path / "metadata.pkl").write_bytes(payload)
        with pytest.raises(pickle.UnpicklingError):
            load_metadata(tmp_path, 0)
        assert invoked == []
    finally:
        copyreg._extension_cache.pop(extension_code)

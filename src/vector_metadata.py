"""JSON metadata storage with a data-only reader for legacy FAISS indexes."""

import io
import json
import logging
import pickle
import pickletools
from pathlib import Path

logger = logging.getLogger(__name__)

_DATA_OPCODES = frozenset(
    "PROTO FRAME STOP MARK POP POP_MARK DUP FLOAT BINFLOAT INT BININT BININT1 BININT2 "
    "LONG LONG1 LONG4 NONE NEWTRUE NEWFALSE STRING BINSTRING SHORT_BINSTRING UNICODE "
    "BINUNICODE SHORT_BINUNICODE BINUNICODE8 EMPTY_LIST APPEND APPENDS LIST EMPTY_TUPLE "
    "TUPLE TUPLE1 TUPLE2 TUPLE3 EMPTY_DICT DICT SETITEM SETITEMS PUT BINPUT LONG_BINPUT "
    "GET BINGET LONG_BINGET MEMOIZE".split()
)


class _DataOnlyUnpickler(pickle.Unpickler):
    """Legacy metadata contains only builtin containers and scalar values."""

    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"Metadata cannot reference global {module}.{name}")

    def persistent_load(self, pid):
        raise pickle.UnpicklingError("Metadata cannot reference persistent objects")


def _load_legacy(path: Path):
    data = path.read_bytes()
    # EXT opcodes can bypass find_class through the process-wide extension cache.
    for opcode, _argument, _position in pickletools.genops(data):
        if opcode.name not in _DATA_OPCODES:
            raise pickle.UnpicklingError(f"Metadata cannot use global/object opcode {opcode.name}")
    return _DataOnlyUnpickler(io.BytesIO(data)).load()


def load_metadata(directory: Path, expected_count: int) -> dict:
    path = directory / "metadata.json"
    if path.exists():
        with path.open() as stream:
            metadata = json.load(stream)
    else:
        legacy = directory / "metadata.pkl"
        logger.info("Reading legacy metadata from %s; next save writes JSON", legacy)
        metadata = _load_legacy(legacy)
    if not isinstance(metadata, dict):
        raise ValueError("Vector metadata must be an object")
    for key, item_type in (("ids", str), ("documents", str), ("metas", dict)):
        values = metadata.get(key)
        if not isinstance(values, list) or len(values) != expected_count:
            raise ValueError(f"Metadata {key} must match the FAISS index length")
        if not all(isinstance(value, item_type) for value in values):
            raise ValueError(f"Metadata {key} contains invalid values")
    return metadata


def save_metadata(path: Path, metadata: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    with temporary.open("w") as stream:
        json.dump(metadata, stream)
    temporary.replace(path)

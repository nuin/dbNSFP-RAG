"""SQLite-backed variant database for fast coordinate lookups.

Replaces FAISS vector store for production use. All lookups are indexed
O(1) or O(log n) instead of O(n) linear scans. No ML dependencies required.
"""

import json
import sqlite3
from pathlib import Path
from typing import Optional

from .config import DATA_DIR


# Default database path
DEFAULT_SQLITE_PATH = DATA_DIR / "sqlite" / "grch37-all-panels.db"


class VariantDatabase:
    """SQLite-backed variant store matching VariantVectorStore's public API.

    Returns the same ``{"id": str, "document": str, "metadata": dict}`` format
    so all downstream code (ACMG scoring, API, export) works unchanged.
    """

    SCHEMA = """
    CREATE TABLE IF NOT EXISTS variants (
        variant_id TEXT PRIMARY KEY,
        chr TEXT NOT NULL,
        pos INTEGER NOT NULL,
        ref TEXT NOT NULL,
        alt TEXT NOT NULL,
        gene TEXT,
        transcript TEXT,
        aa_ref TEXT,
        aa_alt TEXT,
        aa_pos INTEGER,
        cadd_phred REAL,
        revel_score REAL,
        clinvar_id TEXT,
        clinvar_sig TEXT,
        clinvar_review TEXT,
        clinvar_trait TEXT,
        gnomad_af REAL,
        gnomad_hom REAL,
        sift_pred TEXT,
        polyphen_pred TEXT,
        alphamissense_pred TEXT,
        mutationtaster_pred TEXT,
        bayesdel_pred TEXT,
        provean_pred TEXT,
        consequence TEXT,
        gnomad_pli REAL,
        gnomad_mis_oe REAL,
        loeuf REAL,
        interpro_domain TEXT,
        full_annotation TEXT
    );

    CREATE TABLE IF NOT EXISTS variant_panels (
        variant_id TEXT NOT NULL,
        panel_name TEXT NOT NULL,
        PRIMARY KEY (variant_id, panel_name),
        FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
    );

    CREATE TABLE IF NOT EXISTS db_metadata (
        key TEXT PRIMARY KEY,
        value TEXT
    );
    """

    INDEXES = """
    CREATE INDEX IF NOT EXISTS idx_coordinates ON variants(chr, pos, ref, alt);
    CREATE INDEX IF NOT EXISTS idx_gene ON variants(gene);
    CREATE INDEX IF NOT EXISTS idx_clinvar_sig ON variants(clinvar_sig);
    CREATE INDEX IF NOT EXISTS idx_chr_pos ON variants(chr, pos);
    CREATE INDEX IF NOT EXISTS idx_variant_panels_panel ON variant_panels(panel_name);
    """

    # Column names that map directly from metadata dict to DB columns
    METADATA_COLUMNS = [
        "chr", "pos", "ref", "alt", "gene", "transcript",
        "aa_ref", "aa_alt", "aa_pos",
        "cadd_phred", "revel_score",
        "clinvar_id", "clinvar_sig", "clinvar_review", "clinvar_trait",
        "gnomad_af", "gnomad_hom",
        "sift_pred", "polyphen_pred", "alphamissense_pred",
        "mutationtaster_pred", "bayesdel_pred", "provean_pred",
        "consequence", "gnomad_pli", "gnomad_mis_oe", "loeuf",
        "interpro_domain",
    ]

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = Path(db_path or DEFAULT_SQLITE_PATH)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        self._conn = sqlite3.connect(
            str(self.db_path),
            check_same_thread=False,
        )
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.execute("PRAGMA cache_size=-64000")  # 64 MB cache

        self._init_schema()

    def _init_schema(self):
        """Create tables and indexes if they don't exist."""
        self._conn.executescript(self.SCHEMA)
        self._conn.executescript(self.INDEXES)
        self._conn.commit()

    def close(self):
        """Close the database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def __del__(self):
        self.close()

    # ------------------------------------------------------------------
    # Metadata helpers
    # ------------------------------------------------------------------

    def get_meta(self, key: str) -> Optional[str]:
        """Get a value from db_metadata."""
        row = self._conn.execute(
            "SELECT value FROM db_metadata WHERE key = ?", (key,)
        ).fetchone()
        return row["value"] if row else None

    def set_meta(self, key: str, value: str):
        """Set a value in db_metadata."""
        self._conn.execute(
            "INSERT OR REPLACE INTO db_metadata (key, value) VALUES (?, ?)",
            (key, value),
        )
        self._conn.commit()

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def add_variants(
        self,
        ids: list[str],
        texts: list[str],
        metadatas: list[dict],
        show_progress: bool = True,
        batch_size: int = 1000,
    ):
        """Insert variants into the database.

        Matches VariantVectorStore.add_variants() signature.
        The ``batch_size`` param here controls SQLite commit batching
        (not embedding batches — there are no embeddings).
        """
        col_placeholders = ", ".join(["?"] * (len(self.METADATA_COLUMNS) + 2))
        col_names = ", ".join(["variant_id"] + self.METADATA_COLUMNS + ["full_annotation"])

        sql = f"INSERT OR REPLACE INTO variants ({col_names}) VALUES ({col_placeholders})"

        rows = []
        for var_id, text, meta in zip(ids, texts, metadatas):
            values = [var_id]
            for col in self.METADATA_COLUMNS:
                values.append(meta.get(col))
            values.append(text)
            rows.append(values)

        # Batch insert
        cursor = self._conn.cursor()
        for i in range(0, len(rows), batch_size):
            batch = rows[i : i + batch_size]
            cursor.executemany(sql, batch)
        self._conn.commit()

    def add_panel_membership(self, variant_ids: list[str], panel_name: str):
        """Record which panel(s) a set of variants belong to."""
        sql = "INSERT OR IGNORE INTO variant_panels (variant_id, panel_name) VALUES (?, ?)"
        self._conn.executemany(sql, [(vid, panel_name) for vid in variant_ids])
        self._conn.commit()

    def save(self):
        """Explicit commit (for API compatibility with VariantVectorStore)."""
        self._conn.commit()

    def delete_all(self):
        """Delete all variants."""
        self._conn.execute("DELETE FROM variant_panels")
        self._conn.execute("DELETE FROM variants")
        self._conn.execute("DELETE FROM db_metadata")
        self._conn.commit()
        print("Deleted all variants")

    def optimize(self):
        """Run VACUUM and ANALYZE for optimal query performance."""
        self._conn.execute("ANALYZE")
        self._conn.commit()
        # VACUUM cannot run inside a transaction
        self._conn.execute("VACUUM")

    # ------------------------------------------------------------------
    # Read operations (same return format as VariantVectorStore)
    # ------------------------------------------------------------------

    def _row_to_result(self, row: sqlite3.Row) -> dict:
        """Convert a DB row to the standard result dict."""
        meta = {}
        for col in self.METADATA_COLUMNS:
            val = row[col]
            if val is not None:
                meta[col] = val
            else:
                meta[col] = None

        return {
            "id": row["variant_id"],
            "document": row["full_annotation"] or "",
            "metadata": meta,
        }

    def get_by_id(self, variant_id: str) -> Optional[dict]:
        """Get a variant by its exact ID.  O(1) via primary key."""
        row = self._conn.execute(
            "SELECT * FROM variants WHERE variant_id = ?", (variant_id,)
        ).fetchone()
        if row is None:
            return None
        return self._row_to_result(row)

    def get_by_coordinates(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Optional[dict]:
        """Get a variant by genomic coordinates.  Indexed lookup."""
        chrom = str(chrom).replace("chr", "")
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        return self.get_by_id(variant_id)

    def search_by_gene(self, gene: str, k: int = 100) -> list[dict]:
        """Get variants for a gene.  Indexed query, ordered by position."""
        rows = self._conn.execute(
            "SELECT * FROM variants WHERE UPPER(gene) = ? ORDER BY pos LIMIT ?",
            (gene.upper(), k),
        ).fetchall()
        return [self._row_to_result(r) for r in rows]

    def search_by_region(
        self, chrom: str, start: int, end: int, k: int = 1000
    ) -> list[dict]:
        """Get variants in a genomic region.  Indexed range query."""
        chrom = str(chrom).replace("chr", "")
        rows = self._conn.execute(
            "SELECT * FROM variants WHERE chr = ? AND pos BETWEEN ? AND ? ORDER BY pos LIMIT ?",
            (chrom, start, end, k),
        ).fetchall()
        return [self._row_to_result(r) for r in rows]

    def search_by_panel(self, panel_name: str, k: int = 10000) -> list[dict]:
        """Get all variants belonging to a specific panel."""
        rows = self._conn.execute(
            """SELECT v.* FROM variants v
               JOIN variant_panels vp ON v.variant_id = vp.variant_id
               WHERE vp.panel_name = ?
               LIMIT ?""",
            (panel_name, k),
        ).fetchall()
        return [self._row_to_result(r) for r in rows]

    def list_panels(self) -> list[str]:
        """Return distinct panel names stored in the database."""
        rows = self._conn.execute(
            "SELECT DISTINCT panel_name FROM variant_panels ORDER BY panel_name"
        ).fetchall()
        return [r["panel_name"] for r in rows]

    def count(self) -> int:
        """Return total number of variants."""
        row = self._conn.execute("SELECT COUNT(*) AS cnt FROM variants").fetchone()
        return row["cnt"]

    def count_by_panel(self, panel_name: str) -> int:
        """Return number of variants in a specific panel."""
        row = self._conn.execute(
            "SELECT COUNT(*) AS cnt FROM variant_panels WHERE panel_name = ?",
            (panel_name,),
        ).fetchone()
        return row["cnt"]

    # ------------------------------------------------------------------
    # Browse / navigation helpers
    # ------------------------------------------------------------------

    # Natural chromosome sort order
    _CHR_ORDER = {str(i): i for i in range(1, 23)}
    _CHR_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})

    @staticmethod
    def _normalize_chr(raw: str) -> str:
        """Normalize chromosome strings like '2.0' → '2'."""
        raw = str(raw).replace("chr", "")
        # Strip trailing '.0' from float-casted values
        if raw.endswith(".0"):
            raw = raw[:-2]
        return raw

    def list_chromosomes(self) -> list[dict]:
        """Return distinct chromosomes with variant counts, naturally sorted.

        Merges float-casted duplicates (e.g. '2' + '2.0') into a single entry.
        """
        rows = self._conn.execute(
            "SELECT chr, COUNT(*) AS count FROM variants GROUP BY chr"
        ).fetchall()
        merged: dict[str, int] = {}
        for r in rows:
            norm = self._normalize_chr(r["chr"])
            merged[norm] = merged.get(norm, 0) + r["count"]
        results = [{"chr": k, "count": v} for k, v in merged.items()]
        results.sort(key=lambda x: self._CHR_ORDER.get(x["chr"], 99))
        return results

    def list_genes_by_chromosome(self, chrom: str) -> list[dict]:
        """Return genes on a chromosome with variant counts, sorted by count desc.

        Queries both '2' and '2.0' forms to capture all variants.
        """
        chrom = self._normalize_chr(chrom)
        float_form = f"{chrom}.0"
        rows = self._conn.execute(
            "SELECT gene, COUNT(*) AS count FROM variants "
            "WHERE (chr = ? OR chr = ?) AND gene IS NOT NULL "
            "GROUP BY gene ORDER BY COUNT(*) DESC",
            (chrom, float_form),
        ).fetchall()
        return [{"gene": r["gene"], "count": r["count"]} for r in rows]

    # ------------------------------------------------------------------
    # Checkpoint helpers (for build resume)
    # ------------------------------------------------------------------

    def get_completed_chromosomes(self) -> list[str]:
        """Get list of completed chromosomes from db_metadata."""
        val = self.get_meta("completed_chromosomes")
        if val:
            return json.loads(val)
        return []

    def set_completed_chromosomes(self, chroms: list[str]):
        """Save completed chromosomes to db_metadata."""
        self.set_meta("completed_chromosomes", json.dumps(chroms))

    def get_total_variants_processed(self) -> int:
        """Get total variants processed count from db_metadata."""
        val = self.get_meta("total_variants_processed")
        return int(val) if val else 0

    def set_total_variants_processed(self, count: int):
        """Save total variants processed count."""
        self.set_meta("total_variants_processed", str(count))

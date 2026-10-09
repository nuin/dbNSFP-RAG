"""Configure SQLite and FAISS builds without loading unused ML backends."""

import logging
from pathlib import Path

from .build_pipeline import BuildOptions, run_build
from .config import SQLITE_DB_PATH, VECTORDB_DIR, VECTORDB_GRCH37_NGSGENES
from .gene_data import load_gene_data
from .panels import PANELS, get_all_panel_genes, get_panel

logger = logging.getLogger(__name__)


def build_index(source: Path, chromosomes=None, limit=None, resume=True):
    from .vectorstore import VariantVectorStore

    store = VariantVectorStore()
    return run_build(store, BuildOptions(source, chromosomes, limit, resume))


def build_panel_index(
    source: Path, panel_name: str, use_grch37=True, chromosomes=None, limit=None, resume=True
):
    from .vectorstore import VariantVectorStore

    try:
        genes = get_panel(panel_name)
    except ValueError:
        logger.exception("Unknown gene panel %s", panel_name)
        return 1
    if use_grch37 and panel_name.lower() == "ngsgenes":
        path = VECTORDB_GRCH37_NGSGENES
    else:
        path = VECTORDB_DIR / f"{panel_name.lower()}-{'grch37' if use_grch37 else 'grch38'}"
    options = BuildOptions(source, chromosomes, limit, resume, use_grch37, genes, load_gene_data())
    return run_build(VariantVectorStore(db_path=path), options)


def _populate_panels(db):
    for name, genes in PANELS.items():
        ids = [result["id"] for gene in genes for result in db.search_by_gene(gene, k=1000000)]
        if ids:
            db.add_panel_membership(ids, name)
            logger.info("%s: %s variants", name, len(ids))


def build_sqlite_index(
    source: Path, use_grch37=True, chromosomes=None, limit=None, resume=True, db_path=None
):
    from .variantdb import VariantDatabase

    db = VariantDatabase(db_path=db_path or SQLITE_DB_PATH)
    try:
        options = BuildOptions(
            source, chromosomes, limit, resume, use_grch37, get_all_panel_genes(), load_gene_data()
        )
        status = run_build(db, options)
        if status:
            return status
        _populate_panels(db)
        db.optimize()
        logger.info("SQLite build complete: %s variants in %s", db.count(), db.db_path)
        return 0
    finally:
        db.close()

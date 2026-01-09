"""Configuration for dbNSFP RAG pipeline."""

from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
VECTORDB_DIR = DATA_DIR / "vectordb"

# dbNSFP zip file (stream directly without extracting)
DBNSFP_ZIP = Path.home() / "dbNSFP5.3.1a.zip"

# dbNSFP columns to keep (from 600+)
KEEP_COLUMNS = [
    "#chr",
    "pos(1-based)",
    "ref",
    "alt",
    "aaref",
    "aaalt",
    "aapos",
    "genename",
    "Ensembl_geneid",
    "Ensembl_transcriptid",
    "Ensembl_proteinid",
    "Uniprot_acc",
    "HGVSc_snpEff",
    "HGVSp_snpEff",
    # SIFT
    "SIFT_score",
    "SIFT_pred",
    # PolyPhen2
    "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_score",
    "Polyphen2_HVAR_pred",
    # CADD
    "CADD_raw",
    "CADD_phred",
    # REVEL
    "REVEL_score",
    # AlphaMissense
    "AlphaMissense_score",
    "AlphaMissense_pred",
    # ClinPred
    "ClinPred_score",
    "ClinPred_pred",
    # DANN
    "DANN_score",
    # MetaSVM/MetaLR
    "MetaSVM_score",
    "MetaSVM_pred",
    "MetaLR_score",
    "MetaLR_pred",
    # Conservation
    "phyloP100way_vertebrate",
    "phyloP30way_mammalian",
    "phastCons100way_vertebrate",
    "GERP++_RS",
    # Population frequencies
    "gnomAD_exomes_AF",
    "gnomAD_genomes_AF",
    "1000Gp3_AF",
    # ClinVar
    "clinvar_id",
    "clinvar_clnsig",
    "clinvar_review",
    "clinvar_trait",
    # Functional
    "Interpro_domain",
    "GTEx_V8_gene",
    "GTEx_V8_tissue",
]

# Embedding model (all-MiniLM-L6-v2 is 5x faster than PubMedBERT)
EMBEDDING_MODEL = "all-MiniLM-L6-v2"
EMBEDDING_DIMENSION = 384

# ChromaDB
COLLECTION_NAME = "dbnsfp_variants"

# Processing
CHUNK_SIZE = 10000  # Rows per batch when parsing
EMBEDDING_BATCH_SIZE = 1024  # Variants per embedding batch (larger for GPU)

# LLM
OLLAMA_URL = "http://localhost:11434"
DEFAULT_MODEL = "llama3.1:8b"

"""Configuration for dbNSFP RAG pipeline."""

import os
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
VECTORDB_DIR = DATA_DIR / "vectordb"

# dbNSFP data location (extracted folder or zip file)
# Can be overridden with DBNSFP_DIR environment variable
_default_dbnsfp = Path.home() / "Downloads" / "dbNSFP5.3.1a"
DBNSFP_DIR = Path(os.environ.get("DBNSFP_DIR", str(_default_dbnsfp)))
DBNSFP_ZIP = None  # Not using zip file

# Panel database paths by genome build
VECTORDB_GRCH37_NGSGENES = DATA_DIR / "vectordb" / "grch37-ngsgenes"
VECTORDB_GRCH38_NGSGENES = DATA_DIR / "vectordb" / "grch38-ngsgenes"

# dbNSFP columns to keep (from 600+)
KEEP_COLUMNS = [
    "#chr",
    "pos(1-based)",
    "hg19_chr",
    "hg19_pos(1-based)",
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
    # MutationTaster
    "MutationTaster_score",
    "MutationTaster_pred",
    # BayesDel
    "BayesDel_addAF_score",
    "BayesDel_addAF_pred",
    # PROVEAN
    "PROVEAN_score",
    "PROVEAN_pred",
    # Conservation
    "phyloP100way_vertebrate",
    "phyloP470way_mammalian",  # Updated from phyloP30way in dbNSFP 5.x
    "phastCons100way_vertebrate",
    "GERP++_RS",
    # Population frequencies (dbNSFP 5.x column names)
    "gnomAD4.1_joint_AF",  # gnomAD v4.1 combined exomes+genomes
    "gnomAD2.1.1_exomes_controls_AF",  # gnomAD v2 exomes (controls)
    "1000Gp3_AF",
    # ClinVar
    "clinvar_id",
    "clinvar_clnsig",
    "clinvar_review",
    "clinvar_trait",
    # Functional
    "Interpro_domain",
    # Gene constraint (gnomAD) - for PP2, BP1, PVS1 evaluation
    "gnomAD_pLI",  # Probability of LoF intolerance (≥0.9 = intolerant)
    "gnomAD_mis_z",  # Missense Z-score (≥3.09 = constrained)
    "gnomAD_lof_z",  # LoF Z-score (higher = more constrained)
    "LOEUF",  # LoF observed/expected upper bound (v4: <0.6 = constrained)
    # Consequence type (for PVS1, PM4, PP2, BP1, BP7)
    "Ensembl_consequence",  # VEP consequence (frameshift, stop_gained, synonymous, etc.)
    # SpliceAI scores (for BP7 - synonymous splice impact evaluation)
    "SpliceAI_pred_DS_AG",  # Acceptor gain delta score
    "SpliceAI_pred_DS_AL",  # Acceptor loss delta score
    "SpliceAI_pred_DS_DG",  # Donor gain delta score
    "SpliceAI_pred_DS_DL",  # Donor loss delta score
    # gnomAD homozygote count (for BS2 - observed in healthy adults)
    "gnomAD4.1_joint_AC_hom",  # Homozygote allele count in gnomAD v4.1
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

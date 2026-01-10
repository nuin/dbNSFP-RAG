"""Gene panel management for targeted variant annotation."""

from pathlib import Path

# Define your gene panels here
# Each panel is a set of gene symbols

PANELS = {
    # Example: Hereditary cancer panel
    "hereditary_cancer": {
        "BRCA1", "BRCA2", "TP53", "PTEN", "CDH1", "STK11", "CHEK2",
        "PALB2", "ATM", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",
        "APC", "MUTYH", "BMPR1A", "SMAD4", "RAD51C", "RAD51D",
        "BRIP1", "NBN", "BARD1", "CDK4", "CDKN2A",
    },

    # Example: Cardiac panel
    "cardiac": {
        "MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3",
        "ACTC1", "LMNA", "SCN5A", "KCNQ1", "KCNH2", "RYR2", "PKP2",
        "DSP", "DSC2", "DSG2", "TMEM43", "PLN",
    },

    # Example: Neurological panel
    "neurological": {
        "SOD1", "FUS", "TARDBP", "C9orf72", "MAPT", "GRN", "APP",
        "PSEN1", "PSEN2", "HTT", "PARK2", "PINK1", "LRRK2", "SNCA",
    },

    # NGS panel (~301 genes) - cardiovascular and cancer plus combined
    "NGSgenes": {
        "ABCA1", "ABCC9", "ABCG5", "ABCG8", "ACTA1", "ACTA2", "ACTC1", "ACTN2",
        "ACVRL1", "AIP", "AKAP9", "ALK", "ALMS1", "ALPK3", "ANGPTL4", "ANK2",
        "ANKRD1", "APC", "APOA4", "APOA5", "APOB", "APOC2", "APOC3", "APOE",
        "ATM", "BAG3", "BAP1", "BARD1", "BLM", "BMPR1A", "BRAF", "BRCA1",
        "BRCA2", "BRIP1", "BUB1B", "CACNA1C", "CACNA1S", "CACNA2D1", "CACNB2",
        "CALM1", "CALM2", "CALM3", "CALR3", "CASQ2", "CASR", "CAV3", "CAVIN4",
        "CBL", "CBS", "CDC73", "CDH1", "CDK4", "CDKN1B", "CDKN1C", "CDKN2A",
        "CEBPA", "CEP57", "CETP", "CHEK2", "COL1A1", "COL3A1", "COL5A1",
        "COL5A2", "COX15", "CREB3L3", "CRELD1", "CRYAB", "CSRP3", "CTF1",
        "CYLD", "CYP7A1", "DDB2", "DES", "DICER1", "DIS3L2", "DMD", "DNAJC19",
        "DOLK", "DPP6", "DSC2", "DSG2", "DSP", "DTNA", "EFEMP2", "EGFR", "ELN",
        "EMD", "ENG", "EPCAM", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "EXT1",
        "EXT2", "EYA4", "EZH2", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE",
        "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "FBN1", "FBN2", "FH",
        "FHL1", "FHL2", "FKRP", "FKTN", "FLCN", "FLNC", "FXN", "GAA", "GATA2",
        "GATAD1", "GCKR", "GJA5", "GLA", "GNAS", "GPC3", "GPD1L", "GPIHBP1",
        "HADHA", "HCN4", "HFE", "HNF1A", "HRAS", "HSPB8", "ILK", "JAG1", "JPH2",
        "JUP", "KCNA5", "KCND3", "KCNE1", "KCNE2", "KCNE3", "KCNH2", "KCNJ2",
        "KCNJ5", "KCNJ8", "KCNQ1", "KIT", "KLF10", "KRAS", "LAMA2", "LAMA4",
        "LAMP2", "LCAT", "LDB3", "LDLR", "LDLRAP1", "LEMD2", "LIPA", "LIPC",
        "LIPG", "LMF1", "LMNA", "LPL", "LTBP2", "MAP2K1", "MAP2K2", "MAT2A",
        "MAX", "MEN1", "MET", "MIB1", "MLH1", "MRE11", "MRE11A", "MSH2", "MSH6",
        "MURC", "MUTYH", "MYBPC3", "MYH11", "MYH6", "MYH7", "MYL2", "MYL3",
        "MYLK", "MYLK2", "MYO6", "MYOZ2", "MYPN", "NBN", "NEXN", "NF1", "NF2",
        "NKX2", "NKX2-1", "NKX2-2", "NKX2-3", "NKX2-4", "NKX2-5", "NKX2-6",
        "NKX2-8", "NKX4", "NODAL", "NOS1AP", "NOTCH1", "NPPA", "NRAS", "NSD1",
        "PALB2", "PCSK9", "PDE4D", "PDLIM3", "PHOX2B", "PKP2", "PLN", "PMS1",
        "PMS2", "PMS2ALL", "POLD1", "POLE", "PPA2", "PRDM16", "PRF1", "PRKAG2",
        "PRKAR1A", "PRKG1", "PSEN1", "PSEN2", "PTCH1", "PTEN", "PTPN11",
        "RAD50", "RAD51C", "RAD51D", "RAF1", "RANGRF", "RB1", "RBM20", "RECQL4",
        "RET", "RHBDF2", "RIT1", "RUNX1", "RYR1", "RYR2", "SALL4", "SBDS",
        "SCARB1", "SCN10A", "SCN1B", "SCN2B", "SCN3B", "SCN4B", "SCN5A", "SCO2",
        "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEPN1", "SGCB", "SGCD",
        "SGCG", "SHOC2", "SLC25A4", "SLC2A10", "SLMAP", "SLX4", "SMAD3",
        "SMAD4", "SMARCB1", "SNTA1", "SOS1", "SREBF2", "STAP1", "STK11", "SUFU",
        "TAZ", "TBX20", "TBX3", "TBX5", "TCAP", "TECRL", "TGFB2", "TGFB3",
        "TGFBR1", "TGFBR2", "TMEM127", "TMEM43", "TMPO", "TNNC1", "TNNI3",
        "TNNT2", "TP53", "TPM1", "TRDN", "TRIM63", "TRPM4", "TSC1", "TSC2",
        "TTN", "TTR", "TXNRD2", "VCL", "VHL", "WRN", "WT1", "XPA", "XPC",
        "ZBTB17", "ZHX3", "ZIC3",
    },
}


def get_panel(name: str) -> set[str]:
    """Get genes for a specific panel."""
    if name not in PANELS:
        available = ", ".join(PANELS.keys())
        raise ValueError(f"Unknown panel: {name}. Available: {available}")
    return PANELS[name]


def get_all_panel_genes() -> set[str]:
    """Get all genes across all panels."""
    all_genes = set()
    for genes in PANELS.values():
        all_genes.update(genes)
    return all_genes


def list_panels() -> dict[str, int]:
    """List all panels with gene counts."""
    return {name: len(genes) for name, genes in PANELS.items()}


def load_panel_from_file(filepath: Path) -> set[str]:
    """
    Load gene list from a file (one gene per line).

    Args:
        filepath: Path to gene list file

    Returns:
        Set of gene symbols
    """
    genes = set()
    with open(filepath) as f:
        for line in f:
            gene = line.strip().upper()
            if gene and not gene.startswith("#"):
                genes.add(gene)
    return genes


def add_panel(name: str, genes: set[str]) -> None:
    """Add a new panel dynamically."""
    PANELS[name] = genes


if __name__ == "__main__":
    print("Available panels:")
    for name, count in list_panels().items():
        print(f"  {name}: {count} genes")

    print(f"\nTotal unique genes: {len(get_all_panel_genes())}")

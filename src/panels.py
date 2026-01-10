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

    # Add your custom panels below
    # "your_panel_name": {
    #     "GENE1", "GENE2", "GENE3", ...
    # },
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

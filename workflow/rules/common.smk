def get_peptide_sequence(wildcards):
    """Retrieve the amino acid sequence or file path for a given sample wildcard."""
    return config["peptides"][wildcards.sample]

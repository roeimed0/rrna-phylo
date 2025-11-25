"""
Configuration module for rRNA-Phylo default parameters.

This module centralizes all default configuration values used throughout
the phylogenetic pipeline, making it easy to adjust defaults in one place.
"""

# Visualization defaults
VISUALIZATION_DEFAULTS = {
    'dpi': 300,
    'branch_line_width': 1.5,
    'scale': 120,
    'branch_vertical_margin': 10,
    'bootstrap_threshold': 70.0,
}

# Tree building defaults
TREE_DEFAULTS = {
    'gamma_alpha': 1.0,  # Gamma shape parameter
    'bootstrap_replicates': 100,
}

# Distance calculation defaults
DISTANCE_DEFAULTS = {
    'dna_model': 'jukes-cantor',
    'protein_model': 'poisson',
}

# Output defaults
OUTPUT_DEFAULTS = {
    'output_format': 'both',  # ascii, newick, or both
    'show_tree_ascii': True,
}

# Alignment defaults
ALIGNMENT_DEFAULTS = {
    'muscle_threads': 1,
}

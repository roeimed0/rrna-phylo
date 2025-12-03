"""
Utility functions for rRNA phylogenetic analysis.

This package now serves as a compatibility layer. Most functionality has been
moved to more appropriate packages:
- preprocessing/ - Data preparation (strain handling, sampling, outgroups)
- analysis/ - Analysis tools (bootstrap, dataset analysis)
- visualization/ - Tree visualization

Import from those packages directly for new code.
"""

# Backward compatibility imports - deprecated, use new packages
from rrna_phylo.preprocessing import (
    remove_exact_duplicates,
    smart_dereplicate,
    dereplicate_strains,
    get_strain_summary,
    stratified_sample,
    detect_bias,
    suggest_outgroup,
    get_outgroup_sequences
)

from rrna_phylo.analysis import (
    calculate_bootstrap_support,
    bootstrap_tree,
    bootstrap_consensus,
    analyze_distance_matrix,
    recommend_methods
)

from rrna_phylo.visualization import (
    print_tree_ascii,
    visualize_all_trees
)

# Actual utilities
from .console import *

__all__ = [
    # Deprecated - use rrna_phylo.preprocessing
    'remove_exact_duplicates',
    'smart_dereplicate',
    'dereplicate_strains',
    'get_strain_summary',
    'stratified_sample',
    'detect_bias',
    'suggest_outgroup',
    'get_outgroup_sequences',
    # Deprecated - use rrna_phylo.analysis
    'calculate_bootstrap_support',
    'bootstrap_tree',
    'bootstrap_consensus',
    'analyze_distance_matrix',
    'recommend_methods',
    # Deprecated - use rrna_phylo.visualization
    'print_tree_ascii',
    'visualize_all_trees',
]

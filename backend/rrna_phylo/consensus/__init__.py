"""
Consensus tree methods for combining multiple phylogenetic trees.

This module provides:
- Tree comparison metrics (Robinson-Foulds distance)
- Bipartition extraction and manipulation
- Consensus tree construction (majority-rule, strict)
- Support value calculation
"""

from rrna_phylo.consensus.tree_distance import (
    robinson_foulds_distance,
    informative_rf_distance,
    tree_similarity,
    compare_trees
)
from rrna_phylo.consensus.bipartitions import (
    get_bipartitions,
    get_informative_bipartitions,
    are_compatible,
    get_taxa_from_trees
)

__all__ = [
    # Tree distance metrics
    "robinson_foulds_distance",
    "informative_rf_distance",
    "tree_similarity",
    "compare_trees",

    # Bipartition utilities
    "get_bipartitions",
    "get_informative_bipartitions",
    "are_compatible",
    "get_taxa_from_trees",
]

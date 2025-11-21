"""
Tree distance metrics for comparing phylogenetic trees.

Includes Robinson-Foulds distance and related metrics.
"""

from typing import Tuple
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.consensus.bipartitions import (
    get_bipartitions,
    get_informative_bipartitions,
    get_leaf_names
)


def robinson_foulds_distance(tree1: TreeNode, tree2: TreeNode, normalize: bool = False) -> float:
    """
    Calculate Robinson-Foulds (RF) distance between two trees.

    The RF distance is the number of bipartitions that differ between trees.
    RF = |splits in tree1 not in tree2| + |splits in tree2 not in tree1|

    Args:
        tree1: First tree
        tree2: Second tree
        normalize: If True, return normalized RF distance (0-1 scale)

    Returns:
        RF distance (integer if normalize=False, float 0-1 if normalize=True)

    Raises:
        ValueError: If trees have different taxa

    Example:
        tree1 = parse_newick("((A,B),(C,D));")
        tree2 = parse_newick("((A,C),(B,D));")
        rf = robinson_foulds_distance(tree1, tree2)
        # rf > 0 because topologies differ
    """
    # Get bipartitions from both trees
    bp1 = get_bipartitions(tree1)
    bp2 = get_bipartitions(tree2)

    # Verify same taxa
    taxa1 = get_leaf_names(tree1)
    taxa2 = get_leaf_names(tree2)

    if taxa1 != taxa2:
        raise ValueError(
            f"Trees have different taxa. "
            f"Tree1: {sorted(taxa1)}, Tree2: {sorted(taxa2)}"
        )

    # Calculate symmetric difference
    rf_distance = len(bp1.symmetric_difference(bp2))

    if normalize:
        # Maximum possible RF distance for n taxa
        # For unrooted binary tree: max_rf = 2(n-3)
        n = len(taxa1)
        if n < 3:
            return 0.0
        max_rf = 2 * (n - 3)
        return rf_distance / max_rf if max_rf > 0 else 0.0

    return rf_distance


def informative_rf_distance(tree1: TreeNode, tree2: TreeNode, normalize: bool = False) -> float:
    """
    Calculate RF distance using only informative bipartitions.

    Excludes trivial splits (single taxon vs rest), focusing on
    internal tree structure.

    Args:
        tree1: First tree
        tree2: Second tree
        normalize: If True, return normalized distance (0-1)

    Returns:
        RF distance based on informative bipartitions only
    """
    # Get informative bipartitions (exclude trivial splits)
    bp1 = get_informative_bipartitions(tree1)
    bp2 = get_informative_bipartitions(tree2)

    # Calculate symmetric difference
    rf_distance = len(bp1.symmetric_difference(bp2))

    if normalize:
        n = len(get_leaf_names(tree1))
        if n < 4:
            return 0.0
        # Max informative splits for n taxa
        max_rf = 2 * (n - 3)
        return rf_distance / max_rf if max_rf > 0 else 0.0

    return rf_distance


def tree_similarity(tree1: TreeNode, tree2: TreeNode) -> float:
    """
    Calculate similarity between trees (inverse of normalized RF).

    Args:
        tree1: First tree
        tree2: Second tree

    Returns:
        Similarity score from 0.0 (completely different) to 1.0 (identical)
    """
    rf_norm = robinson_foulds_distance(tree1, tree2, normalize=True)
    return 1.0 - rf_norm


def compare_trees(tree1: TreeNode, tree2: TreeNode) -> dict:
    """
    Comprehensive comparison of two trees.

    Args:
        tree1: First tree
        tree2: Second tree

    Returns:
        Dictionary with multiple distance metrics:
        - 'rf_distance': Raw RF distance
        - 'rf_normalized': Normalized RF (0-1)
        - 'similarity': Tree similarity (0-1)
        - 'shared_splits': Number of shared bipartitions
        - 'total_splits': Total unique bipartitions
        - 'identical': Boolean, True if trees are identical

    Example:
        result = compare_trees(upgma_tree, bionj_tree)
        print(f"Similarity: {result['similarity']:.2%}")
        print(f"Shared splits: {result['shared_splits']}/{result['total_splits']}")
    """
    bp1 = get_bipartitions(tree1)
    bp2 = get_bipartitions(tree2)

    shared = bp1.intersection(bp2)
    total = bp1.union(bp2)

    rf = len(bp1.symmetric_difference(bp2))
    n = len(get_leaf_names(tree1))
    max_rf = 2 * (n - 3) if n >= 3 else 0

    return {
        'rf_distance': rf,
        'rf_normalized': rf / max_rf if max_rf > 0 else 0.0,
        'similarity': 1.0 - (rf / max_rf if max_rf > 0 else 0.0),
        'shared_splits': len(shared),
        'total_splits': len(total),
        'identical': rf == 0,
        'num_taxa': n
    }

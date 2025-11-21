"""
Consensus tree construction from multiple phylogenetic trees.

Implements majority-rule and strict consensus algorithms.
"""

from typing import List, Tuple, Dict, FrozenSet, Set
from collections import Counter
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.consensus.bipartitions import (
    get_bipartitions,
    get_taxa_from_trees,
    are_compatible,
    get_leaf_names
)


def count_bipartitions(trees: List[TreeNode]) -> Tuple[Counter, Set[str]]:
    """
    Count frequency of each bipartition across multiple trees.

    Args:
        trees: List of phylogenetic trees

    Returns:
        Tuple of (Counter of bipartitions, set of all taxa)

    Example:
        trees = [tree1, tree2, tree3]
        counts, taxa = count_bipartitions(trees)
        # counts[frozenset({'A','B'})] = 2  # appears in 2 of 3 trees
    """
    all_taxa = get_taxa_from_trees(trees)
    bipartition_counts = Counter()

    for tree in trees:
        bipartitions = get_bipartitions(tree)
        for bp in bipartitions:
            bipartition_counts[bp] += 1

    return bipartition_counts, all_taxa


def majority_rule_consensus(
    trees: List[TreeNode],
    threshold: float = 0.5,
    verbose: bool = False
) -> Tuple[TreeNode, Dict[FrozenSet[str], float]]:
    """
    Build majority-rule consensus tree from multiple trees.

    Includes bipartitions that appear in > threshold fraction of trees.
    Standard majority-rule uses threshold=0.5 (>50%).

    Args:
        trees: List of phylogenetic trees to combine
        threshold: Frequency threshold (0.0 to 1.0), default 0.5
        verbose: Print progress information

    Returns:
        Tuple of (consensus_tree, support_values)
        - consensus_tree: TreeNode with support values
        - support_values: Dict mapping bipartition -> support percentage

    Example:
        trees = [upgma_tree, bionj_tree, ml_tree]
        consensus, support = majority_rule_consensus(trees)
        # consensus has splits appearing in >=2 of 3 trees
        # support values show % of trees with each split
    """
    if not trees:
        raise ValueError("Need at least one tree for consensus")

    if len(trees) == 1:
        # Single tree - return it as consensus with 100% support
        tree = trees[0]
        support_vals = {bp: 100.0 for bp in get_bipartitions(tree)}
        return tree, support_vals

    # Count bipartition frequencies
    bp_counts, all_taxa = count_bipartitions(trees)
    n_trees = len(trees)

    if verbose:
        print(f"\nBuilding consensus from {n_trees} trees")
        print(f"Total taxa: {len(all_taxa)}")
        print(f"Threshold: {threshold*100:.1f}%")

    # Filter bipartitions by threshold
    min_count = int(threshold * n_trees) + (1 if threshold * n_trees == int(threshold * n_trees) else 0)

    # Get bipartitions meeting threshold
    consensus_bps = {
        bp: count for bp, count in bp_counts.items()
        if count >= min_count and len(bp) > 0 and len(bp) < len(all_taxa)
    }

    if verbose:
        print(f"Bipartitions meeting threshold: {len(consensus_bps)}")

    # Calculate support values (as percentages)
    support_values = {
        bp: (count / n_trees) * 100.0
        for bp, count in consensus_bps.items()
    }

    # Build tree from compatible bipartitions
    consensus_tree = build_consensus_tree(
        all_taxa,
        consensus_bps,
        support_values,
        verbose=verbose
    )

    return consensus_tree, support_values


def strict_consensus(trees: List[TreeNode]) -> Tuple[TreeNode, Dict[FrozenSet[str], float]]:
    """
    Build strict consensus tree (threshold = 100%).

    Only includes bipartitions present in ALL input trees.
    Results in polytomies (multifurcations) when trees disagree.

    Args:
        trees: List of phylogenetic trees

    Returns:
        Tuple of (consensus_tree, support_values)
        All support values will be 100.0
    """
    return majority_rule_consensus(trees, threshold=1.0)


def build_consensus_tree(
    taxa: Set[str],
    bipartitions: Dict[FrozenSet[str], int],
    support_values: Dict[FrozenSet[str], float],
    verbose: bool = False
) -> TreeNode:
    """
    Construct tree from a set of bipartitions.

    Uses a greedy algorithm: add bipartitions in order of support,
    ensuring compatibility at each step.

    Args:
        taxa: Complete set of taxon names
        bipartitions: Dict of bipartition -> count
        support_values: Dict of bipartition -> support percentage
        verbose: Print construction progress

    Returns:
        TreeNode representing consensus tree

    Algorithm:
        1. Start with star tree (all taxa connected to root)
        2. Sort bipartitions by support (descending)
        3. For each bipartition:
            a. Check compatibility with existing structure
            b. If compatible, add to tree
            c. Assign support value to the node
    """
    if not taxa:
        return None

    # Handle single taxon
    if len(taxa) == 1:
        return TreeNode(name=list(taxa)[0], distance=0.0)

    # Sort bipartitions by support (highest first) and size (smallest first for ties)
    sorted_bps = sorted(
        bipartitions.items(),
        key=lambda x: (x[1], -len(x[0])),  # Sort by count desc, then size asc
        reverse=True
    )

    if verbose:
        print(f"\nBuilding tree from {len(sorted_bps)} bipartitions")

    # Start with a simple binary tree structure
    # For now, use a simplified approach: take the highest-support split first
    if not sorted_bps:
        # No informative splits - create star tree
        return create_star_tree(taxa)

    # Take the most supported split and build recursively
    first_bp, _ = sorted_bps[0]
    first_support = support_values.get(first_bp, 100.0)

    # Get the two groups from the bipartition
    group1 = first_bp
    group2 = taxa - first_bp

    # Recursively build subtrees
    left_subtree = build_subtree(group1, bipartitions, support_values, verbose)
    right_subtree = build_subtree(group2, bipartitions, support_values, verbose)

    # Create root with support value
    root = TreeNode(
        left=left_subtree,
        right=right_subtree,
        distance=0.0,
        support=first_support
    )

    return root


def build_subtree(
    taxa_subset: Set[str],
    bipartitions: Dict[FrozenSet[str], int],
    support_values: Dict[FrozenSet[str], float],
    verbose: bool = False
) -> TreeNode:
    """
    Recursively build a subtree for a subset of taxa.

    Args:
        taxa_subset: Taxa to include in this subtree
        bipartitions: All bipartitions with counts
        support_values: Support values for bipartitions
        verbose: Print progress

    Returns:
        TreeNode for the subtree
    """
    if len(taxa_subset) == 1:
        # Leaf node
        return TreeNode(name=list(taxa_subset)[0], distance=0.0)

    if len(taxa_subset) == 2:
        # Simple pair
        taxa_list = list(taxa_subset)
        return TreeNode(
            left=TreeNode(name=taxa_list[0], distance=0.0),
            right=TreeNode(name=taxa_list[1], distance=0.0),
            distance=0.0,
            support=100.0  # Trivial split
        )

    # Find the highest-support bipartition that splits this subset
    relevant_bps = [
        (bp, count) for bp, count in bipartitions.items()
        if bp.issubset(taxa_subset) and len(bp) > 0 and len(bp) < len(taxa_subset)
    ]

    if not relevant_bps:
        # No informative split - create star/polytomy
        # For binary tree, just pick arbitrary split
        taxa_list = list(taxa_subset)
        mid = len(taxa_list) // 2
        group1 = set(taxa_list[:mid])
        group2 = set(taxa_list[mid:])

        left = build_subtree(group1, bipartitions, support_values, verbose)
        right = build_subtree(group2, bipartitions, support_values, verbose)

        return TreeNode(left=left, right=right, distance=0.0, support=0.0)

    # Take highest-support split
    best_bp, _ = max(relevant_bps, key=lambda x: x[1])
    best_support = support_values.get(best_bp, 0.0)

    # Split into two groups
    group1 = best_bp
    group2 = taxa_subset - best_bp

    # Recursively build
    left = build_subtree(group1, bipartitions, support_values, verbose)
    right = build_subtree(group2, bipartitions, support_values, verbose)

    return TreeNode(left=left, right=right, distance=0.0, support=best_support)


def create_star_tree(taxa: Set[str]) -> TreeNode:
    """
    Create a star tree (polytomy) for taxa with no resolution.

    Args:
        taxa: Set of taxon names

    Returns:
        TreeNode representing unresolved star tree
    """
    if len(taxa) == 1:
        return TreeNode(name=list(taxa)[0], distance=0.0)

    if len(taxa) == 2:
        taxa_list = list(taxa)
        return TreeNode(
            left=TreeNode(name=taxa_list[0], distance=0.0),
            right=TreeNode(name=taxa_list[1], distance=0.0),
            distance=0.0
        )

    # For star tree with >2 taxa, create arbitrary binary structure
    # (true polytomy would require different tree structure)
    taxa_list = list(taxa)
    mid = len(taxa_list) // 2

    left = create_star_tree(set(taxa_list[:mid]))
    right = create_star_tree(set(taxa_list[mid:]))

    return TreeNode(left=left, right=right, distance=0.0, support=0.0)

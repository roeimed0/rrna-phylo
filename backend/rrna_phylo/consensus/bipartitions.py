"""
Bipartition extraction and manipulation for phylogenetic trees.

A bipartition (or split) is a division of taxa into two groups,
created by removing an internal edge from the tree.
"""

from typing import Set, List, FrozenSet
from rrna_phylo.core.tree import TreeNode


def get_leaf_names(node: TreeNode) -> Set[str]:
    """
    Get all leaf names in the subtree rooted at this node.

    Args:
        node: Root of subtree

    Returns:
        Set of leaf names
    """
    if node.is_leaf():
        return {node.name}

    names = set()
    if node.left:
        names.update(get_leaf_names(node.left))
    if node.right:
        names.update(get_leaf_names(node.right))

    return names


def get_bipartitions(tree: TreeNode) -> Set[FrozenSet[str]]:
    """
    Extract all bipartitions from a phylogenetic tree.

    A bipartition divides taxa into two groups. For each internal edge,
    we get the taxa on one side vs. the other side.

    We use canonical form: store the smaller of the two groups.
    This ensures {A,B}|{C,D} and {C,D}|{A,B} are treated as the same.

    Args:
        tree: Root of the tree

    Returns:
        Set of bipartitions, each as frozenset of taxa names

    Example:
        Tree: ((A,B),(C,D))
        Returns: {frozenset({'A','B'}), frozenset({'A'}), frozenset({'B'}),
                  frozenset({'C'}), frozenset({'D'})}
    """
    bipartitions = set()
    all_taxa = get_leaf_names(tree)

    def extract_splits(node: TreeNode):
        """Recursively extract bipartitions from tree."""
        if node.is_leaf():
            # Trivial split: single taxon vs rest
            bipartitions.add(frozenset({node.name}))
            return

        # Get taxa on left and right sides
        if node.left:
            left_taxa = get_leaf_names(node.left)
            # Store smaller group (canonical form)
            if len(left_taxa) <= len(all_taxa) / 2:
                bipartitions.add(frozenset(left_taxa))
            else:
                complement = all_taxa - left_taxa
                bipartitions.add(frozenset(complement))

            extract_splits(node.left)

        if node.right:
            right_taxa = get_leaf_names(node.right)
            # Store smaller group (canonical form)
            if len(right_taxa) <= len(all_taxa) / 2:
                bipartitions.add(frozenset(right_taxa))
            else:
                complement = all_taxa - right_taxa
                bipartitions.add(frozenset(complement))

            extract_splits(node.right)

    extract_splits(tree)
    return bipartitions


def get_informative_bipartitions(tree: TreeNode) -> Set[FrozenSet[str]]:
    """
    Get only informative (non-trivial) bipartitions.

    Informative bipartitions have at least 2 taxa on each side.
    Trivial bipartitions (single taxon vs rest) are excluded.

    Args:
        tree: Root of the tree

    Returns:
        Set of informative bipartitions

    Example:
        Tree: ((A,B),(C,D))
        Informative: {frozenset({'A','B'})}  # C,D implied
        Trivial (excluded): {A}, {B}, {C}, {D}
    """
    all_bipartitions = get_bipartitions(tree)
    all_taxa = get_leaf_names(tree)

    # Filter out trivial splits (size 1 or n-1)
    informative = {
        bp for bp in all_bipartitions
        if len(bp) > 1 and len(bp) < len(all_taxa) - 1
    }

    return informative


def are_compatible(bp1: FrozenSet[str], bp2: FrozenSet[str], all_taxa: Set[str]) -> bool:
    """
    Check if two bipartitions are compatible (can coexist in same tree).

    Two bipartitions are compatible if one of these is true:
    1. They are disjoint (no overlap)
    2. One is a subset of the other
    3. One is a subset of the complement of the other

    Args:
        bp1: First bipartition
        bp2: Second bipartition
        all_taxa: Complete set of all taxa

    Returns:
        True if compatible, False if conflicting

    Example:
        all_taxa = {A, B, C, D}
        bp1 = {A, B}  # implies {C, D} on other side
        bp2 = {A, C}  # implies {B, D} on other side
        -> NOT compatible (conflict about A and B being together)

        bp1 = {A, B}  # implies {C, D}
        bp2 = {A}     # implies {B, C, D}
        -> Compatible (A is subset of {A,B})
    """
    # Get complements
    comp1 = all_taxa - bp1
    comp2 = all_taxa - bp2

    # Check all compatibility conditions
    # 1. bp1 and bp2 are disjoint
    if bp1.isdisjoint(bp2):
        return True

    # 2. One is subset of the other
    if bp1.issubset(bp2) or bp2.issubset(bp1):
        return True

    # 3. One is subset of the other's complement
    if bp1.issubset(comp2) or bp2.issubset(comp1):
        return True

    # 4. Complements follow same rules
    if comp1.isdisjoint(comp2):
        return True
    if comp1.issubset(bp2) or comp2.issubset(bp1):
        return True

    return False


def get_taxa_from_trees(trees: List[TreeNode]) -> Set[str]:
    """
    Get the complete set of taxa from multiple trees.

    Args:
        trees: List of trees

    Returns:
        Set of all unique taxon names

    Raises:
        ValueError: If trees have different taxa sets
    """
    if not trees:
        return set()

    # Get taxa from first tree
    all_taxa = get_leaf_names(trees[0])

    # Verify all trees have same taxa
    for tree in trees[1:]:
        taxa = get_leaf_names(tree)
        if taxa != all_taxa:
            raise ValueError(
                f"Trees have different taxa sets. "
                f"First tree: {sorted(all_taxa)}, "
                f"Other tree: {sorted(taxa)}"
            )

    return all_taxa

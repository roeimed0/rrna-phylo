"""
Bipartition extraction and manipulation for phylogenetic trees.

A bipartition (or split) is a division of taxa into two groups,
created by removing an internal edge from the tree.
"""

from typing import Set, List, FrozenSet
from rrna_phylo.core.tree import TreeNode


def get_leaf_names(node: TreeNode) -> Set[str]:
    """Get all leaf names in the subtree rooted at this node."""
    if node.is_leaf():
        return {node.name}

    names = set()
    if node.left:
        names.update(get_leaf_names(node.left))
    if node.right:
        names.update(get_leaf_names(node.right))

    return names


def get_bipartitions(tree: TreeNode) -> Set[FrozenSet[str]]:
    """Extract all bipartitions from a phylogenetic tree."""
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
    """Get only informative (non-trivial) bipartitions."""
    all_bipartitions = get_bipartitions(tree)
    all_taxa = get_leaf_names(tree)

    # Filter out trivial splits (size 1 or n-1)
    informative = {
        bp for bp in all_bipartitions
        if len(bp) > 1 and len(bp) < len(all_taxa) - 1
    }

    return informative


def are_compatible(bp1: FrozenSet[str], bp2: FrozenSet[str], all_taxa: Set[str]) -> bool:
    """Check if two bipartitions are compatible (can coexist in same tree)."""
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
    """Get the complete set of taxa from multiple trees."""
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

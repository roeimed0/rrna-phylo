"""
UPGMA (Unweighted Pair Group Method with Arithmetic Mean) tree building.

UPGMA is a simple hierarchical clustering method that assumes a constant
molecular clock (equal evolutionary rates across all lineages).

Note: TreeNode class is now imported from rrna_phylo.core.tree
"""

import numpy as np
from typing import List, Tuple, Dict
from rrna_phylo.core.tree import TreeNode


class UPGMABuilder:
    """Build phylogenetic tree using UPGMA algorithm."""

    def __init__(self):
        """Initialize UPGMA builder."""
        self.nodes = []
        self.cluster_sizes = {}

    def build_tree(self, distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
        """
        Build UPGMA tree from distance matrix.

        Algorithm:
        1. Start with each sequence as a cluster
        2. Repeatedly merge closest clusters
        3. Update distances using arithmetic mean
        4. Continue until one cluster remains

        Args:
            distance_matrix: NxN symmetric distance matrix
            labels: Names of sequences (in matrix order)

        Returns:
            Root node of the tree
        """
        n = len(labels)

        if n < 2:
            raise ValueError("Need at least 2 sequences to build tree")

        # Initialize: each sequence is a cluster
        clusters = {}
        for i, label in enumerate(labels):
            node = TreeNode(name=label)
            clusters[i] = node
            self.cluster_sizes[i] = 1

        # Copy distance matrix (we'll modify it)
        dist = distance_matrix.copy()

        # Active clusters
        active = set(range(n))

        # Keep merging until one cluster remains
        while len(active) > 1:
            # Find minimum distance
            min_dist = float('inf')
            min_i, min_j = -1, -1

            for i in active:
                for j in active:
                    if i < j and dist[i][j] < min_dist:
                        min_dist = dist[i][j]
                        min_i, min_j = i, j

            # Merge clusters i and j
            cluster_i = clusters[min_i]
            cluster_j = clusters[min_j]

            # Calculate height of new cluster
            new_height = min_dist / 2.0

            # Create new internal node
            new_node = TreeNode(
                name=f"({cluster_i.name},{cluster_j.name})",
                left=cluster_i,
                right=cluster_j
            )
            new_node.height = new_height

            # Set branch lengths
            cluster_i.distance = new_height - cluster_i.height
            cluster_j.distance = new_height - cluster_j.height

            # Update distances to new cluster (UPGMA: arithmetic mean)
            new_cluster_idx = min_i  # Reuse index
            size_i = self.cluster_sizes[min_i]
            size_j = self.cluster_sizes[min_j]
            new_size = size_i + size_j

            for k in active:
                if k != min_i and k != min_j:
                    # Average distance
                    new_dist = (dist[min_i][k] * size_i + dist[min_j][k] * size_j) / new_size
                    dist[new_cluster_idx][k] = new_dist
                    dist[k][new_cluster_idx] = new_dist

            # Update cluster
            clusters[new_cluster_idx] = new_node
            self.cluster_sizes[new_cluster_idx] = new_size

            # Remove old cluster j
            active.remove(min_j)

        # Return root (last remaining cluster)
        root_idx = list(active)[0]
        root = clusters[root_idx]
        root.distance = 0.0  # Root has no parent

        return root

    def print_tree(self, node: TreeNode, indent: int = 0):
        """
        Print tree in ASCII format.

        Args:
            node: Tree node to print
            indent: Indentation level
        """
        prefix = "  " * indent

        if node.is_leaf():
            print(f"{prefix}└─ {node.name} ({node.distance:.4f})")
        else:
            print(f"{prefix}└─ Internal ({node.distance:.4f})")
            if node.left:
                self.print_tree(node.left, indent + 1)
            if node.right:
                self.print_tree(node.right, indent + 1)


def build_upgma_tree(distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
    """
    Convenience function to build UPGMA tree.

    Args:
        distance_matrix: NxN symmetric distance matrix
        labels: Sequence names

    Returns:
        Root node of UPGMA tree
    """
    builder = UPGMABuilder()
    return builder.build_tree(distance_matrix, labels)

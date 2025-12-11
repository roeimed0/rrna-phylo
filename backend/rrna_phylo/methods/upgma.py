"""UPGMA tree building using hierarchical clustering."""

import numpy as np
from typing import List, Tuple, Dict
from rrna_phylo.core.tree import TreeNode

class UPGMABuilder:
    """Build phylogenetic tree using UPGMA algorithm."""

    def __init__(self):
        self.cluster_sizes = {}

    def build_tree(self, distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
        n = len(labels)
        if n < 2:
            raise ValueError("Need at least 2 sequences")

        clusters = {}
        for i, label in enumerate(labels):
            node = TreeNode(name=label)
            node.height = 0.0                # <-- IMPORTANT FIX
            clusters[i] = node
            self.cluster_sizes[i] = 1

        dist = distance_matrix.copy()
        active = set(range(n))

        while len(active) > 1:
            min_dist = float('inf')
            min_i = min_j = -1

            for i in active:
                for j in active:
                    if i < j and dist[i][j] < min_dist:
                        min_dist = dist[i][j]
                        min_i, min_j = i, j

            ci = clusters[min_i]
            cj = clusters[min_j]

            new_height = min_dist / 2.0

            new_node = TreeNode(
                name=None,          # cleaner, not fake Newick
                left=ci,
                right=cj
            )
            new_node.height = new_height

            ci.distance = new_height - ci.height
            cj.distance = new_height - cj.height

            size_i = self.cluster_sizes[min_i]
            size_j = self.cluster_sizes[min_j]
            new_size = size_i + size_j

            for k in active:
                if k != min_i and k != min_j:
                    new_dist = (dist[min_i][k] * size_i + dist[min_j][k] * size_j) / new_size
                    dist[min_i][k] = dist[k][min_i] = new_dist

            clusters[min_i] = new_node
            self.cluster_sizes[min_i] = new_size

            active.remove(min_j)

        root = clusters[list(active)[0]]
        root.distance = 0.0
        return root

def build_upgma_tree(distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
    """Build UPGMA tree from distance matrix and sequence labels."""
    builder = UPGMABuilder()
    return builder.build_tree(distance_matrix, labels)

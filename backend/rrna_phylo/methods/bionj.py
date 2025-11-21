"""
BioNJ - Improved Neighbor-Joining with variance weighting.

BioNJ (Gascuel, 1997) is an enhancement of the classic Neighbor-Joining algorithm.
It uses variance-based weighting to produce more accurate trees, especially when
evolutionary rates vary significantly across lineages.

Reference: Gascuel O. (1997) BIONJ: an improved version of the NJ algorithm
based on a simple model of sequence data. Mol Biol Evol. 14(7):685-95.
"""

import numpy as np
from typing import List, Tuple
from rrna_phylo.core.tree import TreeNode


class BioNJBuilder:
    """Build phylogenetic tree using BioNJ algorithm."""

    def __init__(self):
        """Initialize BioNJ builder."""
        self.nodes = {}
        self.variance = None

    def build_tree(self, distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
        """
        Build BioNJ tree from distance matrix.

        BioNJ Algorithm:
        1. Calculate variance matrix (estimates reliability of distances)
        2. Find pair (i,j) that minimizes Q criterion (like NJ)
        3. Create new node joining i and j
        4. Calculate branch lengths using ML estimates
        5. Update distances using VARIANCE-WEIGHTED formula (BioNJ innovation)
        6. Repeat until one node remains

        Args:
            distance_matrix: NxN symmetric distance matrix
            labels: Names of sequences (in matrix order)

        Returns:
            Root node of the tree
        """
        n = len(labels)

        if n < 2:
            raise ValueError("Need at least 2 sequences to build tree")

        # Initialize nodes (each sequence is a leaf)
        for i, label in enumerate(labels):
            self.nodes[i] = TreeNode(name=label)

        # Copy matrices (we'll modify them)
        dist = distance_matrix.copy()

        # Initialize variance matrix (simple model: variance proportional to distance)
        self.variance = distance_matrix.copy()

        # Active nodes
        active = set(range(n))

        # Keep merging until one node remains
        while len(active) > 1:
            # Calculate Q matrix (Neighbor-Joining criterion)
            Q = self._calculate_Q_matrix(dist, active)

            # Find minimum Q value
            min_q = float('inf')
            min_i, min_j = -1, -1

            for i in active:
                for j in active:
                    if i < j and Q[i][j] < min_q:
                        min_q = Q[i][j]
                        min_i, min_j = i, j

            # Calculate branch lengths
            branch_i, branch_j = self._calculate_branch_lengths(
                dist, min_i, min_j, active
            )

            # Create new internal node
            node_i = self.nodes[min_i]
            node_j = self.nodes[min_j]

            new_node = TreeNode(
                name=f"({node_i.name},{node_j.name})",
                left=node_i,
                right=node_j
            )

            # Set branch lengths
            node_i.distance = max(0.0, branch_i)  # Ensure non-negative
            node_j.distance = max(0.0, branch_j)

            # Update distance matrix using BioNJ variance weighting
            new_idx = min_i  # Reuse index
            for k in active:
                if k != min_i and k != min_j:
                    # BioNJ innovation: variance-weighted distance update
                    new_dist = self._bionj_distance_update(
                        dist, self.variance, min_i, min_j, k
                    )
                    dist[new_idx][k] = new_dist
                    dist[k][new_idx] = new_dist

                    # Update variance
                    self.variance[new_idx][k] = new_dist
                    self.variance[k][new_idx] = new_dist

            # Update nodes
            self.nodes[new_idx] = new_node

            # Remove old node j
            active.remove(min_j)

        # Return root (last remaining node)
        root_idx = list(active)[0]
        root = self.nodes[root_idx]
        root.distance = 0.0  # Root has no parent

        return root

    def _calculate_Q_matrix(self, dist: np.ndarray, active: set) -> np.ndarray:
        """
        Calculate Q matrix for Neighbor-Joining criterion.

        Q[i,j] = (n-2) * d[i,j] - sum(d[i,k]) - sum(d[j,k])

        This finds the pair that minimizes total tree length.

        Args:
            dist: Distance matrix
            active: Set of active node indices

        Returns:
            Q matrix
        """
        n_active = len(active)
        active_list = sorted(list(active))
        max_idx = max(active_list) + 1

        Q = np.zeros((max_idx, max_idx))

        # Calculate row sums (needed for Q)
        row_sums = {}
        for i in active:
            row_sums[i] = sum(dist[i][k] for k in active)

        # Calculate Q matrix
        for i in active:
            for j in active:
                if i < j:
                    Q[i][j] = (n_active - 2) * dist[i][j] - row_sums[i] - row_sums[j]
                    Q[j][i] = Q[i][j]

        return Q

    def _calculate_branch_lengths(
        self,
        dist: np.ndarray,
        i: int,
        j: int,
        active: set
    ) -> Tuple[float, float]:
        """
        Calculate branch lengths for nodes i and j.

        Uses Neighbor-Joining formula with correction for unequal rates.

        Args:
            dist: Distance matrix
            i: First node index
            j: Second node index
            active: Set of active nodes

        Returns:
            Tuple of (branch_length_i, branch_length_j)
        """
        n_active = len(active)

        # Calculate sum of distances
        sum_i = sum(dist[i][k] for k in active if k != i and k != j)
        sum_j = sum(dist[j][k] for k in active if k != i and k != j)

        # NJ branch length formula
        if n_active == 2:
            # Special case: only two nodes left
            branch_i = dist[i][j] / 2.0
            branch_j = dist[i][j] / 2.0
        else:
            branch_i = dist[i][j] / 2.0 + (sum_i - sum_j) / (2 * (n_active - 2))
            branch_j = dist[i][j] - branch_i

        return branch_i, branch_j

    def _bionj_distance_update(
        self,
        dist: np.ndarray,
        var: np.ndarray,
        i: int,
        j: int,
        k: int
    ) -> float:
        """
        Update distance using BioNJ variance weighting.

        This is the KEY INNOVATION of BioNJ over standard NJ:
        Uses variance to weight which distance (i or j) to trust more.

        Formula:
        lambda = 0.5 + (var[i,k] - var[j,k]) / (2 * var[i,j])
        d[ij,k] = lambda * d[i,k] + (1-lambda) * d[j,k] - lambda * d[i,j]/2

        If var[i,k] > var[j,k], then d[j,k] is more reliable, so lambda < 0.5

        Args:
            dist: Distance matrix
            var: Variance matrix
            i: First merged node
            j: Second merged node
            k: Node to update distance to

        Returns:
            Updated distance from new cluster (i,j) to k
        """
        d_ij = dist[i][j]
        d_ik = dist[i][k]
        d_jk = dist[j][k]
        v_ij = var[i][j]
        v_ik = var[i][k]
        v_jk = var[j][k]

        # Avoid division by zero
        if abs(v_ij) < 1e-10:
            # Fall back to simple average (standard NJ)
            lambda_val = 0.5
        else:
            # BioNJ variance weighting
            lambda_val = 0.5 + (v_ik - v_jk) / (2 * v_ij)

            # Ensure lambda is in reasonable range [0, 1]
            lambda_val = max(0.0, min(1.0, lambda_val))

        # Calculate new distance
        new_dist = lambda_val * d_ik + (1 - lambda_val) * d_jk - lambda_val * d_ij

        return max(0.0, new_dist)  # Ensure non-negative

    def print_tree(self, node: TreeNode, indent: int = 0):
        """
        Print tree in ASCII format.

        Args:
            node: Tree node to print
            indent: Indentation level
        """
        prefix = "  " * indent

        if node.is_leaf():
            print(f"{prefix}`- {node.name} ({node.distance:.4f})")
        else:
            print(f"{prefix}`- Internal ({node.distance:.4f})")
            if node.left:
                self.print_tree(node.left, indent + 1)
            if node.right:
                self.print_tree(node.right, indent + 1)


def build_bionj_tree(distance_matrix: np.ndarray, labels: List[str]) -> TreeNode:
    """
    Convenience function to build BioNJ tree.

    Args:
        distance_matrix: NxN symmetric distance matrix
        labels: Sequence names

    Returns:
        Root node of BioNJ tree
    """
    builder = BioNJBuilder()
    return builder.build_tree(distance_matrix, labels)

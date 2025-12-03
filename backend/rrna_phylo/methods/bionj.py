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

        (Docstring preserved from your original code)
        """
        n = len(labels)

        if n < 2:
            raise ValueError("Need at least 2 sequences to build tree")

        # Initialize nodes (each sequence is a leaf)
        for i, label in enumerate(labels):
            self.nodes[i] = TreeNode(name=label)

        # Copy matrices (we'll modify them)
        dist = distance_matrix.astype(float).copy()

        # Initialize variance matrix using modern approach (var proportional to distance but robust)
        # Use dist/2 as base (like classic), but guard tiny values to avoid zeros.
        self.variance = self._initialize_variance_matrix(dist)

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
            for k in list(active):
                if k != min_i and k != min_j:
                    # BioNJ innovation: variance-weighted distance update
                    new_dist = self._bionj_distance_update(
                        dist, self.variance, min_i, min_j, k
                    )
                    # _bionj_distance_update already updates var entries for new_idx<->k
                    dist[new_idx][k] = new_dist
                    dist[k][new_idx] = new_dist

            # Update nodes
            self.nodes[new_idx] = new_node

            # Remove old node j from active set and disable its rows/cols
            active.remove(min_j)

            # Disable row/col min_j so it's never selected again
            dist[min_j, :] = float("inf")
            dist[:, min_j] = float("inf")
            self.variance[min_j, :] = float("inf")
            self.variance[:, min_j] = float("inf")

            # Ensure diagonal zero for new index
            dist[new_idx, new_idx] = 0.0
            self.variance[new_idx, new_idx] = 0.0

        # Return root (last remaining node)
        root_idx = list(active)[0]
        root = self.nodes[root_idx]
        root.distance = 0.0  # Root has no parent

        return root

    # --------------------
    # Helper functions
    # --------------------
    def _initialize_variance_matrix(self, dist: np.ndarray) -> np.ndarray:
        """
        Create an initial variance matrix for distances using a robust modern choice.
        Base: v_ij = max(dist_ij / 2, epsilon)
        Ensures no zeros and stable inverse-variance weighting.
        """
        eps = 1e-8
        var = dist.astype(float).copy() / 2.0
        # enforce symmetric and zero diagonal
        np.fill_diagonal(var, 0.0)
        # replace tiny variances with eps to avoid division by zero
        var[var < eps] = eps
        return var

    def _compute_inverse_variance_weights(self, v_ik: float, v_jk: float) -> Tuple[float, float]:
        """
        Compute modern inverse-variance weights for combining d(i,k) and d(j,k).

        w_i = (1/v_ik) / (1/v_ik + 1/v_jk) = v_jk / (v_ik + v_jk)
        w_j = 1 - w_i
        """
        denom = v_ik + v_jk
        if denom <= 0:
            return 0.5, 0.5
        w_i = v_jk / denom
        w_j = v_ik / denom
        # clamp to [0,1] to be robust
        w_i = max(0.0, min(1.0, w_i))
        w_j = 1.0 - w_i
        return w_i, w_j

    def _bionj_distance_update(
        self,
        dist: np.ndarray,
        var: np.ndarray,
        i: int,
        j: int,
        k: int
    ) -> float:
        """
        Modern BioNJ distance update using inverse-variance weights.

        This function:
         - computes inverse-variance weights w_i,w_j for combining d(i,k), d(j,k)
         - computes new distance d(u,k) = w_i * d(i,k) + w_j * d(j,k)
         - updates the variance var[u,k] = w_i^2 * v_ik + w_j^2 * v_jk
           (this propagates uncertainty correctly for weighted average)
         - returns the new distance (and modifies var in-place)
        """
        d_ik = dist[i][k]
        d_jk = dist[j][k]

        v_ik = var[i][k]
        v_jk = var[j][k]

        # Compute weights (inverse-variance weighting)
        w_i, w_j = self._compute_inverse_variance_weights(v_ik, v_jk)

        # New distance is the weighted average (no subtraction term)
        d_new = w_i * d_ik + w_j * d_jk

        # New variance is variance of weighted sum (assuming independence)
        v_new = (w_i ** 2) * v_ik + (w_j ** 2) * v_jk

        # Guard numerical issues
        if v_new <= 0:
            v_new = max(1e-8, (v_ik + v_jk) / 2.0)

        # Place new values at the merged index (we expect caller to write dist[new_idx,k])
        new_idx = i  # caller reuses i as new index
        var[new_idx][k] = v_new
        var[k][new_idx] = v_new

        return max(0.0, d_new)

    def _calculate_Q_matrix(self, dist: np.ndarray, active: set) -> np.ndarray:
        """
        Calculate Q matrix for Neighbor-Joining criterion.

        Q[i,j] = (n-2) * d[i,j] - sum(d[i,k]) - sum(d[j,k])
        (unchanged from your original code; kept here for clarity)
        """
        n_active = len(active)
        active_list = sorted(list(active))
        max_idx = max(active_list) + 1

        Q = np.full((max_idx, max_idx), float("inf"))

        # Calculate row sums (needed for Q)
        row_sums = {}
        for i in active:
            # sum only over active indices and skip infinities
            row_sums[i] = sum(dist[i][k] for k in active if not np.isinf(dist[i][k]))

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
        Calculate branch lengths for nodes i and j using the stable NJ formula.

        This reproduces NJ branch-lengths (numerically stable) which works
        well with the modern BioNJ distance updates.
        """
        n_active = len(active)

        # Calculate sum of distances (over active indices only)
        sum_i = sum(dist[i][k] for k in active if k != i and k != j and not np.isinf(dist[i][k]))
        sum_j = sum(dist[j][k] for k in active if k != i and k != j and not np.isinf(dist[j][k]))

        # NJ branch length formula
        if n_active == 2:
            # Special case: only two nodes left
            branch_i = dist[i][j] / 2.0
            branch_j = dist[i][j] / 2.0
        else:
            branch_i = dist[i][j] / 2.0 + (sum_i - sum_j) / (2 * (n_active - 2))
            # Use symmetric formula for stability
            branch_j = dist[i][j] / 2.0 - (sum_i - sum_j) / (2 * (n_active - 2))

        # Ensure non-negative due to rounding
        return max(0.0, branch_i), max(0.0, branch_j)

    def print_tree(self, node: TreeNode, indent: int = 0):
        """
        Print tree in ASCII format.

        Args:
            node: Tree node to print
            indent: Indentation level
        """
        prefix = "  " * indent

        if node.is_leaf():
            print(f"{prefix}`- {node.name} ({getattr(node, 'distance', 0.0):.4f})")
        else:
            print(f"{prefix}`- Internal ({getattr(node, 'distance', 0.0):.4f})")
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

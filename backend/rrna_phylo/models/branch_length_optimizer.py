"""
Branch length optimization for phylogenetic trees.

After topology changes (like NNI swaps), branch lengths must be re-optimized
to maximize likelihood under the new topology.

Key Insight:
===========
When NNI swaps the tree topology, the old branch lengths were optimized for
the OLD topology. Without re-optimization, branches gradually collapse to zero
as the topology changes, creating invalid star-like trees.

Solution:
=========
After each accepted NNI swap, re-optimize all branch lengths using coordinate
descent (optimize each branch individually, iterate until convergence).
"""

from typing import List, Optional
import numpy as np
from scipy.optimize import minimize_scalar
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood


def optimize_single_branch(
    tree: TreeNode,
    node: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001,
    max_length: float = 10.0,
    calculator = None
) -> float:
    """
    Optimize the branch length leading to a single node.

    Uses Brent's method (scipy.optimize.minimize_scalar) to find the branch
    length that maximizes likelihood.

    Args:
        tree: Full tree (root node)
        node: Node whose incoming branch length to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter for rate heterogeneity
        min_length: Minimum branch length (prevents zero collapse)
        max_length: Maximum branch length (prevents infinite branches)
        calculator: Optional cached LikelihoodCalculatorLevel3 (HUGE speedup!)

    Returns:
        Optimized branch length
    """
    # Save original length (in case optimization fails)
    original_length = node.distance if hasattr(node, 'distance') else 0.01

    def neg_log_likelihood(length):
        """Negative log-likelihood for minimization."""
        # Set branch length (constrained to valid range)
        node.distance = max(min_length, min(length, max_length))

        # Compute likelihood (with cached calculator for speed)
        logL = compute_log_likelihood(tree, sequences, alpha=alpha, calculator=calculator)

        # Return negative (scipy minimizes, we want to maximize logL)
        return -logL

    try:
        # Optimize using Brent's method (efficient for 1D optimization)
        result = minimize_scalar(
            neg_log_likelihood,
            bounds=(min_length, max_length),
            method='bounded',
            options={'xatol': 0.0001}  # Tolerance
        )

        # Set optimized length
        optimal_length = max(min_length, min(result.x, max_length))
        node.distance = optimal_length

        return optimal_length

    except Exception as e:
        # If optimization fails, restore original length
        node.distance = original_length
        return original_length


def optimize_all_branch_lengths(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001,
    max_iterations: int = 3,
    verbose: bool = False
) -> float:
    """
    Optimize all branch lengths in a tree.

    Uses coordinate descent: optimize each branch in turn, repeat until convergence.

    This is the "full" optimization - iterates multiple times through all branches.
    Good for final tree optimization, but slower than single-pass.

    Args:
        tree: Tree with topology to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter
        min_length: Minimum branch length (prevents zero collapse)
        max_iterations: Number of passes through all branches
        verbose: Print progress

    Returns:
        Final log-likelihood after optimization
    """
    if verbose:
        print("  Optimizing all branch lengths (full)...")

    # Get all nodes with incoming branches (exclude root)
    all_nodes = []
    def collect_nodes(node):
        if hasattr(node, 'distance') and node.distance is not None:
            all_nodes.append(node)
        if not node.is_leaf():
            if node.left:
                collect_nodes(node.left)
            if node.right:
                collect_nodes(node.right)
    collect_nodes(tree)

    # Coordinate descent: optimize each branch in turn
    for iteration in range(max_iterations):
        for node in all_nodes:
            optimize_single_branch(tree, node, sequences, alpha, min_length)

        if verbose and iteration < max_iterations - 1:
            logL = compute_log_likelihood(tree, sequences, alpha=alpha)
            print(f"    Iteration {iteration + 1}/{max_iterations}: LogL = {logL:.2f}")

    # Final likelihood
    final_logL = compute_log_likelihood(tree, sequences, alpha=alpha)

    if verbose:
        print(f"    Final LogL after branch optimization: {final_logL:.2f}")

    return final_logL


def optimize_branch_lengths_fast(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001,
    verbose: bool = False
) -> float:
    """
    Fast single-pass branch length optimization.

    Optimizes each branch once (no iteration). Faster but less accurate than
    full optimization. Good enough for NNI search where topology is changing
    frequently.

    PERFORMANCE: Creates calculator once and reuses it for all branches!
    This avoids recreating the site pattern compressor hundreds of times.

    Args:
        tree: Tree to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter
        min_length: Minimum branch length
        verbose: Print progress

    Returns:
        Final log-likelihood
    """
    if verbose:
        print("  Re-optimizing branch lengths (fast)...")

    # Create calculator once and reuse (HUGE speedup!)
    from rrna_phylo.models.ml_tree import GTRModel
    from rrna_phylo.models.ml_tree_level3 import LikelihoodCalculatorLevel3

    model = GTRModel()
    model.estimate_parameters(sequences)

    if alpha is None:
        alpha = 1.0

    calculator = LikelihoodCalculatorLevel3(model, sequences, alpha=alpha)

    # Get all nodes with incoming branches
    all_nodes = []
    def collect_nodes(node):
        if hasattr(node, 'distance') and node.distance is not None:
            all_nodes.append(node)
        if not node.is_leaf():
            if node.left:
                collect_nodes(node.left)
            if node.right:
                collect_nodes(node.right)
    collect_nodes(tree)

    # Single pass optimization (with cached calculator)
    for node in all_nodes:
        optimize_single_branch(tree, node, sequences, alpha, min_length, calculator=calculator)

    # Return final likelihood (using cached calculator)
    final_logL = compute_log_likelihood(tree, sequences, alpha=alpha, calculator=calculator)

    if verbose:
        print(f"    LogL after fast optimization: {final_logL:.2f}")

    return final_logL


def check_for_zero_branches(tree: TreeNode, verbose: bool = False) -> int:
    """
    Count zero-length branches in a tree.

    Useful for debugging NNI issues. A properly optimized tree should have
    very few (ideally zero) branches with exactly 0.0 length.

    Args:
        tree: Tree to check
        verbose: Print details about zero branches

    Returns:
        Number of branches with exactly 0.0 length
    """
    zero_branches = []

    def check_node(node):
        if hasattr(node, 'distance'):
            if node.distance == 0.0:
                zero_branches.append(node.name if node.is_leaf() else "internal")

        if not node.is_leaf():
            if node.left:
                check_node(node.left)
            if node.right:
                check_node(node.right)

    check_node(tree)

    if verbose and zero_branches:
        print(f"  WARNING: Found {len(zero_branches)} zero-length branches!")
        if len(zero_branches) <= 10:
            print(f"  Nodes: {', '.join(zero_branches)}")

    return len(zero_branches)

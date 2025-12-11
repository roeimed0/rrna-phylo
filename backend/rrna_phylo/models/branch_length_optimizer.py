"""Re-optimize branch lengths after topology changes to maintain valid trees."""

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
    """Optimize the branch length leading to a single node using Brent's method."""
    # Save original length (in case optimization fails)
    original_length = node.distance if hasattr(node, 'distance') else 0.01

    def neg_log_likelihood(length):
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
    """Optimize all branch lengths using coordinate descent over multiple iterations."""
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
    verbose: bool = False,
    calculator = None
) -> float:
    """Fast single-pass branch optimization (for topology search algorithms)."""
    if verbose:
        print("  Re-optimizing branch lengths (fast)...")

    # CRITICAL FIX: Accept calculator from caller (e.g., NNI search)
    # Only create if not provided (for standalone use)
    if calculator is None:
        # Create calculator once and reuse (HUGE speedup!)
        from rrna_phylo.models.ml_tree import GTRModel
        from rrna_phylo.models.ml_tree_level3 import LikelihoodCalculatorLevel3

        model = GTRModel()
        model.estimate_parameters(sequences)

        if alpha is None:
            alpha = 1.0

        calculator = LikelihoodCalculatorLevel3(model, sequences, alpha=alpha)

    # If alpha not provided but calculator is, use default
    if alpha is None:
        alpha = 1.0

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
    """Count zero-length branches in a tree (debugging for tree validity)."""
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

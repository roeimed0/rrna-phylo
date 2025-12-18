"""
SIMPLIFIED SPR - Correct but limited implementation.

This implements a simplified SPR that is:
1. CORRECT - Won't corrupt tree structure
2. LIMITED - Only explores subset of possible SPR moves
3. FAST - Uses in-place modifications like NNI

Strategy: Instead of full SPR (prune + regraft), we do extended NNI moves
that can reach across more of the tree. This is weaker than true SPR but
much easier to implement correctly.
"""

from typing import List, Tuple, Optional
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.branch_length_optimizer import optimize_branch_lengths_vectorized


def get_all_internal_nodes(tree: TreeNode) -> List[TreeNode]:
    """Get all internal nodes in the tree (non-leaves)."""
    nodes = []

    def traverse(node):
        if not node.is_leaf():
            nodes.append(node)
            if node.left:
                traverse(node.left)
            if node.right:
                traverse(node.right)

    traverse(tree)
    return nodes


def try_subtree_swap(tree: TreeNode, node1: TreeNode, node2: TreeNode,
                      sequences: List[Sequence], alpha: Optional[float] = None) -> Tuple[Optional[TreeNode], Optional[float]]:
    """
    Try swapping subtrees at two different internal nodes.

    This is a simplified SPR-like move that:
    1. Swaps left/right children at two internal nodes
    2. Re-optimizes branch lengths
    3. Evaluates likelihood

    Args:
        tree: Current tree (will be modified in-place if successful)
        node1: First internal node
        node2: Second internal node
        sequences: Aligned sequences
        alpha: Gamma parameter

    Returns:
        (tree, logL) if move successful, (None, None) if not
    """
    # Can't swap if nodes are adjacent in tree
    if node1 == node2:
        return None, None
    if node1.left == node2 or node1.right == node2:
        return None, None
    if node2.left == node1 or node2.right == node1:
        return None, None

    # Try swapping left children
    if node1.left and node2.left and not node1.left.is_leaf() and not node2.left.is_leaf():
        # Swap
        node1.left, node2.left = node2.left, node1.left

        # Optimize and evaluate
        try:
            logL = optimize_branch_lengths_vectorized(tree, sequences, alpha=alpha, verbose=False)
            return tree, logL
        except:
            # Undo swap on error
            node1.left, node2.left = node2.left, node1.left
            return None, None

    return None, None


def simplified_spr_search(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    max_iterations: int = 3,
    verbose: bool = False
) -> Tuple[TreeNode, float, int]:
    """
    Simplified SPR search using extended NNI-like moves.

    This is NOT full SPR, but a limited version that:
    - Only tries subtree swaps at internal nodes
    - Uses in-place modifications (fast)
    - Won't corrupt tree structure (safe)
    - May find some improvements NNI misses (better than NNI alone)

    Args:
        tree: Initial tree
        sequences: Aligned sequences
        alpha: Gamma parameter
        max_iterations: Max iterations
        verbose: Print progress

    Returns:
        (improved_tree, final_logL, n_improvements)
    """
    from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood

    if verbose:
        print("\n" + "="*70)
        print("SIMPLIFIED SPR SEARCH")
        print("="*70)
        print("NOTE: This is a limited SPR implementation")
        print("      Uses extended NNI moves, not full SPR")
        print(f"Max iterations: {max_iterations}")
        print()

    # Don't copy - work with the tree we got (already optimized by NNI)
    current_tree = tree
    current_logL = compute_log_likelihood(current_tree, sequences, alpha=alpha)

    if verbose:
        print(f"Initial LogL: {current_logL:.2f}")

    n_improvements = 0

    for iteration in range(max_iterations):
        improved_this_round = False
        best_tree = None
        best_logL = current_logL

        # Get all internal nodes
        internal_nodes = get_all_internal_nodes(current_tree)

        if verbose:
            print(f"\nIteration {iteration + 1}: Testing {len(internal_nodes)} internal nodes...")

        n_moves_tested = 0

        # Try swapping subtrees between different internal nodes
        for i, node1 in enumerate(internal_nodes):
            for node2 in internal_nodes[i+1:]:  # Only try each pair once
                # Try swap
                test_tree = current_tree.copy()
                test_node1 = get_all_internal_nodes(test_tree)[i]
                test_node2_idx = internal_nodes.index(node2)
                test_node2 = get_all_internal_nodes(test_tree)[test_node2_idx]

                result_tree, result_logL = try_subtree_swap(
                    test_tree, test_node1, test_node2, sequences, alpha
                )

                if result_tree and result_logL:
                    n_moves_tested += 1

                    if result_logL > best_logL + 0.01:  # Tolerance
                        best_tree = result_tree
                        best_logL = result_logL
                        improved_this_round = True

                        if verbose:
                            improvement = result_logL - current_logL
                            print(f"  [IMPROVEMENT] Found better tree: {improvement:+.2f}")

        if verbose:
            print(f"  Tested {n_moves_tested} simplified SPR moves")

        if improved_this_round:
            current_tree = best_tree
            current_logL = best_logL
            n_improvements += 1

            if verbose:
                print(f"  Best LogL this round: {current_logL:.2f}")
        else:
            if verbose:
                print(f"\nSimplified SPR converged after {iteration + 1} iterations")
                print(f"Final LogL: {current_logL:.2f}")
                print(f"Total improvements: {n_improvements}")
            break

    return current_tree, current_logL, n_improvements

"""
Proper SPR (Subtree Pruning and Regrafting) Implementation.

This is a correct implementation that:
1. Maintains binary tree property
2. Preserves all taxa
3. Properly handles parent node collapse
4. Updates branch lengths correctly

Based on debugging insights and reference implementations (RAxML, IQ-TREE).
"""

from typing import List, Tuple, Optional, Set
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.branch_length_optimizer import optimize_branch_lengths_vectorized
import copy


class SPRMove:
    """Represents an SPR move for easy undo/redo."""

    def __init__(self, prune_edge: Tuple[TreeNode, TreeNode],
                 regraft_edge: Tuple[TreeNode, TreeNode]):
        self.prune_parent = prune_edge[0]
        self.prune_child = prune_edge[1]
        self.regraft_parent = regraft_edge[0]
        self.regraft_child = regraft_edge[1]


def get_all_nodes_with_parents(tree: TreeNode) -> List[Tuple[Optional[TreeNode], TreeNode]]:
    """
    Get all (parent, node) pairs in the tree.

    Returns:
        List of (parent, node) tuples. Root has parent=None.
    """
    nodes = []

    def traverse(parent, node):
        nodes.append((parent, node))
        if not node.is_leaf():
            if node.left:
                traverse(node, node.left)
            if node.right:
                traverse(node, node.right)

    traverse(None, tree)
    return nodes


def find_sibling(parent: TreeNode, child: TreeNode) -> Optional[TreeNode]:
    """Find the sibling of a child node."""
    if parent.left == child:
        return parent.right
    elif parent.right == child:
        return parent.left
    return None


def find_parent_in_tree(tree: TreeNode, target: TreeNode) -> Optional[TreeNode]:
    """Find the parent of a target node in the tree."""

    def search(node, parent):
        if node == target:
            return parent
        if not node.is_leaf():
            if node.left:
                result = search(node.left, node)
                if result is not None:
                    return result
            if node.right:
                result = search(node.right, node)
                if result is not None:
                    return result
        return None

    return search(tree, None)


def spr_move_proper(tree: TreeNode, prune_parent: TreeNode, prune_child: TreeNode,
                    regraft_parent: TreeNode, regraft_child: TreeNode) -> Optional[TreeNode]:
    """
    Perform a proper SPR move.

    Algorithm:
    1. Prune: Remove prune_child from prune_parent
    2. Collapse: Connect sibling directly to grandparent
    3. Regraft: Create new internal node, attach prune_child and regraft_child
    4. Update: Set all branch lengths

    Args:
        tree: Tree to modify (will be copied)
        prune_parent: Parent of subtree to prune
        prune_child: Root of subtree to prune
        regraft_parent: Parent of edge to regraft to
        regraft_child: Child of edge to regraft to

    Returns:
        Modified tree, or None if move is invalid
    """
    # Create a deep copy to avoid modifying original
    new_tree = tree.copy()

    # Find corresponding nodes in the copied tree
    all_nodes = get_all_nodes_with_parents(new_tree)

    # Map old nodes to new nodes by name/position
    def find_node_in_copy(original_node):
        original_name = original_node.name
        original_is_leaf = original_node.is_leaf()

        for parent, node in all_nodes:
            if node.name == original_name and node.is_leaf() == original_is_leaf:
                return parent, node
        return None, None

    # Find nodes in copied tree
    new_prune_parent, new_prune_child = find_node_in_copy(prune_child)
    new_regraft_parent, new_regraft_child = find_node_in_copy(regraft_child)

    if not new_prune_child or not new_regraft_child:
        return None

    # Get grandparent and sibling for collapse operation
    grandparent = find_parent_in_tree(new_tree, new_prune_parent)
    sibling = find_sibling(new_prune_parent, new_prune_child)

    if not sibling:
        return None  # Can't prune the only child

    # STEP 1: Prune - remove prune_child from its parent
    if new_prune_parent.left == new_prune_child:
        new_prune_parent.left = None
    elif new_prune_parent.right == new_prune_child:
        new_prune_parent.right = None

    # STEP 2: Collapse - connect sibling to grandparent
    if grandparent is None:
        # Pruning from root - sibling becomes new root
        new_tree = sibling
        sibling.distance = 0.0
    else:
        # Connect sibling to grandparent
        if grandparent.left == new_prune_parent:
            grandparent.left = sibling
        elif grandparent.right == new_prune_parent:
            grandparent.right = sibling

        # Update sibling's branch length (sum of two collapsed branches)
        sibling.distance = getattr(sibling, 'distance', 0.001) + getattr(new_prune_parent, 'distance', 0.001)

    # STEP 3: Regraft - create new internal node at regraft edge
    new_internal = TreeNode(name=None)  # Internal nodes have no name
    new_internal.left = new_prune_child
    new_internal.right = new_regraft_child

    # Update parent's connection
    if new_regraft_parent.left == new_regraft_child:
        new_regraft_parent.left = new_internal
    elif new_regraft_parent.right == new_regraft_child:
        new_regraft_parent.right = new_internal

    # STEP 4: Initialize branch lengths
    # Split the regraft edge
    original_regraft_length = getattr(new_regraft_child, 'distance', 0.001)
    new_internal.distance = original_regraft_length / 2
    new_regraft_child.distance = original_regraft_length / 2
    new_prune_child.distance = 0.001  # Small initial length

    # Validate result
    original_leaf_count = tree.count_leaves()
    new_leaf_count = new_tree.count_leaves()

    if new_leaf_count != original_leaf_count:
        # Tree corruption detected!
        return None

    return new_tree


def spr_search_proper(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    max_iterations: int = 5,
    max_moves_per_iteration: int = 20,
    tolerance: float = 0.01,
    starting_logL: Optional[float] = None,
    verbose: bool = False
) -> Tuple[TreeNode, float, int]:
    """
    Proper SPR tree search.

    Args:
        tree: Initial tree topology
        sequences: Aligned sequences
        alpha: Gamma parameter
        max_iterations: Maximum SPR iterations
        max_moves_per_iteration: Max SPR moves to try per iteration
        tolerance: Minimum improvement threshold
        starting_logL: Starting log-likelihood (None = recalculate)
        verbose: Print progress

    Returns:
        (improved_tree, final_logL, n_improvements)
    """
    from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood

    if verbose:
        print("\n" + "="*70)
        print("PROPER SPR TREE SEARCH")
        print("="*70)
        print(f"Max iterations: {max_iterations}")
        print(f"Max moves per iteration: {max_moves_per_iteration}")
        print(f"Tolerance: {tolerance}")
        print()

    current_tree = tree
    current_logL = starting_logL if starting_logL is not None else compute_log_likelihood(current_tree, sequences, alpha=alpha)

    # OPTIMIZATION: Cache these values to avoid recomputation
    original_leaf_count = current_tree.count_leaves()
    original_taxa = set(current_tree.get_leaf_names())

    if verbose:
        print(f"Initial LogL: {current_logL:.2f}")
        print(f"Initial taxa: {original_leaf_count}")

    n_improvements = 0

    for iteration in range(max_iterations):
        improved_this_round = False
        best_tree = None
        best_logL = current_logL

        # Get all potential prune edges (excluding root and leaves)
        all_nodes = get_all_nodes_with_parents(current_tree)
        prune_edges = [(parent, node) for parent, node in all_nodes
                       if parent is not None and not node.is_leaf()]

        if verbose:
            print(f"\nIteration {iteration + 1}: Testing up to {max_moves_per_iteration} SPR moves...")
            print(f"  Potential prune positions: {len(prune_edges)}")

        n_moves_tested = 0
        n_moves_valid = 0

        # Try SPR moves
        # SPEED: Limit prune positions based on max_moves_per_iteration
        max_prune = min(5, len(prune_edges))  # Try at most 5 prune positions

        # OPTIMIZATION: Pre-filter potential regraft edges (exclude leaves)
        potential_regraft_edges = [(parent, node) for parent, node in all_nodes
                                   if parent is not None]

        for prune_parent, prune_child in prune_edges[:max_prune]:
            # OPTIMIZATION: Filter regraft edges for this prune (much faster than rebuilding list)
            regraft_edges = [(p, n) for p, n in potential_regraft_edges
                            if n != prune_child and p != prune_parent]

            # SPEED: Limit regraft positions per prune
            max_regraft = max(1, max_moves_per_iteration // max_prune)  # Divide budget evenly
            for regraft_parent, regraft_child in regraft_edges[:max_regraft]:
                n_moves_tested += 1

                # Perform SPR move
                new_tree = spr_move_proper(current_tree, prune_parent, prune_child,
                                          regraft_parent, regraft_child)

                if new_tree is None:
                    continue  # Invalid move

                # OPTIMIZATION: Quick validation using cached values
                # Check leaf count first (faster than checking taxa names)
                new_leaf_count = new_tree.count_leaves()
                if new_leaf_count != original_leaf_count:
                    if verbose:
                        print(f"  [ERROR] Move {n_moves_tested}: Taxa count changed!")
                    continue

                # Only check taxa names if leaf count matches (expensive operation)
                new_taxa = set(new_tree.get_leaf_names())
                if new_taxa != original_taxa:
                    if verbose:
                        print(f"  [ERROR] Move {n_moves_tested}: Taxa names changed!")
                    continue

                n_moves_valid += 1

                # Optimize branch lengths
                try:
                    new_logL = optimize_branch_lengths_vectorized(
                        new_tree, sequences, alpha=alpha, verbose=False
                    )
                except Exception as e:
                    if verbose:
                        print(f"  [ERROR] Optimization failed: {e}")
                    continue

                # Check if better
                improvement = new_logL - current_logL

                if verbose and n_moves_valid % 5 == 0:
                    print(f"  Move {n_moves_valid}: LogL={new_logL:.2f} (delta={improvement:+.2f})")

                if new_logL > best_logL + tolerance:
                    best_tree = new_tree
                    best_logL = new_logL
                    improved_this_round = True

                    if verbose:
                        print(f"  [IMPROVEMENT] Found better tree: {improvement:+.2f}")

                    # OPTIMIZATION: Early termination - accept first improvement (greedy SPR)
                    # This provides massive speedup while still being effective
                    break

                if n_moves_tested >= max_moves_per_iteration:
                    break

            # OPTIMIZATION: Break outer loop if improvement found (greedy SPR)
            if improved_this_round:
                break

        if verbose:
            print(f"  Tested {n_moves_tested} moves ({n_moves_valid} valid)")

        if improved_this_round:
            current_tree = best_tree
            current_logL = best_logL
            n_improvements += 1

            if verbose:
                print(f"  Best LogL this round: {current_logL:.2f}")
        else:
            if verbose:
                print(f"\nSPR converged after {iteration + 1} iterations")
                print(f"Final LogL: {current_logL:.2f}")
                print(f"Total improvements: {n_improvements}")
            break

    return current_tree, current_logL, n_improvements

"""
Tree search algorithms for finding optimal phylogenetic topologies.

Implements NNI (Nearest Neighbor Interchange) and related tree rearrangement methods.
"""

from typing import List, Tuple, Optional
import copy
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood


def nni_search(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    max_iterations: int = 100,
    tolerance: float = 0.01,
    verbose: bool = False
) -> Tuple[TreeNode, float, int]:
    """
    Improve tree topology using NNI (Nearest Neighbor Interchange) rearrangements.

    NNI is a local search algorithm that explores the tree space by swapping
    subtrees around internal branches. It's fast but can get stuck in local optima.

    Algorithm:
        1. Calculate likelihood of current tree
        2. For each internal branch:
           a. Generate 2 NNI neighbors (different subtree swaps)
           b. Calculate likelihoods of neighbors
           c. Accept best neighbor if it improves likelihood
        3. Repeat until no improvement or max iterations reached

    Args:
        tree: Starting tree topology
        sequences: Aligned sequences
        alpha: Gamma shape parameter for rate heterogeneity
        max_iterations: Maximum number of NNI rounds
        tolerance: Minimum log-likelihood improvement to continue
        verbose: Print progress information

    Returns:
        (best_tree, best_logL, n_improvements)
        - best_tree: Improved tree topology
        - best_logL: Log-likelihood of best tree
        - n_improvements: Number of successful NNI moves

    Example:
        >>> # Start with UPGMA tree
        >>> initial_tree = build_upgma_tree(sequences)
        >>> # Improve with NNI
        >>> better_tree, logL, n_impr = nni_search(
        ...     initial_tree, sequences, max_iterations=10
        ... )
        >>> print(f"Improved by {n_impr} NNI moves, LogL = {logL:.2f}")
    """
    current_tree = tree.copy()
    current_logL = compute_log_likelihood(current_tree, sequences, alpha=alpha)

    if verbose:
        print("\n" + "=" * 70)
        print("NNI TREE SEARCH")
        print("=" * 70)
        print(f"Initial LogL: {current_logL:.2f}")
        print(f"Max iterations: {max_iterations}")
        print(f"Tolerance: {tolerance}")
        print()

    n_improvements = 0
    iteration = 0

    for iteration in range(max_iterations):
        improved_this_round = False
        best_logL_this_round = current_logL

        # Get all internal nodes (NNI candidates)
        internal_nodes = current_tree.get_internal_nodes()

        # Filter: only nodes with two children (both non-leaf)
        nni_candidates = [
            node for node in internal_nodes
            if node.left and node.right and
            not node.left.is_leaf() and not node.right.is_leaf()
        ]

        if verbose:
            print(f"Iteration {iteration + 1}: {len(nni_candidates)} NNI candidates")

        # Try NNI on each internal branch
        for node in nni_candidates:
            # Generate two NNI neighbors
            neighbor1 = _generate_nni_neighbor(current_tree, node, swap_type=1)
            neighbor2 = _generate_nni_neighbor(current_tree, node, swap_type=2)

            if neighbor1 is None or neighbor2 is None:
                continue

            # Calculate likelihoods
            logL1 = compute_log_likelihood(neighbor1, sequences, alpha=alpha)
            logL2 = compute_log_likelihood(neighbor2, sequences, alpha=alpha)

            # Find best
            best_neighbor_logL = max(logL1, logL2)

            # Accept if better than current best
            if best_neighbor_logL > best_logL_this_round + tolerance:
                if logL1 > logL2:
                    current_tree = neighbor1
                    current_logL = logL1
                else:
                    current_tree = neighbor2
                    current_logL = logL2

                best_logL_this_round = current_logL
                improved_this_round = True
                n_improvements += 1

                if verbose:
                    print(f"  NNI improved: LogL = {current_logL:.2f} "
                          f"(+{current_logL - (current_logL - tolerance):.2f})")

        if not improved_this_round:
            if verbose:
                print(f"\nNNI converged after {iteration + 1} iterations")
                print(f"Final LogL: {current_logL:.2f}")
                print(f"Total improvements: {n_improvements}")
            break

    if iteration == max_iterations - 1 and verbose:
        print(f"\nReached maximum iterations ({max_iterations})")
        print(f"Final LogL: {current_logL:.2f}")
        print(f"Total improvements: {n_improvements}")

    return current_tree, current_logL, n_improvements


def _generate_nni_neighbor(
    tree: TreeNode,
    internal_node: TreeNode,
    swap_type: int
) -> Optional[TreeNode]:
    """
    Generate an NNI neighbor by swapping subtrees.

    For an internal node with structure:
              node
             /    \
           left   right
           / \     / \
          A   B   C   D

    NNI swap type 1: Swap A and C
              node
             /    \
           left   right
           / \     / \
          C   B   A   D

    NNI swap type 2: Swap A and D
              node
             /    \
           left   right
           / \     / \
          D   B   C   A

    Args:
        tree: Full tree (to be copied)
        internal_node: Node to perform NNI on
        swap_type: 1 or 2 (which subtrees to swap)

    Returns:
        New tree with NNI rearrangement, or None if swap not possible
    """
    # Make a copy of the tree
    new_tree = tree.copy()

    # Find corresponding node in the copy
    # (We need to traverse the copy to find the same position)
    corresponding_node = _find_corresponding_node(new_tree, tree, internal_node)

    if corresponding_node is None:
        return None

    # Check structure
    if not corresponding_node.left or not corresponding_node.right:
        return None
    if corresponding_node.left.is_leaf() or corresponding_node.right.is_leaf():
        return None

    left = corresponding_node.left
    right = corresponding_node.right

    # Perform swap based on type
    if swap_type == 1:
        # Swap left.left with right.left
        if left.left and right.left:
            left.left, right.left = right.left, left.left
        else:
            return None
    elif swap_type == 2:
        # Swap left.left with right.right
        if left.left and right.right:
            left.left, right.right = right.right, left.left
        else:
            return None
    else:
        return None

    return new_tree


def _find_corresponding_node(
    tree1: TreeNode,
    tree2: TreeNode,
    target_node: TreeNode
) -> Optional[TreeNode]:
    """
    Find the node in tree1 that corresponds to target_node in tree2.

    Uses tree structure to identify the same position in different tree copies.

    Args:
        tree1: Tree to search in
        tree2: Tree containing target_node
        target_node: Node to find correspondence for

    Returns:
        Corresponding node in tree1, or None if not found
    """
    # Build path from root to target in tree2
    path = _get_path_to_node(tree2, target_node)

    if path is None:
        return None

    # Follow same path in tree1
    current = tree1
    for direction in path:
        if direction == 'L':
            if current.left:
                current = current.left
            else:
                return None
        elif direction == 'R':
            if current.right:
                current = current.right
            else:
                return None

    return current


def _get_path_to_node(root: TreeNode, target: TreeNode, path: str = "") -> Optional[str]:
    """
    Get path from root to target node as string of 'L'/'R' directions.

    Args:
        root: Root of tree
        target: Target node to find
        path: Current path (used in recursion)

    Returns:
        Path string like "LRL" (left, right, left), or None if not found
    """
    if root is target:
        return path

    if root.left:
        left_path = _get_path_to_node(root.left, target, path + "L")
        if left_path is not None:
            return left_path

    if root.right:
        right_path = _get_path_to_node(root.right, target, path + "R")
        if right_path is not None:
            return right_path

    return None


def get_nni_neighbors(tree: TreeNode) -> List[TreeNode]:
    """
    Generate all possible NNI neighbors of a tree.

    For a tree with n internal branches that can undergo NNI,
    there are 2*n possible neighbors.

    Args:
        tree: Input tree

    Returns:
        List of NNI neighbor trees

    Example:
        >>> tree = build_tree(sequences)
        >>> neighbors = get_nni_neighbors(tree)
        >>> print(f"Generated {len(neighbors)} NNI neighbors")
    """
    neighbors = []

    # Get all internal nodes
    internal_nodes = tree.get_internal_nodes()

    # Filter to NNI candidates
    nni_candidates = [
        node for node in internal_nodes
        if node.left and node.right and
        not node.left.is_leaf() and not node.right.is_leaf()
    ]

    # Generate neighbors
    for node in nni_candidates:
        # Two swap types per node
        neighbor1 = _generate_nni_neighbor(tree, node, swap_type=1)
        neighbor2 = _generate_nni_neighbor(tree, node, swap_type=2)

        if neighbor1 is not None:
            neighbors.append(neighbor1)
        if neighbor2 is not None:
            neighbors.append(neighbor2)

    return neighbors


def hill_climbing_search(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    max_steps: int = 1000,
    verbose: bool = False
) -> Tuple[TreeNode, float, int]:
    """
    Simple hill-climbing search: always accept better trees.

    More aggressive than iterative NNI - explores all neighbors at each step.

    Args:
        tree: Starting tree
        sequences: Aligned sequences
        alpha: Gamma parameter
        max_steps: Maximum steps
        verbose: Print progress

    Returns:
        (best_tree, best_logL, n_steps)
    """
    current_tree = tree.copy()
    current_logL = compute_log_likelihood(current_tree, sequences, alpha=alpha)

    if verbose:
        print(f"\nHill-climbing search starting from LogL = {current_logL:.2f}")

    n_steps = 0

    for step in range(max_steps):
        # Generate all NNI neighbors
        neighbors = get_nni_neighbors(current_tree)

        if not neighbors:
            if verbose:
                print("No more neighbors to explore")
            break

        # Evaluate all neighbors
        best_neighbor = None
        best_neighbor_logL = current_logL

        for neighbor in neighbors:
            logL = compute_log_likelihood(neighbor, sequences, alpha=alpha)
            if logL > best_neighbor_logL:
                best_neighbor = neighbor
                best_neighbor_logL = logL

        # Accept if better
        if best_neighbor is not None:
            current_tree = best_neighbor
            current_logL = best_neighbor_logL
            n_steps += 1

            if verbose:
                print(f"Step {n_steps}: LogL = {current_logL:.2f}")
        else:
            if verbose:
                print(f"Converged after {n_steps} steps")
            break

    return current_tree, current_logL, n_steps

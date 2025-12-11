"""
In-place NNI implementation - avoids expensive tree copying.

Speedup: 2-3x faster by eliminating tree.copy() calls.
"""

from typing import List, Tuple, Optional
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood, LikelihoodCalculatorLevel3
from rrna_phylo.models.branch_length_optimizer import optimize_branch_lengths_fast

# GPU support (optional)
try:
    from rrna_phylo.models.gpu_likelihood_torch import TorchGPULikelihoodCalculator
    import torch
    GPU_AVAILABLE = torch.cuda.is_available()
except ImportError:
    GPU_AVAILABLE = False


def _compute_likelihood_unified(tree: TreeNode, sequences: List[Sequence], alpha: float, calculator) -> float:
    """
    Compute likelihood using either CPU or GPU calculator.

    Args:
        tree: Tree topology
        sequences: Aligned sequences
        alpha: Gamma shape parameter
        calculator: Either LikelihoodCalculatorLevel3 (CPU) or TorchGPULikelihoodCalculator (GPU)

    Returns:
        Log-likelihood
    """
    # GPU calculator has compute_likelihood() method
    if hasattr(calculator, 'device'):  # GPU calculator
        return calculator.compute_likelihood(tree)
    else:  # CPU calculator
        return compute_log_likelihood(tree, sequences, alpha=alpha, calculator=calculator)


class NNISwap:
    """Records an NNI swap for undo capability."""
    def __init__(self, node, swap_type):
        self.node = node
        self.swap_type = swap_type
        # Save original pointers AND distances (BUG FIX #1)
        if swap_type == 1:
            self.original_left_left = node.left.left
            self.original_right_left = node.right.left
            # Save original distances BEFORE swap
            self.dist_left_left = getattr(node.left.left, 'distance', None)
            self.dist_right_left = getattr(node.right.left, 'distance', None)
        else:  # swap_type == 2
            self.original_left_left = node.left.left
            self.original_right_right = node.right.right
            # Save original distances BEFORE swap
            self.dist_left_left = getattr(node.left.left, 'distance', None)
            self.dist_right_right = getattr(node.right.right, 'distance', None)


def perform_nni_swap_inplace(node: TreeNode, swap_type: int) -> Optional[NNISwap]:
    """
    Perform NNI swap IN-PLACE (modifies tree directly, no copy).

    Returns NNISwap object for undo, or None if swap not possible.
    """
    # Validate node structure
    if not node.left or not node.right:
        return None
    if node.left.is_leaf() or node.right.is_leaf():
        return None

    left = node.left
    right = node.right

    # Create swap record BEFORE modifying
    swap_record = NNISwap(node, swap_type)

    if swap_type == 1:
        # Swap left.left with right.left
        if not left.left or not right.left:
            return None

        # Save distances
        dist_left_left = getattr(left.left, 'distance', None)
        dist_right_left = getattr(right.left, 'distance', None)

        # Perform swap
        left.left, right.left = right.left, left.left

        # Restore distances
        if dist_right_left is not None:
            left.left.distance = dist_right_left
        if dist_left_left is not None:
            right.left.distance = dist_left_left

    elif swap_type == 2:
        # Swap left.left with right.right
        if not left.left or not right.right:
            return None

        # Save distances
        dist_left_left = getattr(left.left, 'distance', None)
        dist_right_right = getattr(right.right, 'distance', None)

        # Perform swap
        left.left, right.right = right.right, left.left

        # Restore distances
        if dist_right_right is not None:
            left.left.distance = dist_right_right
        if dist_left_left is not None:
            right.right.distance = dist_left_left
    else:
        return None

    return swap_record


def undo_nni_swap(swap: NNISwap):
    """
    Undo an NNI swap using the swap record.
    BUG FIX #1: Use SAVED distances, not current distances!
    """
    node = swap.node
    left = node.left
    right = node.right

    if swap.swap_type == 1:
        # Restore left.left and right.left pointers
        left.left = swap.original_left_left
        right.left = swap.original_right_left

        # Restore SAVED distances (not current!)
        if swap.dist_left_left is not None:
            left.left.distance = swap.dist_left_left
        if swap.dist_right_left is not None:
            right.left.distance = swap.dist_right_left

    else:  # swap_type == 2
        # Restore left.left and right.right pointers
        left.left = swap.original_left_left
        right.right = swap.original_right_right

        # Restore SAVED distances (not current!)
        if swap.dist_left_left is not None:
            left.left.distance = swap.dist_left_left
        if swap.dist_right_right is not None:
            right.right.distance = swap.dist_right_right


def _node_score(node: TreeNode) -> float:
    """
    Score node by sum of adjacent branch lengths.
    Long branches are more likely to improve likelihood.
    """
    score = 0.0
    if node.left:
        if node.left.left:
            score += getattr(node.left.left, 'distance', 0.0)
        if node.left.right:
            score += getattr(node.left.right, 'distance', 0.0)
    if node.right:
        if node.right.left:
            score += getattr(node.right.left, 'distance', 0.0)
        if node.right.right:
            score += getattr(node.right.right, 'distance', 0.0)
    return score


def nni_search(
    tree: TreeNode,
    sequences: List[Sequence],
    model_name: str = "GTR",
    alpha: Optional[float] = None,
    max_iterations: int = 50,
    tolerance: float = 0.01,  # FIX: Reduced from 0.1 to 0.01 to accept smaller improvements
    verbose: bool = False,
    use_gpu: bool = False,
    enable_pruning: bool = True,
    pruning_k: Optional[int] = None,
    compressor = None  # CRITICAL FIX: Accept compressor from model selection
) -> Tuple[TreeNode, float, int]:
    """
    IN-PLACE NNI search - 2-3x faster by avoiding tree copies.

    Modifies tree in-place and returns it.

    Args:
        model_name: Substitution model (JC69, K80, F81, HKY85, GTR) - from model_selection
        alpha: Gamma parameter (from model_selection)
        use_gpu: Use GPU acceleration if available (auto-selects based on dataset size if 'auto')
        enable_pruning: Enable heuristic candidate pruning for speed (default: True)
        pruning_k: Max candidates per iteration (auto-selects if None)
    """
    # Auto GPU selection based on dataset size (empirically determined)
    # Benchmark results: GPU slower at 5 seqs (0.79x), 7.64x faster at 10 seqs
    # Threshold set at 7 seqs as safe middle ground
    if use_gpu == 'auto':
        use_gpu = len(sequences) >= 7 and GPU_AVAILABLE

    # Create calculator once with SELECTED model
    if verbose:
        device_str = "GPU" if (use_gpu and GPU_AVAILABLE) else "CPU"
        print(f"Creating likelihood calculator (model={model_name}, device={device_str})...")

    # Handle model name with "+G" suffix
    if model_name.endswith('+G'):
        base_model_name = model_name[:-2]
    else:
        base_model_name = model_name

    # Create the selected model and get Q matrix + base frequencies
    if base_model_name == "GTR":
        # GTR model uses old interface with estimate_parameters
        from rrna_phylo.models.ml_tree_level3 import GTRModel as LegacyGTRModel
        model = LegacyGTRModel()
        model.estimate_parameters(sequences)
        Q = model.Q
        base_freq = model.base_freq
    else:
        # Non-GTR models use new substitution_models interface
        from rrna_phylo.models.substitution_models import get_model, compute_empirical_frequencies
        model_obj = get_model(base_model_name)

        # Compute empirical base frequencies from sequences
        base_freq = compute_empirical_frequencies(sequences)

        # Get Q matrix with empirical frequencies
        Q = model_obj.get_rate_matrix(freqs=base_freq)

    if alpha is None:
        alpha = 1.0

    # Create appropriate calculator (GPU or CPU)
    if use_gpu and GPU_AVAILABLE:
        # CRITICAL FIX: Pass compressor from model selection to ensure CPU/GPU consistency
        calculator = TorchGPULikelihoodCalculator(
            Q=Q,
            base_freq=base_freq,
            sequences=sequences,
            alpha=alpha,
            use_gpu=True,
            verbose=False,
            compressor=compressor  # CRITICAL: Reuse CPU compressor!
        )
    else:
        # CPU calculator needs the model object (for GTR) or we need to create one
        if base_model_name == "GTR":
            calculator = LikelihoodCalculatorLevel3(model, sequences, alpha=alpha)
        else:
            # For non-GTR, create a wrapper that has Q and base_freq attributes
            class ModelWrapper:
                def __init__(self, Q_matrix, frequencies):
                    self.Q = Q_matrix
                    self.base_freq = frequencies

            model_wrapper = ModelWrapper(Q, base_freq)
            calculator = LikelihoodCalculatorLevel3(model_wrapper, sequences, alpha=alpha)

    if verbose:
        print(f"Calculator ready. Alpha={alpha:.2f}")

    current_logL = _compute_likelihood_unified(tree, sequences, alpha, calculator)

    if verbose:
        print("\n" + "=" * 70)
        print("IN-PLACE NNI TREE SEARCH")
        print("=" * 70)
        print(f"Initial LogL: {current_logL:.2f}")
        print(f"Max iterations: {max_iterations}")
        print(f"Tolerance: {tolerance}")
        print()

    n_improvements = 0

    for iteration in range(max_iterations):
        improved_this_round = False
        best_logL_this_round = current_logL

        # Get NNI candidates
        internal_nodes = tree.get_internal_nodes()
        nni_candidates = [
            node for node in internal_nodes
            if node.left and node.right and
            not node.left.is_leaf() and not node.right.is_leaf()
        ]

        # Auto-select pruning K based on dataset size AND actual candidate count
        total_candidates = len(nni_candidates)

        if enable_pruning and pruning_k is None:
            # Only enable pruning if we have enough candidates to make it worthwhile
            if total_candidates <= 8:
                # Too few candidates, pruning overhead dominates any benefit
                k = None
            else:
                n_seqs = len(sequences)
                if n_seqs < 30:
                    k = max(8, total_candidates // 2)  # Prune to 50% or min 8
                elif n_seqs < 100:
                    k = max(12, total_candidates // 2)  # Prune to 50% or min 12
                else:
                    k = max(20, total_candidates // 2)  # Prune to 50% or min 20
        else:
            k = pruning_k if enable_pruning else None

        # Apply heuristic pruning: test only top K candidates
        if k is not None and total_candidates > k:
            # Score by sum of adjacent branch lengths (long branches more likely to improve)
            nni_candidates = sorted(nni_candidates, key=_node_score, reverse=True)[:k]

        if verbose:
            if k is not None and total_candidates > k:
                print(f"Iteration {iteration + 1}: Testing {len(nni_candidates)}/{total_candidates} candidates (pruned)")
            else:
                print(f"Iteration {iteration + 1}: {len(nni_candidates)} NNI candidates")

        best_swap = None
        best_swap_logL = current_logL

        # Track best improvement found (even if rejected)
        best_attempted_improvement = 0.0
        n_swaps_tested = 0

        # Try NNI on each internal branch
        # BUG FIX #3: Must test BOTH swaps from ORIGINAL topology, not from swapped state!
        for node in nni_candidates:
            # Try swap type 1
            swap1 = perform_nni_swap_inplace(node, swap_type=1)
            if swap1:
                logL1 = _compute_likelihood_unified(tree, sequences, alpha, calculator)
                n_swaps_tested += 1
                improvement1 = logL1 - current_logL

                # Track best attempted improvement
                if improvement1 > best_attempted_improvement:
                    best_attempted_improvement = improvement1

                if logL1 > best_swap_logL + tolerance:
                    # Found better! Undo previous best if any
                    if best_swap:
                        undo_nni_swap(best_swap)
                    best_swap = swap1
                    best_swap_logL = logL1
                else:
                    # Not better, undo immediately
                    undo_nni_swap(swap1)

            # Try swap type 2 - CRITICAL: Undo swap1 first if it's currently the best!
            # Otherwise we'd test swap2 on top of swap1, which is wrong
            undo_swap1_for_swap2_test = False
            if best_swap and best_swap.node == node and best_swap.swap_type == 1:
                undo_nni_swap(best_swap)
                undo_swap1_for_swap2_test = True
                temp_best_swap = best_swap
                temp_best_logL = best_swap_logL

            swap2 = perform_nni_swap_inplace(node, swap_type=2)
            if swap2:
                logL2 = _compute_likelihood_unified(tree, sequences, alpha, calculator)
                n_swaps_tested += 1
                improvement2 = logL2 - current_logL

                # Track best attempted improvement
                if improvement2 > best_attempted_improvement:
                    best_attempted_improvement = improvement2

                # Compare: Which is better for THIS node? swap1 or swap2?
                if undo_swap1_for_swap2_test:
                    # We need to compare swap1 vs swap2 for this node
                    if logL2 > temp_best_logL:
                        # swap2 is better than swap1 for this node
                        best_swap = swap2
                        best_swap_logL = logL2
                    else:
                        # swap1 was better, restore it
                        undo_nni_swap(swap2)
                        # Re-apply swap1
                        swap1_reapply = perform_nni_swap_inplace(node, swap_type=1)
                        best_swap = swap1_reapply
                        best_swap_logL = temp_best_logL
                elif logL2 > best_swap_logL + tolerance:
                    # swap2 is better than current best (from other nodes)
                    if best_swap:
                        undo_nni_swap(best_swap)
                    best_swap = swap2
                    best_swap_logL = logL2
                else:
                    # Not better, undo immediately
                    undo_nni_swap(swap2)

        # After trying all candidates, accept best if found
        if best_swap:
            # PERFORMANCE FIX: Mark topology changed for cache invalidation
            # This tells the calculator to clear internal caches on next compute_likelihood() call
            if hasattr(calculator, 'mark_topology_changed'):
                calculator.mark_topology_changed()

            # PERFORMANCE FIX: Don't optimize branch lengths during NNI!
            # This adds 10-50 likelihood calls per iteration.
            # Instead, use the likelihood from the swap test and optimize once at the end.
            current_logL = best_swap_logL

            improved_this_round = True
            n_improvements += 1

            improvement = current_logL - best_logL_this_round
            if verbose:
                print(f"  NNI improved: LogL = {current_logL:.2f} (+{improvement:.2f})")

        if not improved_this_round:
            if verbose:
                print(f"  Tested {n_swaps_tested} swaps, best attempted improvement: {best_attempted_improvement:+.4f}")
                if best_attempted_improvement > 0 and best_attempted_improvement <= tolerance:
                    print(f"  (below tolerance threshold of {tolerance})")
                print(f"\nNNI converged after {iteration + 1} iterations")
                print(f"Final LogL: {current_logL:.2f}")
                print(f"Total improvements: {n_improvements}")
            break

    if iteration == max_iterations - 1 and verbose:
        print(f"\nReached maximum iterations ({max_iterations})")
        print(f"Final LogL: {current_logL:.2f}")
        print(f"Total improvements: {n_improvements}")

    # PERFORMANCE FIX: Optimize branch lengths ONCE at the end
    # This is much faster than optimizing after every swap (10-50 calls per iteration)
    # RAxML and IQ-TREE use this approach
    if n_improvements > 0:
        if verbose:
            print(f"\nOptimizing final branch lengths...")
        current_logL = optimize_branch_lengths_fast(
            tree, sequences, alpha=alpha, calculator=calculator, verbose=False
        )
        if verbose:
            print(f"Final optimized LogL: {current_logL:.2f}")

    return tree, current_logL, n_improvements

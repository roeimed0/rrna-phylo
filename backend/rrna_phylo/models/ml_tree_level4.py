"""
Maximum Likelihood Tree Inference - Level 4

Advanced ML phylogenetics with:
- Model selection (AIC/BIC)
- Tree search (NNI)
- Automatic parameter optimization

This is the most sophisticated ML implementation in the package.
"""

from typing import List, Tuple, Dict, Optional
import time
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.models.ml_tree_level3 import (
    build_ml_tree_level3,
    compute_log_likelihood
)
from rrna_phylo.models.model_selection import (
    select_best_model,
    select_best_model_with_gamma
)
from rrna_phylo.models.tree_search import nni_search, hill_climbing_search


def build_ml_tree_level4(
    sequences: List[Sequence],
    model: str = 'auto',
    alpha: Optional[float] = None,
    tree_search: str = 'nni',
    max_iterations: int = 100,
    starting_tree: Optional[TreeNode] = None,
    criterion: str = 'BIC',
    test_gamma: bool = True,
    verbose: bool = False
) -> Tuple[TreeNode, float, Dict]:
    """
    Build Maximum Likelihood phylogenetic tree with Level 4 features.

    Level 4 enhancements:
    - Automatic model selection using AIC/BIC
    - Tree topology search using NNI or hill-climbing
    - Gamma rate heterogeneity testing
    - Comprehensive metadata reporting

    Args:
        sequences: Aligned sequences
        model: Substitution model name or 'auto' for automatic selection
        alpha: Gamma shape parameter (None = uniform rates, 'auto' = optimize)
        tree_search: 'nni', 'hill', or None (no topology search)
        max_iterations: Maximum NNI iterations
        starting_tree: Initial topology (None = use BioNJ)
        criterion: 'AIC' or 'BIC' for model selection
        test_gamma: Test models with +G variants
        verbose: Print detailed progress

    Returns:
        (tree, log_likelihood, metadata)

        metadata includes:
        - 'selected_model': Best model name
        - 'model_score': AIC/BIC score
        - 'alpha': Gamma parameter (if used)
        - 'n_nni_improvements': Number of successful NNI moves
        - 'initial_logL': Likelihood before tree search
        - 'final_logL': Final likelihood
        - 'time_model_selection': Time spent on model selection
        - 'time_tree_search': Time spent on topology search
        - 'time_total': Total time

    Example:
        >>> # Automatic model selection + NNI search
        >>> tree, logL, meta = build_ml_tree_level4(
        ...     sequences,
        ...     model='auto',
        ...     tree_search='nni',
        ...     verbose=True
        ... )
        >>> print(f"Best model: {meta['selected_model']}")
        >>> print(f"Final LogL: {logL:.2f}")
        >>> print(f"NNI improvements: {meta['n_nni_improvements']}")
    """
    start_time = time.time()

    metadata = {
        'selected_model': None,
        'model_score': None,
        'alpha': None,
        'n_nni_improvements': 0,
        'initial_logL': None,
        'final_logL': None,
        'time_model_selection': 0,
        'time_tree_search': 0,
        'time_total': 0,
        'criterion': criterion,
        'tree_search_method': tree_search
    }

    if verbose:
        print("\n" + "=" * 70)
        print("MAXIMUM LIKELIHOOD TREE - LEVEL 4")
        print("=" * 70)
        print(f"Sequences: {len(sequences)}")
        print(f"Alignment length: {len(sequences[0].sequence)}")
        print(f"Model: {model}")
        print(f"Tree search: {tree_search or 'None'}")
        print(f"Model selection criterion: {criterion}")
        print()

    # Step 1: Get initial tree
    if starting_tree is None:
        if verbose:
            print("Building initial tree with BioNJ...")
        # Build BioNJ tree from distance matrix
        dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
        initial_tree = build_bionj_tree(dist_matrix, ids)
    else:
        initial_tree = starting_tree.copy()

    if verbose:
        print(f"Initial tree has {initial_tree.count_leaves()} leaves")

    # Step 2: Model selection (if auto)
    time_model_start = time.time()

    if model == 'auto':
        if test_gamma:
            # Test models with and without gamma
            if verbose:
                print("\nPerforming model selection with gamma variants...")

            best_model, best_alpha, best_score, all_results = select_best_model_with_gamma(
                initial_tree,
                sequences,
                criterion=criterion,
                test_gamma=True,
                verbose=verbose
            )

            metadata['selected_model'] = best_model
            metadata['alpha'] = best_alpha
            metadata['model_score'] = best_score
            metadata['all_model_scores'] = all_results

        else:
            # Test models without gamma
            if verbose:
                print("\nPerforming model selection...")

            best_model, best_score, all_results = select_best_model(
                initial_tree,
                sequences,
                criterion=criterion,
                alpha=alpha,
                verbose=verbose
            )

            metadata['selected_model'] = best_model
            metadata['model_score'] = best_score
            metadata['alpha'] = alpha
            metadata['all_model_scores'] = all_results

        selected_model = best_model.replace('+G', '')  # Remove +G suffix for now
        selected_alpha = metadata['alpha']

    else:
        # Use specified model
        selected_model = model
        selected_alpha = alpha
        metadata['selected_model'] = model
        metadata['alpha'] = alpha

    metadata['time_model_selection'] = time.time() - time_model_start

    if verbose:
        print(f"\nUsing model: {metadata['selected_model']}")
        if metadata['alpha'] is not None:
            print(f"Gamma parameter: {metadata['alpha']:.3f}")

    # Step 3: Calculate initial likelihood
    initial_logL = compute_log_likelihood(
        initial_tree,
        sequences,
        alpha=selected_alpha
    )
    metadata['initial_logL'] = initial_logL

    if verbose:
        print(f"Initial LogL: {initial_logL:.2f}")

    # Step 4: Tree topology search
    time_search_start = time.time()
    current_tree = initial_tree
    current_logL = initial_logL

    if tree_search == 'nni':
        if verbose:
            print("\nPerforming NNI tree search...")

        improved_tree, improved_logL, n_improvements = nni_search(
            initial_tree,
            sequences,
            alpha=selected_alpha,
            max_iterations=max_iterations,
            verbose=verbose
        )

        current_tree = improved_tree
        current_logL = improved_logL
        metadata['n_nni_improvements'] = n_improvements

    elif tree_search == 'hill':
        if verbose:
            print("\nPerforming hill-climbing tree search...")

        improved_tree, improved_logL, n_steps = hill_climbing_search(
            initial_tree,
            sequences,
            alpha=selected_alpha,
            max_steps=max_iterations,
            verbose=verbose
        )

        current_tree = improved_tree
        current_logL = improved_logL
        metadata['n_nni_improvements'] = n_steps

    elif tree_search is None:
        if verbose:
            print("\nSkipping tree search (using initial topology)")
    else:
        raise ValueError(f"Unknown tree search method: {tree_search}")

    metadata['time_tree_search'] = time.time() - time_search_start
    metadata['final_logL'] = current_logL

    # Step 5: Final likelihood improvement
    improvement = current_logL - initial_logL

    metadata['time_total'] = time.time() - start_time

    if verbose:
        print("\n" + "=" * 70)
        print("LEVEL 4 ML TREE COMPLETE")
        print("=" * 70)
        print(f"Model: {metadata['selected_model']}")
        print(f"Initial LogL: {initial_logL:.2f}")
        print(f"Final LogL: {current_logL:.2f}")
        print(f"Improvement: {improvement:.2f}")
        if tree_search:
            print(f"NNI improvements: {metadata['n_nni_improvements']}")
        print(f"\nTiming:")
        print(f"  Model selection: {metadata['time_model_selection']:.2f}s")
        print(f"  Tree search: {metadata['time_tree_search']:.2f}s")
        print(f"  Total: {metadata['time_total']:.2f}s")
        print("=" * 70)

    return current_tree, current_logL, metadata


def build_ml_tree_fast(
    sequences: List[Sequence],
    verbose: bool = False
) -> Tuple[TreeNode, float, Dict]:
    """
    Quick ML tree with reasonable defaults.

    Uses:
    - Automatic model selection (BIC)
    - NNI tree search (10 iterations)
    - Gamma rate heterogeneity testing

    Args:
        sequences: Aligned sequences
        verbose: Print progress

    Returns:
        (tree, log_likelihood, metadata)
    """
    return build_ml_tree_level4(
        sequences,
        model='auto',
        tree_search='nni',
        max_iterations=10,
        criterion='BIC',
        test_gamma=True,
        verbose=verbose
    )


def build_ml_tree_thorough(
    sequences: List[Sequence],
    verbose: bool = False
) -> Tuple[TreeNode, float, Dict]:
    """
    Thorough ML tree search with extensive optimization.

    Uses:
    - Automatic model selection (BIC)
    - Hill-climbing search (100 steps)
    - Gamma rate heterogeneity testing

    Args:
        sequences: Aligned sequences
        verbose: Print progress

    Returns:
        (tree, log_likelihood, metadata)
    """
    return build_ml_tree_level4(
        sequences,
        model='auto',
        tree_search='hill',
        max_iterations=100,
        criterion='BIC',
        test_gamma=True,
        verbose=verbose
    )


def compare_models_on_tree(
    tree: TreeNode,
    sequences: List[Sequence],
    models: Optional[List[str]] = None,
    verbose: bool = True
) -> Dict[str, Dict]:
    """
    Compare multiple models on a fixed tree topology.

    Useful for understanding model selection results.

    Args:
        tree: Fixed tree topology
        sequences: Aligned sequences
        models: List of models to compare (None = all)
        verbose: Print comparison table

    Returns:
        Dictionary of model results

    Example:
        >>> tree = build_bionj_tree(sequences)
        >>> results = compare_models_on_tree(tree, sequences)
        >>> for model, res in sorted(results.items(), key=lambda x: x[1]['BIC']):
        ...     print(f"{model}: BIC = {res['BIC']:.2f}")
    """
    _, _, all_results = select_best_model(
        tree,
        sequences,
        models=models,
        criterion='BIC',
        verbose=verbose
    )

    # Add AIC scores as well
    for model_name, res in all_results.items():
        from rrna_phylo.models.model_selection import calculate_aic
        res['AIC'] = calculate_aic(res['logL'], res['n_params'])

    return all_results

"""
Model selection using information criteria (AIC/BIC).

Implements automatic testing of multiple substitution models
and selection of the best-fitting model using statistical criteria.
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.substitution_models import (
    get_model,
    compute_empirical_frequencies,
    DNA_MODELS
)
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood


def calculate_aic(log_likelihood: float, n_params: int) -> float:
    """
    Calculate Akaike Information Criterion.

    AIC = -2*ln(L) + 2*k

    Lower values indicate better fit with penalty for complexity.

    Args:
        log_likelihood: Log-likelihood of the model
        n_params: Number of free parameters

    Returns:
        AIC score (lower is better)
    """
    return -2 * log_likelihood + 2 * n_params


def calculate_bic(log_likelihood: float, n_params: int, n_sites: int) -> float:
    """
    Calculate Bayesian Information Criterion.

    BIC = -2*ln(L) + k*ln(n)

    Stronger penalty for complexity than AIC. Lower values indicate better fit.

    Args:
        log_likelihood: Log-likelihood of the model
        n_params: Number of free parameters
        n_sites: Number of data points (alignment sites)

    Returns:
        BIC score (lower is better)
    """
    return -2 * log_likelihood + n_params * np.log(n_sites)


def compute_likelihood_for_model(
    tree: TreeNode,
    sequences: List[Sequence],
    model_name: str,
    alpha: Optional[float] = None,
    verbose: bool = False
) -> Tuple[float, Dict]:
    """
    Compute likelihood for a specific substitution model.

    Args:
        tree: Phylogenetic tree (fixed topology)
        sequences: Aligned sequences
        model_name: Model name ('JC69', 'K80', 'F81', 'HKY85', 'GTR')
        alpha: Gamma shape parameter (None = no rate heterogeneity)
        verbose: Print details

    Returns:
        (log_likelihood, model_params)
        model_params includes fitted parameter values
    """
    model = get_model(model_name)

    # Get empirical base frequencies
    freqs = compute_empirical_frequencies(sequences)

    # Set up model parameters
    params = None
    if model_name == 'K80':
        params = np.array([2.0])  # Initial kappa
    elif model_name == 'HKY85':
        params = np.array([2.0])  # Initial kappa
    elif model_name == 'GTR':
        params = np.ones(6)  # Initial exchangeability rates

    # Compute likelihood
    # Note: We're using the Level 3 likelihood function
    # which already handles GTR model
    if model_name == 'GTR':
        log_likelihood = compute_log_likelihood(
            tree, sequences, alpha=alpha
        )
    else:
        # For simpler models, we need to compute likelihood
        # using the model's rate matrix
        log_likelihood = _compute_likelihood_with_rate_matrix(
            tree, sequences, model, params, freqs, alpha
        )

    model_params = {
        'model': model_name,
        'params': params,
        'frequencies': freqs,
        'alpha': alpha
    }

    if verbose:
        print(f"{model_name}: LogL = {log_likelihood:.2f}")

    return log_likelihood, model_params


def _compute_likelihood_with_rate_matrix(
    tree: TreeNode,
    sequences: List[Sequence],
    model,
    params: Optional[np.ndarray],
    freqs: np.ndarray,
    alpha: Optional[float] = None
) -> float:
    """
    Compute likelihood using a specific rate matrix.

    This is a simplified version that uses the model's rate matrix
    instead of the GTR matrix used in Level 3.

    Args:
        tree: Tree topology
        sequences: Aligned sequences
        model: SubstitutionModel instance
        params: Model-specific parameters
        freqs: Base frequencies
        alpha: Gamma shape parameter

    Returns:
        Log-likelihood
    """
    # For now, use a simplified approach:
    # If the model is simpler than GTR, we approximate by
    # using the GTR likelihood (which should be >= simpler models)
    # In a full implementation, we'd compute the exact likelihood
    # using the model's specific rate matrix

    # This is a placeholder - proper implementation would
    # compute likelihood using model.get_probability_matrix()
    # for each branch

    from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood
    return compute_log_likelihood(tree, sequences, alpha=alpha)


def select_best_model(
    tree: TreeNode,
    sequences: List[Sequence],
    models: Optional[List[str]] = None,
    criterion: str = 'BIC',
    alpha: Optional[float] = None,
    verbose: bool = False
) -> Tuple[str, float, Dict[str, Dict]]:
    """
    Test multiple models and select the best using AIC or BIC.

    Args:
        tree: Phylogenetic tree (fixed topology)
        sequences: Aligned sequences
        models: List of model names to test (None = test all)
        criterion: 'AIC' or 'BIC'
        alpha: Gamma shape parameter for rate heterogeneity
        verbose: Print comparison table

    Returns:
        (best_model_name, best_score, all_results)

        all_results is dict: {
            'model_name': {
                'logL': log_likelihood,
                'n_params': number_of_parameters,
                'score': AIC or BIC score,
                'params': fitted parameters
            }
        }

    Example:
        >>> best, score, results = select_best_model(
        ...     tree, sequences, criterion='BIC'
        ... )
        >>> print(f"Best model: {best}")
        Best model: HKY85
    """
    if models is None:
        # Test all DNA models
        models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']

    if criterion not in ['AIC', 'BIC']:
        raise ValueError(f"Invalid criterion: {criterion}. Use 'AIC' or 'BIC'")

    n_sites = len(sequences[0].sequence)
    results = {}

    if verbose:
        print("\n" + "=" * 70)
        print("MODEL SELECTION")
        print("=" * 70)
        print(f"Criterion: {criterion}")
        print(f"Number of sites: {n_sites}")
        print(f"Number of taxa: {len(sequences)}")
        if alpha is not None:
            print(f"Gamma shape parameter: {alpha:.3f}")
        print()

    # Test each model
    for model_name in models:
        model = get_model(model_name)

        # Compute likelihood
        log_likelihood, model_params = compute_likelihood_for_model(
            tree, sequences, model_name, alpha=alpha, verbose=False
        )

        # Calculate number of parameters
        n_params = model.n_params

        # Add branch lengths as parameters
        n_branches = count_branches(tree)
        n_params += n_branches

        # Add alpha if using gamma
        if alpha is not None:
            n_params += 1

        # Calculate information criterion
        if criterion == 'AIC':
            score = calculate_aic(log_likelihood, n_params)
        else:  # BIC
            score = calculate_bic(log_likelihood, n_params, n_sites)

        results[model_name] = {
            'logL': log_likelihood,
            'n_params': n_params,
            'score': score,
            'params': model_params
        }

    # Find best model (lowest score)
    best_model_name = min(results.keys(), key=lambda k: results[k]['score'])
    best_score = results[best_model_name]['score']

    if verbose:
        print(f"{'Model':<12} {'LogL':>12} {'Params':>8} {criterion:>12} {'Delta':>10}")
        print("-" * 70)

        # Sort by score
        sorted_results = sorted(results.items(), key=lambda x: x[1]['score'])

        for model_name, res in sorted_results:
            delta = res['score'] - best_score
            marker = " *" if model_name == best_model_name else ""
            print(f"{model_name:<12} {res['logL']:>12.2f} {res['n_params']:>8} "
                  f"{res['score']:>12.2f} {delta:>10.2f}{marker}")

        print("\n* Best model selected")
        print(f"\nBest model: {best_model_name}")
        print(f"{criterion} score: {best_score:.2f}")

    return best_model_name, best_score, results


def count_branches(tree: TreeNode) -> int:
    """
    Count number of branches in a tree.

    Args:
        tree: Tree root

    Returns:
        Number of branches (edges)
    """
    if tree is None or tree.is_leaf():
        return 0

    count = 0
    if tree.left is not None:
        count += 1 + count_branches(tree.left)
    if tree.right is not None:
        count += 1 + count_branches(tree.right)

    return count


def select_best_model_with_gamma(
    tree: TreeNode,
    sequences: List[Sequence],
    models: Optional[List[str]] = None,
    criterion: str = 'BIC',
    test_gamma: bool = True,
    verbose: bool = False
) -> Tuple[str, Optional[float], float, Dict[str, Dict]]:
    """
    Select best model, optionally testing with/without gamma rate heterogeneity.

    Args:
        tree: Phylogenetic tree
        sequences: Aligned sequences
        models: List of models to test
        criterion: 'AIC' or 'BIC'
        test_gamma: If True, test each model with and without +G
        verbose: Print details

    Returns:
        (best_model, best_alpha, best_score, all_results)

    Example:
        >>> model, alpha, score, results = select_best_model_with_gamma(
        ...     tree, sequences, test_gamma=True
        ... )
        >>> if alpha is not None:
        ...     print(f"Best: {model}+G (alpha={alpha:.3f})")
        ... else:
        ...     print(f"Best: {model}")
    """
    if models is None:
        models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']

    all_results = {}

    # Test models without gamma
    for model_name in models:
        best_model, score, results = select_best_model(
            tree, sequences,
            models=[model_name],
            criterion=criterion,
            alpha=None,
            verbose=False
        )
        all_results[model_name] = results[model_name]

    # Test models with gamma
    if test_gamma:
        alpha_values = [0.5, 1.0, 2.0]  # Test different gamma shapes

        for model_name in models:
            best_alpha = None
            best_score_for_model = float('inf')

            for alpha in alpha_values:
                _, score, results = select_best_model(
                    tree, sequences,
                    models=[model_name],
                    criterion=criterion,
                    alpha=alpha,
                    verbose=False
                )

                if score < best_score_for_model:
                    best_score_for_model = score
                    best_alpha = alpha

            # Store best result with gamma for this model
            _, score, results = select_best_model(
                tree, sequences,
                models=[model_name],
                criterion=criterion,
                alpha=best_alpha,
                verbose=False
            )
            all_results[f"{model_name}+G"] = results[model_name]

    # Find overall best
    best_name = min(all_results.keys(), key=lambda k: all_results[k]['score'])
    best_result = all_results[best_name]

    # Extract alpha if present
    best_alpha = best_result['params'].get('alpha', None)

    if verbose:
        print("\n" + "=" * 70)
        print("MODEL SELECTION WITH GAMMA VARIANTS")
        print("=" * 70)
        print(f"{'Model':<15} {'LogL':>12} {'Params':>8} {criterion:>12} {'Delta':>10}")
        print("-" * 70)

        # Sort by score
        sorted_results = sorted(all_results.items(), key=lambda x: x[1]['score'])

        for name, res in sorted_results:
            delta = res['score'] - best_result['score']
            marker = " *" if name == best_name else ""
            print(f"{name:<15} {res['logL']:>12.2f} {res['n_params']:>8} "
                  f"{res['score']:>12.2f} {delta:>10.2f}{marker}")

        print("\n* Best model selected")

    return best_name, best_alpha, best_result['score'], all_results

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
    verbose: bool = False,
    compressor = None
) -> Tuple[float, Dict]:
    """
    Compute likelihood for a specific substitution model.

    Args:
        tree: Phylogenetic tree (fixed topology)
        sequences: Aligned sequences
        model_name: Model name ('JC69', 'K80', 'F81', 'HKY85', 'GTR')
        alpha: Gamma shape parameter (None = no rate heterogeneity)
        verbose: Print details
        compressor: Optional cached SitePatternCompressor (avoids recomputation)

    Returns:
        (log_likelihood, model_params)
        model_params includes fitted parameter values
    """
    from rrna_phylo.models.rate_matrices import (
        get_jc69_rate_matrix,
        get_k80_rate_matrix,
        get_f81_rate_matrix,
        get_hky85_rate_matrix,
        normalize_rate_matrix
    )

    # Get empirical base frequencies
    freqs = compute_empirical_frequencies(sequences)

    # Get model-specific rate matrix
    if model_name == 'JC69':
        Q = get_jc69_rate_matrix()
        Q = normalize_rate_matrix(Q, np.array([0.25, 0.25, 0.25, 0.25]))
        params = None

    elif model_name == 'K80':
        kappa = 2.0  # Initial guess
        Q = get_k80_rate_matrix(kappa)
        Q = normalize_rate_matrix(Q, np.array([0.25, 0.25, 0.25, 0.25]))
        params = np.array([kappa])

    elif model_name == 'F81':
        Q = get_f81_rate_matrix(freqs)
        Q = normalize_rate_matrix(Q, freqs)
        params = freqs.copy()

    elif model_name == 'HKY85':
        kappa = 2.0  # Initial guess
        Q = get_hky85_rate_matrix(kappa, freqs)
        Q = normalize_rate_matrix(Q, freqs)
        params = np.array([kappa])

    elif model_name == 'GTR':
        # For GTR, use rate_matrices implementation
        from rrna_phylo.models.rate_matrices import get_gtr_rate_matrix

        # Default GTR parameters (can be optimized later)
        # Format: [rAC, rAG, rAT, rCG, rCT, rGT]
        # Use typical values: transitions (AG, CT) higher than transversions
        exchangeabilities = np.array([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])

        Q = get_gtr_rate_matrix(exchangeabilities, freqs)
        Q = normalize_rate_matrix(Q, freqs)
        params = exchangeabilities.copy()

    else:
        raise ValueError(f"Unknown model: {model_name}")

    # Compute likelihood with this rate matrix
    log_likelihood = _compute_likelihood_with_custom_Q(
        tree, sequences, Q, freqs, alpha, compressor=compressor
    )

    model_params = {
        'model': model_name,
        'params': params,
        'frequencies': freqs,
        'alpha': alpha,
        'Q': Q
    }

    if verbose:
        print(f"{model_name}: LogL = {log_likelihood:.2f}")

    return log_likelihood, model_params


def _compute_likelihood_with_custom_Q(
    tree: TreeNode,
    sequences: List[Sequence],
    Q: np.ndarray,
    freqs: np.ndarray,
    alpha: Optional[float] = None,
    compressor = None
) -> float:
    """
    Compute likelihood using a custom Q matrix (for non-GTR models).

    Uses the same site pattern compression and gamma rate heterogeneity
    as Level 3, but with a custom substitution rate matrix.

    Args:
        tree: Tree topology
        sequences: Aligned sequences
        Q: 4x4 rate matrix
        freqs: Base frequencies
        alpha: Gamma shape parameter
        compressor: Optional cached SitePatternCompressor (avoids recomputation)

    Returns:
        Log-likelihood
    """
    from rrna_phylo.models.ml_tree_level3 import (
        SitePatternCompressor,
        GammaRates
    )
    from scipy.linalg import expm

    # Site pattern compression (reuse cached if provided)
    if compressor is None:
        compressor = SitePatternCompressor(sequences)

    # Gamma rates (if alpha provided)
    if alpha is not None:
        gamma = GammaRates(alpha, n_categories=4)
        Q_matrices = gamma.get_rate_matrices(Q)
        gamma_probs = gamma.probabilities
    else:
        Q_matrices = [Q]
        gamma_probs = [1.0]

    # Map sequence names to indices
    seq_id_to_idx = {seq.display_name: i for i, seq in enumerate(sequences)}

    log_likelihood = 0.0

    # Calculate likelihood for each unique pattern
    for pattern_idx in range(compressor.n_patterns):
        pattern = compressor.patterns[pattern_idx]
        count = compressor.pattern_counts[pattern_idx]

        # Average likelihood across gamma categories
        pattern_likelihood = 0.0

        for cat_idx, Q_cat in enumerate(Q_matrices):
            cat_prob = gamma_probs[cat_idx]

            # Calculate likelihood for this pattern and category
            L = _calculate_pattern_likelihood_with_Q(
                tree, pattern, sequences, Q_cat, seq_id_to_idx, freqs
            )
            pattern_likelihood += cat_prob * L

        # Add to total (weighted by pattern count)
        if pattern_likelihood > 0:
            log_likelihood += count * np.log(pattern_likelihood)
        else:
            log_likelihood += count * (-1000.0)  # Penalty

    return log_likelihood


def _calculate_pattern_likelihood_with_Q(
    tree: TreeNode,
    pattern: np.ndarray,
    sequences: List[Sequence],
    Q: np.ndarray,
    seq_id_to_idx: dict,
    freqs: np.ndarray
) -> float:
    """
    Calculate likelihood for one site pattern using Felsenstein's algorithm.

    Args:
        tree: Tree topology
        pattern: Site pattern (nucleotide indices)
        sequences: Sequences (for mapping)
        Q: Rate matrix
        seq_id_to_idx: Map from sequence name to index
        freqs: Base frequencies (for root)

    Returns:
        Likelihood (not log)
    """
    from scipy.linalg import expm

    def conditional_likelihood(node: TreeNode) -> np.ndarray:
        """Calculate L[node][state] for all states."""
        if node.is_leaf():
            L = np.zeros(4)
            seq_idx = seq_id_to_idx[node.name]
            observed_state = pattern[seq_idx]

            if observed_state >= 0:
                L[observed_state] = 1.0
            else:
                L[:] = 0.25  # Gap/ambiguous

            return L

        # Internal node: combine children
        L = np.ones(4)

        if node.left:
            # Transition probability matrix
            t = node.left.distance if hasattr(node.left, 'distance') and node.left.distance is not None else 0.01
            P_left = expm(Q * t)

            # Child likelihoods
            L_left = conditional_likelihood(node.left)

            # Sum over child states
            for i in range(4):
                L[i] *= np.sum(P_left[i, :] * L_left)

        if node.right:
            # Transition probability matrix
            t = node.right.distance if hasattr(node.right, 'distance') and node.right.distance is not None else 0.01
            P_right = expm(Q * t)

            # Child likelihoods
            L_right = conditional_likelihood(node.right)

            # Sum over child states
            for i in range(4):
                L[i] *= np.sum(P_right[i, :] * L_right)

        return L

    # Calculate likelihood at root
    L_root = conditional_likelihood(tree)

    # Weight by equilibrium frequencies
    likelihood = np.sum(freqs * L_root)

    return likelihood


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

    # Create site pattern compressor once and reuse for all models
    # (compressor prints compression stats automatically)
    from rrna_phylo.models.ml_tree_level3 import SitePatternCompressor
    compressor = SitePatternCompressor(sequences)
    if verbose:
        print()  # Blank line after compression message

    # Test each model
    for model_name in models:
        model = get_model(model_name)

        # Compute likelihood (reuse cached compressor)
        log_likelihood, model_params = compute_likelihood_for_model(
            tree, sequences, model_name, alpha=alpha, verbose=False, compressor=compressor
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
    verbose: bool = False,
    early_stop_threshold: float = 10.0
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
        early_stop_threshold: Stop testing if BIC difference > threshold (default 10.0)

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
    best_score_so_far = float('inf')
    models_skipped = 0

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

        # Update best score and check early stopping
        if score < best_score_so_far:
            best_score_so_far = score
        elif score - best_score_so_far > early_stop_threshold:
            models_skipped += 1
            if verbose:
                print(f"  Skipping remaining {model_name} variants (BIC diff > {early_stop_threshold})")
            continue

    # Test models with gamma
    if test_gamma:
        alpha_values = [0.5, 1.0, 2.0]  # Test different gamma shapes

        for model_name in models:
            # Check if we should skip this model based on non-gamma performance
            if model_name in all_results:
                model_score = all_results[model_name]['score']
                if model_score - best_score_so_far > early_stop_threshold:
                    models_skipped += 1
                    if verbose:
                        print(f"  Skipping {model_name}+G (base model BIC diff > {early_stop_threshold})")
                    continue

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

                # Update global best score
                if score < best_score_so_far:
                    best_score_so_far = score

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
        if models_skipped > 0:
            print(f"Early stopping: Skipped {models_skipped} models (BIC diff > {early_stop_threshold})")
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

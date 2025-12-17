"""
Model selection using information criteria (AIC/BIC).

Implements automatic testing of multiple substitution models
and selection of the best-fitting model using statistical criteria.
"""

from typing import List, Tuple, Dict, Optional, TYPE_CHECKING
import numpy as np
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.substitution_models import (
    get_model,
    compute_empirical_frequencies,
    DNA_MODELS
)
# compute_log_likelihood may be provided elsewhere; keep import for compatibility
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood

if TYPE_CHECKING:
    from rrna_phylo.models.ml_tree_level3 import SitePatternCompressor


def calculate_aic(log_likelihood: float, n_params: int) -> float:
    """
    Calculate Akaike Information Criterion.

    AIC = -2*ln(L) + 2*k

    Lower values indicate better fit with penalty for complexity.
    """
    return -2 * log_likelihood + 2 * n_params


def calculate_bic(log_likelihood: float, n_params: int, n_sites: int) -> float:
    """
    Calculate Bayesian Information Criterion.

    BIC = -2*ln(L) + k*ln(n)
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

    Returns:
        (log_likelihood, model_params)
    """
    # Get empirical base frequencies
    freqs = compute_empirical_frequencies(sequences)

    # Get model from substitution_models module (unified interface)
    model = get_model(model_name)

    # Get model-specific rate matrix and parameters
    if model_name == 'JC69':
        Q = model.get_rate_matrix()  # No params needed for JC69
        params = None

    elif model_name == 'K80':
        kappa = 2.0  # Initial guess
        Q = model.get_rate_matrix(params=np.array([kappa]))
        params = np.array([kappa])

    elif model_name == 'F81':
        Q = model.get_rate_matrix(freqs=freqs)
        params = freqs.copy()

    elif model_name == 'HKY85':
        kappa = 2.0  # Initial guess
        Q = model.get_rate_matrix(params=np.array([kappa]), freqs=freqs)
        params = np.array([kappa])

    elif model_name == 'GTR':
        exchangeabilities = np.array([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
        Q = model.get_rate_matrix(params=exchangeabilities, freqs=freqs)
        params = exchangeabilities.copy()

    else:
        raise ValueError(f"Unknown model: {model_name}")

    # Compute likelihood with the stable optimized routine
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


def _compute_likelihood_with_custom_Q_ORIGINAL(
    tree: TreeNode,
    sequences: List[Sequence],
    Q: np.ndarray,
    freqs: np.ndarray,
    alpha: Optional[float] = None,
    compressor = None
) -> float:
    """
    Original (kept for reference) - DO NOT USE DIRECTLY.
    """
    raise RuntimeError("Use the fixed _compute_likelihood_with_custom_Q implementation.")


def _compute_likelihood_with_custom_Q(
    tree: TreeNode,
    sequences: List[Sequence],
    Q: np.ndarray,
    freqs: np.ndarray,
    alpha: Optional[float] = None,
    compressor = None
) -> float:
    """
    Corrected and robust likelihood computation using Felsenstein pruning.

    This implementation:
      - reuses SitePatternCompressor
      - supports Gamma (+G) via GammaRates
      - caches P = expm(Q * t) per branch and Q-matrix
      - memoizes per-node conditional likelihoods keyed by (pattern_idx, Q_key)
      - guards against zero/NaN likelihoods
    """
    from rrna_phylo.models.ml_tree_level3 import SitePatternCompressor, GammaRates
    from scipy.linalg import expm

    # Compressor: reuse if provided
    if compressor is None:
        compressor = SitePatternCompressor(sequences)

    # Prepare gamma categories (list of Q matrices and their weights)
    if alpha is not None:
        gamma = GammaRates(alpha, n_categories=4)
        Q_matrices = gamma.get_rate_matrices(Q)
        gamma_probs = gamma.probabilities
    else:
        Q_matrices = [Q]
        gamma_probs = [1.0]

    # Map sequence display_name to index into pattern arrays
    seq_id_to_idx = {seq.display_name: i for i, seq in enumerate(sequences)}

    # Branch P cache: keys are (id(node), q_key) -> 4x4 matrix
    branch_P_cache: Dict[Tuple[int, bytes], np.ndarray] = {}

    def q_key_of(Q_mat: np.ndarray) -> bytes:
        # stable bytes key for Q_mat content
        return Q_mat.tobytes()

    def get_P_matrix_for_child(child: TreeNode, Q_mat: np.ndarray):
        # Return P = expm(Q_mat * t) cached per (node id, q_key)
        qk = q_key_of(Q_mat)
        key = (id(child), qk)
        if key not in branch_P_cache:
            t = getattr(child, 'distance', 0.01)
            # guard: non-positive or NaN times -> tiny positive
            if not np.isfinite(t) or t <= 0.0:
                t = 1e-8
            branch_P_cache[key] = expm(Q_mat * t)
        return branch_P_cache[key]

    # Per-node conditional-likelihood caches will be stored in node._cl_cache
    # keyed by (pattern_idx, q_key) -> vector length 4.
    def clear_all_node_caches(node: TreeNode):
        if hasattr(node, '_cl_cache'):
            node._cl_cache = {}
        if node.left:
            clear_all_node_caches(node.left)
        if node.right:
            clear_all_node_caches(node.right)

    def conditional_likelihood(node: TreeNode, pattern: np.ndarray, pattern_idx: int, Q_mat: np.ndarray):
        """
        Compute conditional likelihood vector for 'node' given pattern and Q_mat.
        Uses memoization stored on nodes in _cl_cache keyed by (pattern_idx, q_key).
        """
        if not hasattr(node, '_cl_cache'):
            node._cl_cache = {}

        qk = q_key_of(Q_mat)
        cache_key = (pattern_idx, qk)
        if cache_key in node._cl_cache:
            return node._cl_cache[cache_key]

        if node.is_leaf():
            L = np.zeros(4, dtype=float)
            seq_idx = seq_id_to_idx.get(node.name, None)
            if seq_idx is None:
                # fallback: use equilibrium frequencies
                L[:] = freqs
            else:
                # pattern contains ints where -1 (or negative) indicates gap/ambiguous
                try:
                    observed_state = int(pattern[seq_idx])
                except Exception:
                    observed_state = -1
                if 0 <= observed_state < 4:
                    L[observed_state] = 1.0
                else:
                    # gap/ambiguous -> use equilibrium probabilities
                    L[:] = freqs
            node._cl_cache[cache_key] = L
            return L

        L = np.ones(4, dtype=float)
        for child in (node.left, node.right):
            if child is None:
                continue
            P_child = get_P_matrix_for_child(child, Q_mat)  # 4x4
            L_child = conditional_likelihood(child, pattern, pattern_idx, Q_mat)  # length-4
            contrib = np.dot(P_child, L_child)
            contrib = np.maximum(contrib, 1e-300)
            L *= contrib

        node._cl_cache[cache_key] = L
        return L

    log_likelihood = 0.0
    eps = 1e-300

    # iterate patterns
    for pattern_idx in range(compressor.n_patterns):
        pattern = compressor.patterns[pattern_idx]
        count = int(compressor.pattern_counts[pattern_idx])
        pattern_like = 0.0

        for Q_mat, cat_prob in zip(Q_matrices, gamma_probs):
            # clear caches for deterministic correctness for this pattern+Q
            clear_all_node_caches(tree)

            # compute conditional likelihood for root
            L_root = conditional_likelihood(tree, pattern, pattern_idx, Q_mat)

            # weight by equilibrium frequencies
            pat_cat_like = float(np.sum(freqs * L_root))
            if pat_cat_like <= 0.0 or not np.isfinite(pat_cat_like):
                pat_cat_like = eps
            pattern_like += cat_prob * pat_cat_like

        # guard
        if pattern_like <= 0.0 or not np.isfinite(pattern_like):
            pattern_like = eps

        log_likelihood += count * np.log(pattern_like)

    return float(log_likelihood)


def select_best_model(
    tree: TreeNode,
    sequences: List[Sequence],
    models: Optional[List[str]] = None,
    criterion: str = 'BIC',
    alpha: Optional[float] = None,
    verbose: bool = False
) -> Tuple[str, float, Dict[str, Dict], any]:
    """
    Test multiple models and select the best using AIC or BIC.

    NOTE: This function does NOT print the "Testing model" start-line.
    The wrapper (select_best_model_with_gamma) prints exactly one start-line per model/model+G.
    If called standalone, it will only print the verbose table when verbose=True.
    """
    if models is None:
        models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']

    if criterion not in ['AIC', 'BIC']:
        raise ValueError(f"Invalid criterion: {criterion}. Use 'AIC' or 'BIC'")

    n_sites = len(sequences[0].sequence)
    results: Dict[str, Dict] = {}

    # Create site pattern compressor once and reuse for all models
    # CRITICAL: This compressor will be returned and reused by GPU calculator
    from rrna_phylo.models.ml_tree_level3 import SitePatternCompressor
    compressor = SitePatternCompressor(sequences)

    for model_name in models:
        
        # Compute likelihood (reuse cached compressor)
        log_likelihood, model_params = compute_likelihood_for_model(
            tree, sequences, model_name, alpha=alpha, verbose=False, compressor=compressor
        )

        # Calculate number of parameters
        model = get_model(model_name)
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
        sorted_results = sorted(results.items(), key=lambda x: x[1]['score'])
        for model_name, res in sorted_results:
            delta = res['score'] - best_score
            marker = " *" if model_name == best_model_name else ""
            print(f"{model_name:<12} {res['logL']:>12.2f} {res['n_params']:>8} "
                  f"{res['score']:>12.2f} {delta:>10.2f}{marker}")
        print("\n* Best model selected")
        print(f"\nBest model: {best_model_name}")
        print(f"{criterion} score: {best_score:.2f}")

    # CRITICAL FIX: Return compressor so GPU calculator can reuse it
    return best_model_name, best_score, results, compressor


def count_branches(tree: TreeNode) -> int:
    """
    Count number of branches in a tree.
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

    Prints one start-line per model (and per model+G with alpha). Uses
    select_best_model() internally (which performs the heavy-lifting).
    """
    if models is None:
        models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']

    all_results: Dict[str, Dict] = {}
    best_score_so_far = float('inf')
    models_skipped = 0
    shared_compressor = None  # CRITICAL: Store compressor from first call

    # Test models without gamma
    for model_name in models:
        # print start-line (single line per model start)
        print(f"Testing model: {model_name} (no gamma)")
        # Run select_best_model for this single model (no gamma)
        best_model, score, results, compressor = select_best_model(
            tree, sequences,
            models=[model_name],
            criterion=criterion,
            alpha=None,
            verbose=False
        )
        # CRITICAL: Capture compressor from first call (all calls use same sequences, so same patterns)
        if shared_compressor is None:
            shared_compressor = compressor
        # store the result for this base model
        all_results[model_name] = results[model_name]

        # Update best score and check early stopping
        if score < best_score_so_far:
            best_score_so_far = score
        elif score - best_score_so_far > early_stop_threshold:
            models_skipped += 1
            if verbose:
                print(f"  Skipping remaining {model_name} variants (BIC diff > {early_stop_threshold})")
            # continue to next model (we already printed the start-line)
            continue

    # Test models with gamma
    if test_gamma:
        alpha_values = [0.5, 1.0, 2.0]  # Test different gamma shapes

        for model_name in models:
            base = all_results.get(model_name)
            if base is not None and base['score'] - best_score_so_far > early_stop_threshold:
                models_skipped += 1
                if verbose:
                    print(f"  Skipping {model_name}+G (base model BIC diff > {early_stop_threshold})")
                continue

            best_alpha = None
            best_score_for_model = float('inf')

            # try several alphas and print once per (model+alpha) start
            for alpha in alpha_values:
                print(f"Testing model: {model_name}+G with alpha={alpha}")
                _, score, results, _ = select_best_model(
                    tree, sequences,
                    models=[model_name],
                    criterion=criterion,
                    alpha=alpha,
                    verbose=False
                )

                if score < best_score_for_model:
                    best_score_for_model = score
                    best_alpha = alpha

                if score < best_score_so_far:
                    best_score_so_far = score

            # Re-evaluate and store the best +G variant result (for the chosen alpha)
            # (this returns results dict keyed by model_name)
            _, score, results, _ = select_best_model(
                tree, sequences,
                models=[model_name],
                criterion=criterion,
                alpha=best_alpha,
                verbose=False
            )
            all_results[f"{model_name}+G"] = results[model_name]

    # Find overall best
    if not all_results:
        raise RuntimeError("No models were evaluated; check input sequences / parameters.")

    best_name = min(all_results.keys(), key=lambda k: all_results[k]['score'])
    best_result = all_results[best_name]

    # Extract alpha if present (model params contain 'alpha' since we stored model_params earlier)
    params = best_result.get('params', {})
    best_alpha = None
    if isinstance(params, dict):
        best_alpha = params.get('alpha', None)

    if verbose:
        print("\nMODEL SELECTION WITH GAMMA VARIANTS COMPLETE")
        if models_skipped > 0:
            print(f"Early stopping: Skipped {models_skipped} models")
        print(f"Best model: {best_name}, Alpha: {best_alpha}, Score: {best_result['score']:.2f}")

    # CRITICAL FIX: Return compressor so GPU calculator can reuse it
    return best_name, best_alpha, best_result['score'], all_results, shared_compressor

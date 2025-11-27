"""
Substitution model rate matrices for phylogenetic likelihood.

Implements proper Q matrices for JC69, K80, F81, HKY85, and GTR models.

References:
- Jukes & Cantor (1969) - "Evolution of protein molecules"
- Kimura (1980) - "A simple method for estimating evolutionary rates"
- Felsenstein (1981) - "Evolutionary trees from DNA sequences"
- Hasegawa et al. (1985) - "Dating of human-ape splitting by molecular clock"
"""

import numpy as np
from typing import Optional


def get_jc69_rate_matrix() -> np.ndarray:
    """
    Jukes-Cantor 1969 model - all substitutions equally likely.

    Assumes:
    - All base frequencies equal (0.25 each)
    - All substitution rates equal

    Returns:
        4x4 rate matrix Q where Q[i,j] is rate from base i to base j
        Order: A(0), C(1), G(2), T(3)
    """
    # All off-diagonal elements = 1.0
    # Diagonal = -3.0 (sum of row)
    Q = np.array([
        [-3.0,  1.0,  1.0,  1.0],
        [ 1.0, -3.0,  1.0,  1.0],
        [ 1.0,  1.0, -3.0,  1.0],
        [ 1.0,  1.0,  1.0, -3.0]
    ], dtype=float)
    return Q


def get_k80_rate_matrix(kappa: float = 2.0) -> np.ndarray:
    """
    Kimura 1980 model - different rates for transitions vs transversions.

    Assumes:
    - Equal base frequencies (0.25 each)
    - Transitions (A<->G, C<->T) have rate kappa
    - Transversions (all others) have rate 1.0

    Args:
        kappa: Transition/transversion rate ratio (typically 2-4 for DNA)

    Returns:
        4x4 rate matrix Q
    """
    # Transitions get rate kappa, transversions get rate 1.0
    Q = np.array([
        [-(kappa + 2.0),  1.0,           kappa,          1.0          ],
        [ 1.0,           -(kappa + 2.0),  1.0,           kappa        ],
        [ kappa,          1.0,           -(kappa + 2.0),  1.0          ],
        [ 1.0,            kappa,          1.0,           -(kappa + 2.0)]
    ], dtype=float)
    return Q


def get_f81_rate_matrix(freqs: np.ndarray) -> np.ndarray:
    """
    Felsenstein 1981 model - equal rates but empirical base frequencies.

    Assumes:
    - Substitution rates proportional to target base frequency
    - All substitutions equally likely (no ti/tv bias)

    Args:
        freqs: Base frequencies [πA, πC, πG, πT]
               Must sum to 1.0

    Returns:
        4x4 rate matrix Q
    """
    assert abs(np.sum(freqs) - 1.0) < 1e-6, "Frequencies must sum to 1.0"

    Q = np.zeros((4, 4), dtype=float)

    # Off-diagonal: rate to base j is πj
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i, j] = freqs[j]

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def get_hky85_rate_matrix(kappa: float, freqs: np.ndarray) -> np.ndarray:
    """
    Hasegawa-Kishino-Yano 1985 model - K80 with empirical frequencies.

    Combines:
    - Unequal base frequencies (like F81)
    - Transition/transversion bias (like K80)

    Args:
        kappa: Transition/transversion ratio
        freqs: Base frequencies [πA, πC, πG, πT]

    Returns:
        4x4 rate matrix Q
    """
    assert abs(np.sum(freqs) - 1.0) < 1e-6, "Frequencies must sum to 1.0"

    Q = np.zeros((4, 4), dtype=float)

    # Define transitions: A(0) <-> G(2), C(1) <-> T(3)
    transitions = [(0, 2), (2, 0), (1, 3), (3, 1)]

    for i in range(4):
        for j in range(4):
            if i != j:
                if (i, j) in transitions:
                    # Transition: rate = kappa * πj
                    Q[i, j] = kappa * freqs[j]
                else:
                    # Transversion: rate = πj
                    Q[i, j] = freqs[j]

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def get_gtr_rate_matrix(
    exchangeabilities: np.ndarray,
    freqs: np.ndarray
) -> np.ndarray:
    """
    General Time Reversible model - most complex reversible model.

    Most flexible model - allows different rates for all 6 possible substitutions.

    Args:
        exchangeabilities: Six exchangeability rates [rAC, rAG, rAT, rCG, rCT, rGT]
                          Typically normalized so rGT = 1.0
        freqs: Base frequencies [πA, πC, πG, πT]

    Returns:
        4x4 rate matrix Q satisfying reversibility: Qij * πi = Qji * πj
    """
    assert len(exchangeabilities) == 6, "Need 6 exchangeability rates"
    assert abs(np.sum(freqs) - 1.0) < 1e-6, "Frequencies must sum to 1.0"

    Q = np.zeros((4, 4), dtype=float)
    r = exchangeabilities

    # Fill matrix using reversibility constraint
    # Upper triangle (and corresponding lower triangle)
    Q[0, 1] = r[0] * freqs[1]  # A->C
    Q[1, 0] = r[0] * freqs[0]  # C->A (reversible)

    Q[0, 2] = r[1] * freqs[2]  # A->G
    Q[2, 0] = r[1] * freqs[0]  # G->A

    Q[0, 3] = r[2] * freqs[3]  # A->T
    Q[3, 0] = r[2] * freqs[0]  # T->A

    Q[1, 2] = r[3] * freqs[2]  # C->G
    Q[2, 1] = r[3] * freqs[1]  # G->C

    Q[1, 3] = r[4] * freqs[3]  # C->T
    Q[3, 1] = r[4] * freqs[1]  # T->C

    Q[2, 3] = r[5] * freqs[3]  # G->T
    Q[3, 2] = r[5] * freqs[2]  # T->G

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def normalize_rate_matrix(Q: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    """
    Normalize rate matrix so mean rate = 1.0 substitution/site/time.

    Standard practice in phylogenetics to make branch lengths interpretable
    as expected number of substitutions per site.

    Args:
        Q: Unnormalized rate matrix
        freqs: Equilibrium frequencies

    Returns:
        Normalized rate matrix Q' where mean rate = 1.0
    """
    # Mean rate = -Σ(πi * Qii)
    # This is the expected rate of substitution
    mean_rate = -np.sum(freqs * np.diag(Q))

    if mean_rate > 0:
        Q_normalized = Q / mean_rate
    else:
        # Should not happen with valid Q matrix
        Q_normalized = Q

    return Q_normalized


def get_model_parameters(model_name: str) -> dict:
    """
    Get metadata about a substitution model.

    Args:
        model_name: One of 'JC69', 'K80', 'F81', 'HKY85', 'GTR'

    Returns:
        Dict with 'n_free_params' and 'description'
    """
    models = {
        'JC69': {
            'n_free_params': 0,
            'description': 'Jukes-Cantor - equal rates, equal frequencies'
        },
        'K80': {
            'n_free_params': 1,  # kappa
            'description': 'Kimura - transition/transversion ratio'
        },
        'F81': {
            'n_free_params': 3,  # 3 frequencies (4th constrained)
            'description': 'Felsenstein - empirical frequencies, equal rates'
        },
        'HKY85': {
            'n_free_params': 4,  # kappa + 3 frequencies
            'description': 'Hasegawa-Kishino-Yano - ti/tv ratio + frequencies'
        },
        'GTR': {
            'n_free_params': 8,  # 5 rates (6th normalized) + 3 frequencies
            'description': 'General Time Reversible - most general model'
        }
    }

    if model_name not in models:
        raise ValueError(f"Unknown model: {model_name}")

    return models[model_name]


# Example usage and testing
if __name__ == "__main__":
    # Test JC69
    print("JC69 Rate Matrix:")
    Q_jc = get_jc69_rate_matrix()
    print(Q_jc)
    print()

    # Test K80
    print("K80 Rate Matrix (kappa=2.0):")
    Q_k80 = get_k80_rate_matrix(kappa=2.0)
    print(Q_k80)
    print()

    # Test F81
    print("F81 Rate Matrix (empirical freqs):")
    freqs = np.array([0.25, 0.25, 0.30, 0.20])  # Example frequencies
    Q_f81 = get_f81_rate_matrix(freqs)
    print(Q_f81)
    print()

    # Test HKY85
    print("HKY85 Rate Matrix (kappa=2.0, empirical freqs):")
    Q_hky = get_hky85_rate_matrix(kappa=2.0, freqs=freqs)
    print(Q_hky)
    print()

    # Test GTR
    print("GTR Rate Matrix:")
    rates = np.array([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])  # Example: high A<->G and C<->T
    Q_gtr = get_gtr_rate_matrix(rates, freqs)
    print(Q_gtr)
    print()

    # Test normalization
    print("Normalized JC69:")
    Q_jc_norm = normalize_rate_matrix(Q_jc, np.array([0.25, 0.25, 0.25, 0.25]))
    print(Q_jc_norm)
    mean_rate = -np.sum(np.array([0.25, 0.25, 0.25, 0.25]) * np.diag(Q_jc_norm))
    print(f"Mean rate after normalization: {mean_rate:.4f} (should be 1.0)")

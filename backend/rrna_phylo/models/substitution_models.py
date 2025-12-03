"""
Substitution model classes for phylogenetic inference.

Implements various DNA substitution models with increasing complexity:
- JC69: Jukes-Cantor (equal rates, equal frequencies)
- K80: Kimura 2-parameter (transition/transversion ratio)
- F81: Felsenstein 1981 (equal rates, empirical frequencies)
- HKY85: Hasegawa-Kishino-Yano (K80 + empirical frequencies)
- GTR: General Time Reversible (full 6-parameter model)
"""

from typing import Tuple, Optional
import numpy as np
from scipy.linalg import expm


class SubstitutionModel:
    """Base class for DNA substitution models."""

    def __init__(self, name: str, n_params: int):
        """
        Initialize substitution model.

        Args:
            name: Model name (e.g., 'JC69', 'GTR')
            n_params: Number of free parameters (for AIC/BIC calculation)
        """
        self.name = name
        self.n_params = n_params
        self.state_order = ['A', 'C', 'G', 'T']

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Get instantaneous rate matrix Q.

        Args:
            params: Model-specific parameters
            freqs: Base frequencies [A, C, G, T]

        Returns:
            4x4 rate matrix Q where Q[i,j] is rate from state i to j

        The Q matrix satisfies:
        - Q[i,j] >= 0 for i != j (non-negative off-diagonal)
        - Sum of each row is 0 (conservation of probability)
        - Scaled so mean rate = 1.0
        """
        raise NotImplementedError

    def get_probability_matrix(
        self,
        branch_length: float,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Get transition probability matrix P(t) = exp(Qt).

        Args:
            branch_length: Time (evolutionary distance)
            params: Model-specific parameters
            freqs: Base frequencies

        Returns:
            4x4 probability matrix P where P[i,j] = Prob(j at time t | i at time 0)
        """
        Q = self.get_rate_matrix(params, freqs)
        P = expm(Q * branch_length)
        return P


class JC69Model(SubstitutionModel):
    """
    Jukes-Cantor 1969 model.

    Simplest model: all substitution rates equal, all frequencies equal (0.25).
    No free parameters.

    Rate matrix:
        Q = [[-3μ,  μ,   μ,   μ  ]
             [ μ,  -3μ,  μ,   μ  ]
             [ μ,   μ,  -3μ,  μ  ]
             [ μ,   μ,   μ,  -3μ ]]

    where μ is scaled so mean rate = 1.0
    """

    def __init__(self):
        super().__init__('JC69', n_params=0)

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """Build JC69 rate matrix (equal rates, equal frequencies)."""
        # All frequencies equal
        pi = np.array([0.25, 0.25, 0.25, 0.25])

        # All off-diagonal rates equal
        Q = np.ones((4, 4)) * 0.25
        np.fill_diagonal(Q, 0)

        # Set diagonal (row sums to 0)
        np.fill_diagonal(Q, -Q.sum(axis=1))

        # Scale so mean rate = 1.0
        # Mean rate = -sum(pi_i * Q[i,i])
        mean_rate = -np.dot(pi, np.diag(Q))
        Q = Q / mean_rate

        return Q


class K80Model(SubstitutionModel):
    """
    Kimura 1980 two-parameter model.

    Distinguishes transitions (A<->G, C<->T) from transversions.
    Equal base frequencies (0.25 each).

    Parameters:
        kappa: transition/transversion rate ratio

    Rate matrix:
        Q = [[-2α-β,   β,    α,    β  ]
             [  β,   -2α-β,  β,    α  ]
             [  α,     β,  -2α-β,  β  ]
             [  β,     α,    β,  -2α-β]]

    where α = transition rate, β = transversion rate
    """

    def __init__(self):
        super().__init__('K80', n_params=1)

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Build K80 rate matrix.

        Args:
            params: [kappa] where kappa = transition/transversion ratio
        """
        if params is None or len(params) < 1:
            kappa = 2.0  # Default transition/transversion ratio
        else:
            kappa = params[0]

        # Equal frequencies
        pi = np.array([0.25, 0.25, 0.25, 0.25])

        # Build rate matrix
        # Transitions: A<->G (indices 0<->2), C<->T (indices 1<->3)
        # Transversions: all other changes
        Q = np.array([
            [0,      1,      kappa,  1     ],  # A -> C(tv), G(ts), T(tv)
            [1,      0,      1,      kappa ],  # C -> A(tv), G(tv), T(ts)
            [kappa,  1,      0,      1     ],  # G -> A(ts), C(tv), T(tv)
            [1,      kappa,  1,      0     ]   # T -> A(tv), C(ts), G(tv)
        ])

        # Multiply by frequencies
        for i in range(4):
            for j in range(4):
                if i != j:
                    Q[i, j] *= pi[j]

        # Set diagonal
        np.fill_diagonal(Q, -Q.sum(axis=1))

        # Scale to mean rate = 1
        mean_rate = -np.dot(pi, np.diag(Q))
        Q = Q / mean_rate

        return Q


class F81Model(SubstitutionModel):
    """
    Felsenstein 1981 model.

    All substitution rates equal, but allows unequal base frequencies.

    Parameters:
        None (uses empirical frequencies from data)

    Rate matrix:
        Q[i,j] = π_j for i != j
        Q[i,i] = -sum_j(Q[i,j])
    """

    def __init__(self):
        super().__init__('F81', n_params=3)  # 3 free frequency parameters (4th constrained)

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Build F81 rate matrix.

        Args:
            freqs: Base frequencies [πA, πC, πG, πT]
        """
        if freqs is None:
            # Default to equal frequencies
            freqs = np.array([0.25, 0.25, 0.25, 0.25])

        # Normalize frequencies
        pi = freqs / freqs.sum()

        # All off-diagonal rates proportional to target frequency
        Q = np.tile(pi, (4, 1))
        np.fill_diagonal(Q, 0)

        # Set diagonal
        np.fill_diagonal(Q, -Q.sum(axis=1))

        # Scale to mean rate = 1
        mean_rate = -np.dot(pi, np.diag(Q))
        Q = Q / mean_rate

        return Q


class HKY85Model(SubstitutionModel):
    """
    Hasegawa-Kishino-Yano 1985 model.

    Combines K80 (transition/transversion distinction) with empirical frequencies.

    Parameters:
        kappa: transition/transversion rate ratio
        freqs: base frequencies [πA, πC, πG, πT]

    Rate matrix:
        Q[i,j] = { kappa * π_j   if i,j are transition
                 { π_j           if i,j are transversion
                 { -sum_k(Q[i,k]) if i == j
    """

    def __init__(self):
        super().__init__('HKY85', n_params=4)  # kappa + 3 frequencies

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Build HKY85 rate matrix.

        Args:
            params: [kappa] transition/transversion ratio
            freqs: Base frequencies [πA, πC, πG, πT]
        """
        if params is None or len(params) < 1:
            kappa = 2.0
        else:
            kappa = params[0]

        if freqs is None:
            freqs = np.array([0.25, 0.25, 0.25, 0.25])

        # Normalize frequencies
        pi = freqs / freqs.sum()

        # Build rate matrix
        # Transitions: A<->G (purines), C<->T (pyrimidines)
        Q = np.zeros((4, 4))

        for i in range(4):
            for j in range(4):
                if i == j:
                    continue

                # Check if transition
                is_transition = (
                    (i == 0 and j == 2) or (i == 2 and j == 0) or  # A <-> G
                    (i == 1 and j == 3) or (i == 3 and j == 1)     # C <-> T
                )

                if is_transition:
                    Q[i, j] = kappa * pi[j]
                else:
                    Q[i, j] = pi[j]

        # Set diagonal
        np.fill_diagonal(Q, -Q.sum(axis=1))

        # Scale to mean rate = 1
        mean_rate = -np.dot(pi, np.diag(Q))
        Q = Q / mean_rate

        return Q


class GTRModel(SubstitutionModel):
    """
    General Time Reversible model.

    Most general reversible model: 6 exchangeability parameters + 4 frequencies.

    Parameters:
        rates: [a, b, c, d, e, f] exchangeability parameters
               A<->C, A<->G, A<->T, C<->G, C<->T, G<->T
        freqs: [πA, πC, πG, πT] base frequencies

    Rate matrix:
        Q[i,j] = r_{ij} * π_j  for i != j
        Q[i,i] = -sum_k(Q[i,k])

    where r_{ij} are exchangeability parameters (symmetric: r_{ij} = r_{ji})
    """

    def __init__(self):
        super().__init__('GTR', n_params=9)  # 6 rates + 3 frequencies

    def get_rate_matrix(
        self,
        params: Optional[np.ndarray] = None,
        freqs: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Build GTR rate matrix.

        Args:
            params: [a, b, c, d, e, f] exchangeability rates
                    A<->C, A<->G, A<->T, C<->G, C<->T, G<->T
            freqs: Base frequencies [πA, πC, πG, πT]
        """
        if params is None or len(params) < 6:
            # Default: all rates equal (reduces to F81)
            params = np.ones(6)

        if freqs is None:
            freqs = np.array([0.25, 0.25, 0.25, 0.25])

        # Normalize frequencies
        pi = freqs / freqs.sum()

        # Unpack exchangeability parameters
        a, b, c, d, e, f = params[:6]

        # Build symmetric exchangeability matrix
        R = np.array([
            [0, a, b, c],  # A -> C, G, T
            [a, 0, d, e],  # C -> A, G, T
            [b, d, 0, f],  # G -> A, C, T
            [c, e, f, 0]   # T -> A, C, G
        ])

        # Build Q = R * diag(π)
        Q = R * pi

        # Set diagonal
        np.fill_diagonal(Q, -Q.sum(axis=1))

        # Scale to mean rate = 1
        mean_rate = -np.dot(pi, np.diag(Q))
        if mean_rate > 0:
            Q = Q / mean_rate

        return Q


# Model registry for easy access
DNA_MODELS = {
    'JC69': JC69Model(),
    'K80': K80Model(),
    'F81': F81Model(),
    'HKY85': HKY85Model(),
    'GTR': GTRModel(),
}


def get_model(name: str) -> SubstitutionModel:
    """
    Get substitution model by name.

    Args:
        name: Model name ('JC69', 'K80', 'F81', 'HKY85', 'GTR')

    Returns:
        SubstitutionModel instance

    Raises:
        ValueError: If model name not recognized
    """
    if name not in DNA_MODELS:
        raise ValueError(
            f"Unknown model: {name}. "
            f"Available models: {list(DNA_MODELS.keys())}"
        )
    return DNA_MODELS[name]


def compute_empirical_frequencies(sequences: list) -> np.ndarray:
    """
    Calculate empirical base frequencies from sequences.

    Args:
        sequences: List of Sequence objects

    Returns:
        Array of frequencies [πA, πC, πG, πT]
    """
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    for seq in sequences:
        sequence = seq.sequence.upper().replace('U', 'T')
        for base in sequence:
            if base in counts:
                counts[base] += 1
                total += 1

    if total == 0:
        return np.array([0.25, 0.25, 0.25, 0.25])

    freqs = np.array([
        counts['A'] / total,
        counts['C'] / total,
        counts['G'] / total,
        counts['T'] / total
    ])

    # VALIDATION: Base frequencies must sum to 1.0
    freq_sum = np.sum(freqs)
    assert np.abs(freq_sum - 1.0) < 1e-10, \
        f"Base frequencies sum ({freq_sum}) != 1.0!"

    return freqs

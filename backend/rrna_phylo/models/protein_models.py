"""
Protein substitution models for phylogenetic inference.

Implements empirical amino acid substitution matrices:
- WAG (Whelan and Goldman, 2001)
- LG (Le and Gascuel, 2008)
- JTT (Jones, Taylor, Thornton, 1992)

Unlike DNA (GTR), protein models use empirically-determined rate matrices
derived from large databases of protein alignments.
"""

import numpy as np
from typing import List
from scipy.linalg import expm
from rrna_phylo.io.fasta_parser import Sequence


class ProteinModel:
    """
    Base class for protein substitution models.

    Protein models are fundamentally different from DNA models:
    - 20 amino acids (vs 4 nucleotides)
    - 20x20 rate matrices (vs 4x4)
    - Empirical rates (not estimated from data)
    - Based on millions of observed substitutions
    """

    def __init__(self):
        """Initialize protein model."""
        # 20 standard amino acids
        self.amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
        self.aa_to_idx = {aa: i for i, aa in enumerate(self.amino_acids)}

        # Amino acid frequencies (will be set by subclass)
        self.aa_freq = np.ones(20) / 20  # Default: uniform

        # Rate matrix Q (will be built by subclass)
        self.Q = None

        # Exchangeability matrix S (symmetric, will be set by subclass)
        self.S = None

    def estimate_frequencies(self, sequences: List[Sequence]):
        """
        Estimate amino acid frequencies from sequences.

        Args:
            sequences: List of aligned protein sequences
        """
        counts = {aa: 0 for aa in self.amino_acids}
        total = 0

        for seq in sequences:
            for aa in seq.sequence.upper():
                if aa in counts:
                    counts[aa] += 1
                    total += 1

        if total > 0:
            self.aa_freq = np.array([counts[aa] / total for aa in self.amino_acids])

        print(f"Estimated amino acid frequencies (top 5):")
        freq_sorted = sorted(zip(self.amino_acids, self.aa_freq), key=lambda x: x[1], reverse=True)
        for aa, freq in freq_sorted[:5]:
            print(f"  {aa}: {freq:.3f}")

    def _build_rate_matrix(self):
        """
        Build rate matrix Q from exchangeability matrix S and frequencies.

        For proteins:
        Q[i,j] = S[i,j] * freq[j]  (for i ≠ j)
        Q[i,i] = -sum(Q[i,:])

        Then normalize so average rate = 1.
        """
        if self.S is None:
            raise ValueError("Exchangeability matrix S not set")

        Q = np.zeros((20, 20))

        # Off-diagonal elements
        for i in range(20):
            for j in range(20):
                if i != j:
                    Q[i, j] = self.S[i, j] * self.aa_freq[j]

        # Diagonal elements
        for i in range(20):
            Q[i, i] = -np.sum(Q[i, :])

        # Normalize: average rate = 1
        mu = -np.sum(np.diagonal(Q) * self.aa_freq)
        if mu > 0:
            Q = Q / mu

        self.Q = Q

    def probability_matrix(self, branch_length: float) -> np.ndarray:
        """
        Calculate probability matrix P(t) = e^(Qt).

        Args:
            branch_length: Evolutionary time

        Returns:
            20x20 probability matrix
        """
        if self.Q is None:
            raise ValueError("Must build rate matrix first")

        # P(t) = e^(Q*t)
        P = expm(self.Q * branch_length)

        return P


class WAGModel(ProteinModel):
    """
    WAG (Whelan and Goldman, 2001) protein substitution model.

    Derived from globular proteins in PFAM database.
    General-purpose model, works well for most proteins.

    Reference:
    Whelan, S. and Goldman, N. (2001)
    "A general empirical model of protein evolution derived from multiple
    protein families using a maximum-likelihood approach"
    Mol Biol Evol 18:691-699
    """

    def __init__(self):
        """Initialize WAG model with empirical parameters."""
        super().__init__()

        # WAG amino acid frequencies (empirical)
        self.aa_freq = np.array([
            0.0866279,  # A
            0.043972,   # C
            0.0390894,  # D
            0.0570451,  # E
            0.0193078,  # F
            0.0367281,  # G
            0.0580589,  # H
            0.0832518,  # I
            0.0244313,  # K
            0.048466,   # L
            0.086209,   # M
            0.0395639,  # N
            0.0428693,  # P
            0.0514261,  # Q
            0.0620634,  # R
            0.0698339,  # S
            0.0582827,  # T
            0.0134815,  # V
            0.0250624,  # W
            0.0315382,  # Y
        ])

        # WAG exchangeability matrix S (symmetric 20x20)
        # This is a simplified version - full matrix has 190 parameters
        # For now, use JTT-like structure (will implement full WAG later)
        self._set_simplified_exchangeabilities()

        # Build rate matrix
        self._build_rate_matrix()

    def _set_simplified_exchangeabilities(self):
        """
        Set simplified exchangeability matrix.

        Full WAG has 190 parameters. This is a placeholder
        that uses a simplified structure for now.

        TODO: Implement full WAG matrix
        """
        # Simplified: all exchanges have rate 1.0
        # (Real WAG has empirically-determined rates)
        self.S = np.ones((20, 20))
        np.fill_diagonal(self.S, 0)

        print("Note: Using simplified exchangeability matrix")
        print("TODO: Implement full WAG empirical matrix (190 parameters)")


class LGModel(ProteinModel):
    """
    LG (Le and Gascuel, 2008) protein substitution model.

    Derived from recent phylogenomic databases.
    Better than WAG for closely-related proteins.

    Reference:
    Le, S.Q. and Gascuel, O. (2008)
    "An improved general amino acid replacement matrix"
    Mol Biol Evol 25:1307-1320
    """

    def __init__(self):
        """Initialize LG model with empirical parameters."""
        super().__init__()

        # LG amino acid frequencies (empirical)
        self.aa_freq = np.array([
            0.079066,   # A
            0.055941,   # C
            0.041977,   # D
            0.053052,   # E
            0.012937,   # F
            0.073152,   # G
            0.022944,   # H
            0.062029,   # I
            0.064718,   # K
            0.0914  ,   # L
            0.022815,   # M
            0.042645,   # N
            0.04074 ,   # P
            0.046089,   # Q
            0.051769,   # R
            0.071586,   # S
            0.05877 ,   # T
            0.066039,   # V
            0.014157,   # W
            0.031726,   # Y
        ])

        # LG exchangeability matrix (simplified for now)
        self._set_simplified_exchangeabilities()

        # Build rate matrix
        self._build_rate_matrix()

    def _set_simplified_exchangeabilities(self):
        """Set simplified exchangeability matrix (placeholder)."""
        self.S = np.ones((20, 20))
        np.fill_diagonal(self.S, 0)

        print("Note: Using simplified exchangeability matrix")
        print("TODO: Implement full LG empirical matrix (190 parameters)")


class JTTModel(ProteinModel):
    """
    JTT (Jones, Taylor, Thornton, 1992) protein substitution model.

    Classic model, derived from analysis of 16,130 substitutions
    in 50 families of protein sequences.

    Reference:
    Jones, D.T., Taylor, W.R., and Thornton, J.M. (1992)
    "The rapid generation of mutation data matrices from protein sequences"
    Comput Appl Biosci 8:275-282
    """

    def __init__(self):
        """Initialize JTT model with empirical parameters."""
        super().__init__()

        # JTT amino acid frequencies (empirical)
        self.aa_freq = np.array([
            0.077,   # A
            0.019,   # C
            0.052,   # D
            0.063,   # E
            0.041,   # F
            0.074,   # G
            0.022,   # H
            0.052,   # I
            0.057,   # K
            0.091,   # L
            0.024,   # M
            0.044,   # N
            0.051,   # P
            0.041,   # Q
            0.051,   # R
            0.069,   # S
            0.059,   # T
            0.066,   # V
            0.013,   # W
            0.033,   # Y
        ])

        # JTT exchangeability matrix (simplified for now)
        self._set_simplified_exchangeabilities()

        # Build rate matrix
        self._build_rate_matrix()

    def _set_simplified_exchangeabilities(self):
        """Set simplified exchangeability matrix (placeholder)."""
        self.S = np.ones((20, 20))
        np.fill_diagonal(self.S, 0)

        print("Note: Using simplified exchangeability matrix")
        print("TODO: Implement full JTT empirical matrix (190 parameters)")


def get_protein_model(model_name: str = "WAG") -> ProteinModel:
    """
    Get a protein substitution model by name.

    Args:
        model_name: Name of model ("WAG", "LG", or "JTT")

    Returns:
        ProteinModel instance

    Raises:
        ValueError: If model name not recognized
    """
    models = {
        "WAG": WAGModel,
        "LG": LGModel,
        "JTT": JTTModel,
    }

    if model_name.upper() not in models:
        raise ValueError(f"Unknown protein model: {model_name}. "
                        f"Available: {', '.join(models.keys())}")

    return models[model_name.upper()]()


# Example usage
if __name__ == "__main__":
    print("=" * 60)
    print("PROTEIN SUBSTITUTION MODELS")
    print("=" * 60)

    # Test each model
    for model_name in ["WAG", "LG", "JTT"]:
        print(f"\n{model_name} Model:")
        print("-" * 60)

        model = get_protein_model(model_name)

        print(f"Amino acids: {len(model.amino_acids)}")
        print(f"Matrix size: {model.Q.shape}")
        print(f"Top 3 amino acid frequencies:")
        top3 = sorted(zip(model.amino_acids, model.aa_freq),
                     key=lambda x: x[1], reverse=True)[:3]
        for aa, freq in top3:
            print(f"  {aa}: {freq:.4f}")

        # Test probability matrix
        P = model.probability_matrix(0.1)
        print(f"P(t=0.1) shape: {P.shape}")
        print(f"P row sums (should be ~1.0): {P.sum(axis=1)[:3]}")

    print("\n" + "=" * 60)
    print("All models initialized successfully!")
    print("\nNote: Currently using simplified exchangeability matrices")
    print("TODO: Implement full empirical matrices (190 parameters each)")
    print("=" * 60)

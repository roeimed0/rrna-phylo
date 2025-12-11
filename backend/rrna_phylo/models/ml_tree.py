"""Maximum Likelihood phylogenetic tree inference with GTR model."""

import numpy as np
from typing import List, Tuple, Dict
from scipy.linalg import expm
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode


class GTRModel:
    """General Time Reversible substitution model for DNA/RNA."""

    def __init__(self):
        """Initialize GTR model with default parameters."""
        # Nucleotide order: A, C, G, T (U mapped to T for RNA)
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.nuc_to_idx = {
            'A': 0, 'C': 1, 'G': 2, 'T': 3,
            'U': 3  # RNA: U -> T mapping
        }

        # Base frequencies (will estimate from data)
        self.base_freq = np.array([0.25, 0.25, 0.25, 0.25])

        # GTR rate parameters (will estimate from data)
        # Order: AC, AG, AT, CG, CT, GT
        self.rates = np.array([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
        # Note: AG and CT are higher (transitions vs transversions)

        # Rate matrix Q (will be computed)
        self.Q = None

    def estimate_parameters(self, sequences: List[Sequence]):
        """Estimate GTR parameters from aligned sequences."""
        # Count nucleotide frequencies (RNA U -> T)
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total = 0

        for seq in sequences:
            for base in seq.sequence.upper():
                # Map RNA U to T
                if base == 'U':
                    base = 'T'
                if base in counts:
                    counts[base] += 1
                    total += 1

        # Convert to frequencies
        if total > 0:
            self.base_freq = np.array([
                counts['A'] / total,
                counts['C'] / total,
                counts['G'] / total,
                counts['T'] / total
            ])


        # Build rate matrix
        self._build_rate_matrix()

    def _build_rate_matrix(self):
        """Build the GTR rate matrix Q."""
        Q = np.zeros((4, 4))

        # Rate parameter indices
        # AC=0, AG=1, AT=2, CG=3, CT=4, GT=5
        rate_map = {
            (0, 1): 0,  # A -> C
            (0, 2): 1,  # A -> G
            (0, 3): 2,  # A -> T
            (1, 2): 3,  # C -> G
            (1, 3): 4,  # C -> T
            (2, 3): 5,  # G -> T
        }

        # Fill off-diagonal elements
        for i in range(4):
            for j in range(4):
                if i != j:
                    key = (min(i, j), max(i, j))
                    rate_idx = rate_map[key]
                    Q[i, j] = self.rates[rate_idx] * self.base_freq[j]

        # Fill diagonal (negative sum of row)
        for i in range(4):
            Q[i, i] = -np.sum(Q[i, :])

        # Normalize so average rate = 1
        mu = -np.sum(np.diagonal(Q) * self.base_freq)
        Q = Q / mu

        self.Q = Q

    def probability_matrix(self, branch_length: float) -> np.ndarray:
        """Calculate probability matrix P(t) = e^(Qt)."""
        if self.Q is None:
            raise ValueError("Must estimate parameters first")

        # P(t) = e^(Q*t)
        P = expm(self.Q * branch_length)

        return P


class MaximumLikelihoodTree:
    """Maximum Likelihood phylogenetic tree inference."""

    def __init__(self):
        """Initialize ML tree builder."""
        self.model = GTRModel()
        self.sequences = None
        self.alignment = None

    def build_tree(self, sequences: List[Sequence]) -> TreeNode:
        """Build Maximum Likelihood tree from aligned sequences."""
        self.sequences = sequences
        self._prepare_alignment()

        self.model.estimate_parameters(sequences)

        from rrna_phylo.methods.bionj import build_bionj_tree
        from rrna_phylo.distance.distance import calculate_distance_matrix

        dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
        initial_tree = build_bionj_tree(dist_matrix, ids)

        likelihood = self._calculate_likelihood(initial_tree)

        return initial_tree

    def _prepare_alignment(self):
        """Convert sequences to numeric alignment matrix."""
        n_seq = len(self.sequences)
        seq_len = self.sequences[0].aligned_length

        alignment = np.zeros((n_seq, seq_len), dtype=int)

        for i, seq in enumerate(self.sequences):
            for j, base in enumerate(seq.sequence.upper()):
                if base in self.model.nuc_to_idx:
                    alignment[i, j] = self.model.nuc_to_idx[base]
                else:
                    alignment[i, j] = -1  # Gap or unknown

        self.alignment = alignment

    def _calculate_likelihood(self, tree: TreeNode) -> float:
        """Calculate likelihood of tree given alignment (placeholder)."""
        return -1000.0  # Placeholder


def build_ml_tree(sequences: List[Sequence]) -> TreeNode:
    """Build Maximum Likelihood tree."""
    builder = MaximumLikelihoodTree()
    return builder.build_tree(sequences)

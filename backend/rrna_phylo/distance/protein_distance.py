"""
Protein distance calculation with correction models.

For proteins, we can't use Jukes-Cantor (which assumes 4 equal-frequency nucleotides).
Instead, we use protein-specific distance corrections.
"""

import math
from typing import List, Tuple
import numpy as np
from rrna_phylo.io.fasta_parser import Sequence


class ProteinDistanceCalculator:
    """Calculate evolutionary distances between protein sequences."""

    # Standard amino acids
    AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

    def __init__(self, model: str = "poisson"):
        """Initialize protein distance calculator."""
        self.model = model.lower()

    def pairwise_distance(self, seq1: Sequence, seq2: Sequence) -> float:
        """Calculate evolutionary distance between two aligned sequences."""
        if seq1.aligned_length != seq2.aligned_length:
            raise ValueError(
                f"Sequences must be same length: {seq1.id}={seq1.aligned_length}, "
                f"{seq2.id}={seq2.aligned_length}"
            )

        # Count differences (excluding gaps)
        differences = 0
        valid_positions = 0

        for i in range(seq1.aligned_length):
            aa1 = seq1.sequence[i].upper()
            aa2 = seq2.sequence[i].upper()

            # Skip positions with gaps
            if aa1 in ('-', '.', 'X') or aa2 in ('-', '.', 'X'):
                continue

            # Only count valid amino acids
            if aa1 not in self.AMINO_ACIDS or aa2 not in self.AMINO_ACIDS:
                continue

            valid_positions += 1

            # Count if different
            if aa1 != aa2:
                differences += 1

        if valid_positions == 0:
            raise ValueError(f"No valid positions to compare between {seq1.id} and {seq2.id}")

        # Calculate p-distance (proportion of differences)
        p = differences / valid_positions

        # Apply correction based on model
        if self.model == "p-distance":
            return p
        elif self.model == "poisson":
            return self._poisson_correction(p)
        elif self.model == "kimura":
            return self._kimura_correction(p)
        else:
            raise ValueError(f"Unknown model: {self.model}")

    def _poisson_correction(self, p: float) -> float:
        """Apply Poisson correction: d = -ln(1 - p)."""
        if p >= 1.0:
            raise ValueError(
                f"Sequences too divergent (p={p:.3f}). "
                f"Poisson correction requires p < 1.0"
            )

        if p == 0:
            return 0.0

        # d = -ln(1 - p)
        d = -math.log(1 - p)

        return d

    def _kimura_correction(self, p: float) -> float:
        """Apply Kimura correction: d = -ln(1 - p - p²/5)."""
        # Check validity
        threshold = 0.85  # Kimura breaks down above ~85% divergence

        if p >= threshold:
            raise ValueError(
                f"Sequences too divergent (p={p:.3f}). "
                f"Kimura correction requires p < {threshold}"
            )

        if p == 0:
            return 0.0

        # d = -ln(1 - p - p²/5)
        inner = 1 - p - (p * p / 5.0)

        if inner <= 0:
            raise ValueError(
                f"Cannot calculate Kimura distance: p={p:.3f} too high"
            )

        d = -math.log(inner)

        return d

    def calculate_distance_matrix(
        self,
        sequences: List[Sequence]
    ) -> Tuple[np.ndarray, List[str]]:
        """Calculate pairwise distance matrix for all sequences."""
        n = len(sequences)
        dist_matrix = np.zeros((n, n))
        ids = [seq.display_name for seq in sequences]

        for i in range(n):
            for j in range(i + 1, n):
                try:
                    dist = self.pairwise_distance(sequences[i], sequences[j])
                    dist_matrix[i, j] = dist
                    dist_matrix[j, i] = dist
                except ValueError as e:
                    # Sequences too divergent - use maximum distance
                    print(f"Warning: {sequences[i].id} vs {sequences[j].id}: {e}")
                    dist_matrix[i, j] = 10.0  # Large but finite
                    dist_matrix[j, i] = 10.0

        return dist_matrix, ids


def calculate_protein_distance_matrix(
    sequences: List[Sequence],
    model: str = "poisson"
) -> Tuple[np.ndarray, List[str]]:
    """Calculate protein distance matrix."""
    calculator = ProteinDistanceCalculator(model=model)
    return calculator.calculate_distance_matrix(sequences)



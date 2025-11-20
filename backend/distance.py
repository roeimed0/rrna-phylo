"""
Calculate pairwise evolutionary distances from aligned sequences.

Uses Jukes-Cantor correction - the gold standard for phylogenetics.
"""

import math
from typing import List, Tuple
import numpy as np
from fasta_parser import Sequence


class DistanceCalculator:
    """Calculate evolutionary distances between sequences."""

    def __init__(self, model: str = "jukes-cantor"):
        """
        Initialize distance calculator.

        Args:
            model: Distance model to use ("jukes-cantor" or "p-distance")
        """
        self.model = model

    def pairwise_distance(self, seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate evolutionary distance between two aligned sequences.

        Args:
            seq1: First aligned sequence
            seq2: Second aligned sequence

        Returns:
            Evolutionary distance (0.0 to infinity)

        Raises:
            ValueError: If sequences have different lengths
        """
        if seq1.aligned_length != seq2.aligned_length:
            raise ValueError(
                f"Sequences must be same length: {seq1.id}={seq1.aligned_length}, "
                f"{seq2.id}={seq2.aligned_length}"
            )

        # Count differences (excluding gaps)
        differences = 0
        valid_positions = 0

        for i in range(seq1.aligned_length):
            base1 = seq1.sequence[i].upper()
            base2 = seq2.sequence[i].upper()

            # Skip positions with gaps
            if base1 in ('-', '.') or base2 in ('-', '.'):
                continue

            valid_positions += 1

            # Count if different
            if base1 != base2:
                differences += 1

        if valid_positions == 0:
            raise ValueError(f"No valid positions to compare between {seq1.id} and {seq2.id}")

        # Calculate p-distance (proportion of differences)
        p = differences / valid_positions

        # Apply correction based on model
        if self.model == "p-distance":
            return p
        elif self.model == "jukes-cantor":
            return self._jukes_cantor_correction(p)
        else:
            raise ValueError(f"Unknown model: {self.model}")

    def _jukes_cantor_correction(self, p: float) -> float:
        """
        Apply Jukes-Cantor correction to p-distance.

        Formula: d = -3/4 * ln(1 - 4p/3)

        This corrects for multiple substitutions at the same site.

        Args:
            p: Observed proportion of differences (0.0 to 1.0)

        Returns:
            Corrected evolutionary distance

        Raises:
            ValueError: If p >= 0.75 (sequences too divergent)
        """
        # Check if correction is possible
        if p >= 0.75:
            raise ValueError(
                f"Sequences too divergent (p={p:.3f}). "
                "Jukes-Cantor correction requires p < 0.75"
            )

        # Jukes-Cantor formula
        d = -0.75 * math.log(1 - (4 * p / 3))

        return d

    def distance_matrix(self, sequences: List[Sequence]) -> Tuple[np.ndarray, List[str]]:
        """
        Calculate pairwise distance matrix for all sequences.

        Args:
            sequences: List of aligned Sequence objects

        Returns:
            Tuple of (distance_matrix, sequence_ids)
            - distance_matrix: NxN numpy array of distances
            - sequence_ids: List of sequence IDs in matrix order

        Raises:
            ValueError: If sequences are not aligned or < 2 sequences
        """
        n = len(sequences)

        if n < 2:
            raise ValueError("Need at least 2 sequences for distance matrix")

        # Check all sequences are same length
        lengths = [seq.aligned_length for seq in sequences]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"Sequences must be aligned (same length). "
                f"Found lengths: {set(lengths)}"
            )

        # Initialize matrix
        matrix = np.zeros((n, n))
        ids = [seq.id for seq in sequences]

        # Calculate pairwise distances
        for i in range(n):
            for j in range(i + 1, n):
                try:
                    dist = self.pairwise_distance(sequences[i], sequences[j])
                    matrix[i][j] = dist
                    matrix[j][i] = dist  # Symmetric
                except ValueError as e:
                    # If sequences too divergent, use maximum distance
                    print(f"Warning: {sequences[i].id} vs {sequences[j].id}: {e}")
                    matrix[i][j] = 10.0  # Large distance
                    matrix[j][i] = 10.0

        return matrix, ids

    def print_matrix(self, matrix: np.ndarray, ids: List[str]):
        """
        Pretty-print distance matrix.

        Args:
            matrix: Distance matrix
            ids: Sequence IDs
        """
        n = len(ids)

        # Print header
        print("\nDistance Matrix:")
        print(f"{'':12}", end='')
        for id in ids:
            print(f"{id:12}", end='')
        print()

        # Print rows
        for i in range(n):
            print(f"{ids[i]:12}", end='')
            for j in range(n):
                print(f"{matrix[i][j]:12.6f}", end='')
            print()


def calculate_distance_matrix(
    sequences: List[Sequence],
    model: str = "jukes-cantor"
) -> Tuple[np.ndarray, List[str]]:
    """
    Convenience function to calculate distance matrix.

    Args:
        sequences: List of aligned sequences
        model: Distance model ("jukes-cantor" or "p-distance")

    Returns:
        Tuple of (distance_matrix, sequence_ids)
    """
    calculator = DistanceCalculator(model=model)
    return calculator.distance_matrix(sequences)

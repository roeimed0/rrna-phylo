"""
Calculate pairwise evolutionary distances from aligned sequences.

Uses Jukes-Cantor correction - the gold standard for phylogenetics.
"""

import math
from typing import List, Tuple
import numpy as np
from rrna_phylo.io.fasta_parser import Sequence


class DistanceCalculator:
    """Calculate evolutionary distances between sequences."""

    def __init__(self, model: str = "jukes-cantor"):
        self.model = model

    def pairwise_distance(self, seq1: Sequence, seq2: Sequence) -> float:
        if seq1.aligned_length != seq2.aligned_length:
            raise ValueError(
                f"Sequences must be same length: {seq1.id}={seq1.aligned_length}, "
                f"{seq2.id}={seq2.aligned_length}"
            )

        differences = 0
        valid_positions = 0

        for i in range(seq1.aligned_length):
            base1 = seq1.sequence[i].upper()
            base2 = seq2.sequence[i].upper()

            if base1 in ('-', '.') or base2 in ('-', '.'):
                continue

            valid_positions += 1
            if base1 != base2:
                differences += 1

        if valid_positions == 0:
            raise ValueError(f"No valid positions to compare between {seq1.id} and {seq2.id}")

        p = differences / valid_positions

        if self.model == "p-distance":
            return p
        elif self.model == "jukes-cantor":
            return self._jukes_cantor_correction(p)
        else:
            raise ValueError(f"Unknown model: {self.model}")

    def _jukes_cantor_correction(self, p: float) -> float:
        if p >= 0.75:
            raise ValueError(
                f"Sequences too divergent (p={p:.3f}). "
                "Jukes-Cantor correction requires p < 0.75"
            )

        # Clamp for floating point safety
        value = 1 - (4 * p / 3)
        if value <= 0:
            value = 1e-12

        return -0.75 * math.log(value)

    def distance_matrix(self, sequences: List[Sequence]) -> Tuple[np.ndarray, List[str]]:
        n = len(sequences)
        if n < 2:
            raise ValueError("Need at least 2 sequences for distance matrix")

        lengths = [seq.aligned_length for seq in sequences]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"Sequences must be aligned (same length). "
                f"Found lengths: {set(lengths)}"
            )

        matrix = np.zeros((n, n))
        ids = [seq.display_name for seq in sequences]

        warnings_count = 0

        for i in range(n):
            for j in range(i + 1, n):
                try:
                    dist = self.pairwise_distance(sequences[i], sequences[j])
                except ValueError:
                    warnings_count += 1
                    dist = 10.0

                matrix[i][j] = matrix[j][i] = dist

        if warnings_count > 0:
            total = (n * (n - 1)) // 2
            pct = 100 * warnings_count / total
            print(f"\nDistance calculation: {warnings_count} of {total} pairs ({pct:.1f}%) too divergent")

        return matrix, ids

    def print_matrix(self, matrix: np.ndarray, ids: List[str]):
        n = len(ids)

        print("\nDistance Matrix:")
        print(f"{'':12}", end='')
        for id in ids:
            print(f"{id:12}", end='')
        print()

        for i in range(n):
            print(f"{ids[i]:12}", end='')
            for j in range(n):
                print(f"{matrix[i][j]:12.6f}", end='')
            print()


def calculate_distance_matrix(
    sequences: List[Sequence],
    model: str = "jukes-cantor"
) -> Tuple[np.ndarray, List[str]]:
    calculator = DistanceCalculator(model=model)
    return calculator.distance_matrix(sequences)

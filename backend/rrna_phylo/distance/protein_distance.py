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
        """
        Initialize protein distance calculator.

        Args:
            model: Distance model to use
                - "p-distance": Simple proportion of differences
                - "poisson": Poisson correction (-ln(1-p))
                - "kimura": Kimura protein distance
        """
        self.model = model.lower()

    def pairwise_distance(self, seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate evolutionary distance between two aligned protein sequences.

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
        """
        Apply Poisson correction to p-distance.

        Poisson model assumes:
        - All amino acid substitutions equally likely
        - Single substitution per site
        - Rate is constant

        Formula: d = -ln(1 - p)

        Args:
            p: Proportion of differences

        Returns:
            Corrected distance

        Raises:
            ValueError: If p >= 1.0 (sequences too divergent)
        """
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
        """
        Apply Kimura protein distance correction.

        Kimura's formula accounts for:
        - Multiple substitutions at same site
        - Assumes equal substitution rates (like Poisson)

        Formula: d = -ln(1 - p - p²/5)

        This is slightly more conservative than Poisson for high divergence.

        Args:
            p: Proportion of differences

        Returns:
            Corrected distance

        Raises:
            ValueError: If sequences too divergent
        """
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
        """
        Calculate pairwise distance matrix for all sequences.

        Args:
            sequences: List of aligned protein sequences

        Returns:
            (distance_matrix, sequence_ids)
            - distance_matrix: NxN numpy array
            - sequence_ids: List of sequence IDs in same order

        Example:
            calculator = ProteinDistanceCalculator(model="poisson")
            dist_matrix, ids = calculator.calculate_distance_matrix(sequences)
        """
        n = len(sequences)
        dist_matrix = np.zeros((n, n))
        ids = [seq.id for seq in sequences]

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
    """
    Convenient function to calculate protein distance matrix.

    Args:
        sequences: Aligned protein sequences
        model: Distance model ("poisson", "kimura", or "p-distance")

    Returns:
        (distance_matrix, sequence_ids)

    Example:
        dist_matrix, ids = calculate_protein_distance_matrix(sequences, model="poisson")
    """
    calculator = ProteinDistanceCalculator(model=model)
    return calculator.calculate_distance_matrix(sequences)


# Example usage and testing
if __name__ == "__main__":
    print("=" * 60)
    print("PROTEIN DISTANCE CALCULATION")
    print("=" * 60)

    # Create test protein sequences (cytochrome c-like)
    # These are short for testing - real proteins are longer
    sequences = [
        Sequence("human", "Human cytochrome c", "GDVEKGKKIFIMKCSQCHTVEK"),
        Sequence("chimp", "Chimp cytochrome c", "GDVEKGKKIFIMKCSQCHTVEK"),  # Identical
        Sequence("mouse", "Mouse cytochrome c", "GDVEKGKKIFVMKCSQCHTVEK"),  # 1 diff (I->V)
        Sequence("chicken", "Chicken cytochrome c", "GDVAKGKKIFIMKCSQCHTVAK"),  # 2 diffs
    ]

    print("\nTest sequences:")
    for seq in sequences:
        print(f"  {seq.id:10} {seq.sequence}")

    # Test all three models
    for model in ["p-distance", "poisson", "kimura"]:
        print(f"\n{model.upper()} Model:")
        print("-" * 60)

        calculator = ProteinDistanceCalculator(model=model)

        # Human vs Chimp (identical)
        d = calculator.pairwise_distance(sequences[0], sequences[1])
        print(f"  Human-Chimp:   {d:.6f} (identical, should be 0)")

        # Human vs Mouse (1 substitution in 22 positions = 4.5%)
        d = calculator.pairwise_distance(sequences[0], sequences[2])
        print(f"  Human-Mouse:   {d:.6f} (1 diff / 22 sites)")

        # Human vs Chicken (2 substitutions = 9.1%)
        d = calculator.pairwise_distance(sequences[0], sequences[3])
        print(f"  Human-Chicken: {d:.6f} (2 diffs / 22 sites)")

    # Calculate full distance matrix
    print("\n\nDistance Matrix (Poisson correction):")
    print("-" * 60)

    dist_matrix, ids = calculate_protein_distance_matrix(sequences, model="poisson")

    print(f"{'':12}", end='')
    for id in ids:
        print(f"{id:12}", end='')
    print()

    for i, id in enumerate(ids):
        print(f"{id:12}", end='')
        for j in range(len(ids)):
            print(f"{dist_matrix[i][j]:12.6f}", end='')
        print()

    print("\n" + "=" * 60)
    print("PROTEIN DISTANCE MODELS WORKING!")
    print("=" * 60)
    print("\nComparison of models:")
    print("  - p-distance: Simple, no correction (underestimates)")
    print("  - Poisson: Standard correction, assumes equal rates")
    print("  - Kimura: More conservative, better for moderate divergence")
    print("\nRecommendation: Use Poisson for general protein phylogenetics")
    print("=" * 60)

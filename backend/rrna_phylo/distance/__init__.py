"""
Distance calculation modules for phylogenetic analysis.

This package contains:
- distance: DNA/RNA distance calculations (Jukes-Cantor, Kimura)
- protein_distance: Protein distance calculations (Poisson, Kimura)
"""

from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.distance.protein_distance import (
    ProteinDistanceCalculator,
    calculate_protein_distance_matrix
)

__all__ = [
    "calculate_distance_matrix",
    "ProteinDistanceCalculator",
    "calculate_protein_distance_matrix",
]

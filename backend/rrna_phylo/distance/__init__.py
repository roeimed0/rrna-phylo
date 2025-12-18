"""
Distance calculation modules for phylogenetic analysis.

This package contains:
- distance: DNA/RNA distance calculations (Jukes-Cantor, Kimura)
"""

from rrna_phylo.distance.distance import calculate_distance_matrix

__all__ = [
    "calculate_distance_matrix",
]

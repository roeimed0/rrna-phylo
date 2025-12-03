"""
Analysis package for phylogenetic analysis tools.
"""

from .bootstrap import (
    calculate_bootstrap_support,
    bootstrap_tree,
    bootstrap_consensus
)
from .dataset_analyzer import (
    analyze_distance_matrix,
    recommend_methods,
    DatasetAnalysis
)

__all__ = [
    'calculate_bootstrap_support',
    'bootstrap_tree',
    'bootstrap_consensus',
    'analyze_distance_matrix',
    'recommend_methods',
    'DatasetAnalysis',
]

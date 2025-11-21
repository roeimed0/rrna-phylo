"""
Phylogenetic models for maximum likelihood inference.

This package contains:
- ml_tree: Basic GTR model implementation
- ml_tree_level2: GTR model with rate heterogeneity
- ml_tree_level3: GTR+Gamma model with site pattern compression (recommended)
- protein_models: Protein substitution models (WAG, LG, JTT)
"""

from rrna_phylo.models.ml_tree_level3 import build_ml_tree_level3, GammaRates
from rrna_phylo.models.protein_models import ProteinModel, get_protein_model

__all__ = [
    "build_ml_tree_level3",
    "GammaRates",
    "ProteinModel",
    "get_protein_model",
]

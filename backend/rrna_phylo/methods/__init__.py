"""
Phylogenetic tree building methods.

This package contains:
- upgma: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- bionj: BioNJ (variance-weighted neighbor-joining)
- protein_ml: Maximum likelihood for protein sequences
"""

from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.protein_ml import build_protein_ml_tree

__all__ = [
    "build_upgma_tree",
    "build_bionj_tree",
    "build_protein_ml_tree",
]

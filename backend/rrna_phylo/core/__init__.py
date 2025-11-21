"""
Core functionality for phylogenetic tree building.

This package contains:
- tree: TreeNode class for representing phylogenetic trees
- sequence_type: Automatic detection of DNA/RNA/Protein sequences
- builder: High-level interface for building phylogenetic trees
"""

from rrna_phylo.core.tree import TreeNode
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
from rrna_phylo.core.builder import PhylogeneticTreeBuilder, build_trees

__all__ = [
    "TreeNode",
    "SequenceType",
    "SequenceTypeDetector",
    "PhylogeneticTreeBuilder",
    "build_trees",
]

"""
Core functionality for phylogenetic tree building.

This package contains:
- tree: TreeNode class for representing phylogenetic trees
- sequence_type: Automatic detection of DNA/RNA/Protein sequences
- builder: High-level interface (will be added during refactoring)
"""

from rrna_phylo.core.tree import TreeNode
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector

__all__ = [
    "TreeNode",
    "SequenceType",
    "SequenceTypeDetector",
]

"""
rRNA-Phylo: General-Purpose Phylogenetic Analysis

A comprehensive phylogenetic tree inference system for DNA, RNA, and Protein sequences.

Quick Start:
    >>> from rrna_phylo import build_trees, Sequence
    >>>
    >>> sequences = [
    ...     Sequence("seq1", "Species 1", "ATGCAT"),
    ...     Sequence("seq2", "Species 2", "ATGCAC"),
    ...     Sequence("seq3", "Species 3", "ATCCAA"),
    ... ]
    >>>
    >>> results = build_trees(sequences)
    >>> upgma = results["upgma"]
    >>> bionj = results["bionj"]
    >>> ml_tree, logL = results["ml"]

Main Components:
    - build_trees(): Main interface for building phylogenetic trees
    - Sequence: Data class for biological sequences
    - TreeNode: Phylogenetic tree data structure
    - SequenceType: Enum for DNA/RNA/Protein types

For more details, see the documentation and examples.
"""

__version__ = "0.1.0"
__author__ = "rRNA-Phylo Contributors"

# Core exports - only import what exists during refactoring
# TEMPORARY: Minimal imports during refactoring
try:
    from rrna_phylo.core.tree import TreeNode
except ImportError:
    pass

try:
    from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
except ImportError:
    pass

try:
    from rrna_phylo.io.fasta_parser import Sequence, FastaParser
except ImportError:
    pass

# NOTE: More imports will be added as refactoring progresses
# from rrna_phylo.core.builder import PhylogeneticTreeBuilder, build_trees

__all__ = [
    "TreeNode",
    "SequenceType",
    "SequenceTypeDetector",
    "Sequence",
    "FastaParser",
    "__version__",
]

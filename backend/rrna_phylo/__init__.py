"""
rRNA-Phylo: Phylogenetic Analysis for rRNA Sequences

A phylogenetic tree inference system for DNA and RNA sequences, optimized for rRNA.

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
    - SequenceType: Enum for DNA/RNA types

For more details, see the documentation and examples.
"""

__version__ = "0.1.0"
__author__ = "rRNA-Phylo Contributors"

# Core API exports
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
from rrna_phylo.core.builder import PhylogeneticTreeBuilder, build_trees
from rrna_phylo.io.fasta_parser import Sequence, FastaParser

# Method exports
from rrna_phylo.methods import build_upgma_tree, build_bionj_tree

# Model exports
from rrna_phylo.models import build_ml_tree_level3, build_ml_tree_level4

# Distance calculation exports
from rrna_phylo.distance import calculate_distance_matrix

__all__ = [
    # Core classes
    "TreeNode",
    "SequenceType",
    "SequenceTypeDetector",
    "PhylogeneticTreeBuilder",
    # Main API functions
    "build_trees",
    # IO classes
    "Sequence",
    "FastaParser",
    # Tree building methods
    "build_upgma_tree",
    "build_bionj_tree",
    "build_ml_tree_level3",
    "build_ml_tree_level4",
    # Distance calculations
    "calculate_distance_matrix",
    # Version
    "__version__",
]

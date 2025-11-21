"""
Input/Output modules for sequence data.

This package handles reading and writing sequence data, including:
- FASTA file parsing
- Multiple sequence alignment (MUSCLE wrapper)
"""

from rrna_phylo.io.fasta_parser import Sequence, FastaParser
from rrna_phylo.io.aligner import MuscleAligner

__all__ = [
    "Sequence",
    "FastaParser",
    "MuscleAligner",
]

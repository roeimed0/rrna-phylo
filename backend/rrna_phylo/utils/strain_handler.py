"""
Strain handling utilities for phylogenetic analysis.

This module provides functions to handle multiple rRNA copies from the same genome,
preventing overrepresentation bias in phylogenetic trees.
"""

from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import re
from rrna_phylo.io.fasta_parser import Sequence


def detect_strain_groups(sequences: List[Sequence]) -> Dict[str, List[Sequence]]:
    """
    Detect which sequences belong to the same strain/genome.

    Groups sequences by:
    1. Genome accession (e.g., U00096 for E. coli K-12)
    2. Species name in description

    Args:
        sequences: List of sequences to group

    Returns:
        Dictionary mapping strain identifier to list of sequences

    Example:
        >>> seqs = [
        ...     Sequence("U00096.123.456", "E. coli K-12", "ATGC"),
        ...     Sequence("U00096.789.012", "E. coli K-12", "ATGC"),
        ...     Sequence("AE004091.111.222", "P. aeruginosa PAO1", "ATGC")
        ... ]
        >>> groups = detect_strain_groups(seqs)
        >>> len(groups["U00096"])  # E. coli copies
        2
    """
    strain_groups = defaultdict(list)

    for seq in sequences:
        # Extract genome accession from sequence ID
        # Pattern: ACCESSION.start.end
        match = re.match(r'^([A-Z]+\d+)', seq.id)
        if match:
            accession = match.group(1)
            strain_groups[accession].append(seq)
        else:
            # Fallback: use full ID if no accession pattern
            strain_groups[seq.id].append(seq)

    return dict(strain_groups)


def select_representative(sequences: List[Sequence], method: str = "longest") -> Sequence:
    """
    Select one representative sequence from multiple copies.

    Args:
        sequences: List of sequences from same strain
        method: Selection method
            - "longest": Select longest sequence (default)
            - "first": Select first sequence
            - "median": Select sequence with median length

    Returns:
        Representative sequence

    Example:
        >>> seqs = [
        ...     Sequence("copy1", "E. coli", "ATGC"),
        ...     Sequence("copy2", "E. coli", "ATGCATGC"),  # longest
        ...     Sequence("copy3", "E. coli", "ATG")
        ... ]
        >>> rep = select_representative(seqs, method="longest")
        >>> rep.id
        'copy2'
    """
    if not sequences:
        raise ValueError("Cannot select representative from empty list")

    if len(sequences) == 1:
        return sequences[0]

    if method == "longest":
        return max(sequences, key=lambda s: len(s.sequence))
    elif method == "first":
        return sequences[0]
    elif method == "median":
        sorted_seqs = sorted(sequences, key=lambda s: len(s.sequence))
        return sorted_seqs[len(sorted_seqs) // 2]
    else:
        raise ValueError(f"Unknown method: {method}")


def create_consensus_sequence(sequences: List[Sequence], strain_id: str) -> Sequence:
    """
    Create consensus sequence from multiple copies using majority voting.

    Args:
        sequences: List of sequences from same strain (must be aligned)
        strain_id: Identifier for the strain

    Returns:
        Consensus sequence

    Raises:
        ValueError: If sequences have different lengths (not aligned)

    Example:
        >>> seqs = [
        ...     Sequence("copy1", "E. coli", "ATGC"),
        ...     Sequence("copy2", "E. coli", "ATGC"),
        ...     Sequence("copy3", "E. coli", "ATCC")
        ... ]
        >>> consensus = create_consensus_sequence(seqs, "E_coli_K12")
        >>> consensus.sequence
        'ATGC'  # Majority at each position
    """
    if not sequences:
        raise ValueError("Cannot create consensus from empty list")

    if len(sequences) == 1:
        return sequences[0]

    # Check all sequences have same length
    lengths = [len(s.sequence) for s in sequences]
    if len(set(lengths)) > 1:
        raise ValueError(
            f"Sequences must be aligned (same length). Got lengths: {set(lengths)}"
        )

    seq_length = lengths[0]
    consensus = []

    for i in range(seq_length):
        # Count nucleotides at position i
        bases = [seq.sequence[i] for seq in sequences]
        # Get most common base (majority vote)
        base_counts = {}
        for base in bases:
            base_counts[base] = base_counts.get(base, 0) + 1

        # Select most common base
        consensus_base = max(base_counts.items(), key=lambda x: x[1])[0]
        consensus.append(consensus_base)

    # Create consensus sequence
    # Extract species name from first sequence
    description = sequences[0].description if sequences[0].description else f"{strain_id} consensus"

    return Sequence(
        id=f"{strain_id}_consensus",
        description=description,
        sequence=''.join(consensus)
    )


def dereplicate_strains(
    sequences: List[Sequence],
    method: str = "representative",
    selection: str = "longest"
) -> Tuple[List[Sequence], Dict[str, List[str]]]:
    """
    Remove redundant sequences from same strain, keeping one per strain.

    This is the main function for handling multi-strain phylogeny.

    Args:
        sequences: List of all sequences
        method: Dereplication method
            - "representative": Keep one representative copy per strain (default)
            - "consensus": Create consensus sequence from all copies
        selection: Method for selecting representative (if method="representative")
            - "longest": Select longest sequence
            - "first": Select first sequence
            - "median": Select median-length sequence

    Returns:
        Tuple of:
        - List of dereplicated sequences (one per strain)
        - Dictionary mapping strain ID to list of original sequence IDs

    Example:
        >>> seqs = [
        ...     Sequence("U00096.1.100", "E. coli K-12", "ATGC"),
        ...     Sequence("U00096.200.300", "E. coli K-12", "ATGC"),
        ...     Sequence("AE004091.1.100", "P. aeruginosa PAO1", "GGCC")
        ... ]
        >>> derep_seqs, mapping = dereplicate_strains(seqs)
        >>> len(derep_seqs)
        2  # One E. coli, one P. aeruginosa
        >>> mapping["U00096"]
        ['U00096.1.100', 'U00096.200.300']
    """
    # Group sequences by strain
    strain_groups = detect_strain_groups(sequences)

    dereplicated = []
    mapping = {}

    for strain_id, strain_seqs in strain_groups.items():
        # Store original IDs
        mapping[strain_id] = [seq.id for seq in strain_seqs]

        # Select or create representative
        if method == "representative":
            representative = select_representative(strain_seqs, method=selection)
        elif method == "consensus":
            representative = create_consensus_sequence(strain_seqs, strain_id)
        else:
            raise ValueError(f"Unknown method: {method}")

        dereplicated.append(representative)

    return dereplicated, mapping


def get_strain_summary(sequences: List[Sequence]) -> str:
    """
    Get human-readable summary of strain groups.

    Args:
        sequences: List of sequences to analyze

    Returns:
        Formatted string with strain summary

    Example:
        >>> summary = get_strain_summary(sequences)
        >>> print(summary)
        Detected 5 strains with 24 total sequences:
          - U00096 (E. coli K-12): 7 copies
          - AE004091 (P. aeruginosa PAO1): 4 copies
          ...
    """
    strain_groups = detect_strain_groups(sequences)

    lines = [
        f"Detected {len(strain_groups)} strain(s) with {len(sequences)} total sequences:"
    ]

    for strain_id, strain_seqs in sorted(strain_groups.items()):
        # Extract species name from first sequence
        if strain_seqs[0].description:
            species = strain_seqs[0].description.split(';')[-1] if ';' in strain_seqs[0].description else strain_seqs[0].description
        else:
            species = "Unknown species"

        lines.append(f"  - {strain_id} ({species}): {len(strain_seqs)} cop{'y' if len(strain_seqs) == 1 else 'ies'}")

    return '\n'.join(lines)


if __name__ == "__main__":
    # Example usage
    test_seqs = [
        Sequence("U00096.1.100", "Escherichia coli K-12", "ATGC"),
        Sequence("U00096.200.300", "Escherichia coli K-12", "ATGCATGC"),
        Sequence("AE004091.1.100", "Pseudomonas aeruginosa PAO1", "GGCC"),
    ]

    print(get_strain_summary(test_seqs))
    print()

    derep, mapping = dereplicate_strains(test_seqs, method="representative", selection="longest")
    print(f"Dereplicated: {len(derep)} sequences")
    for seq in derep:
        print(f"  - {seq.id}: {len(seq.sequence)} bp")

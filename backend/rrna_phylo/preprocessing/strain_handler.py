"""
Strain handling utilities for phylogenetic analysis.

This module provides functions to handle multiple rRNA copies from the same genome,
preventing overrepresentation bias in phylogenetic trees.
"""

from typing import List, Dict, Tuple, Optional, Set
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


def remove_exact_duplicates(sequences: List[Sequence]) -> Tuple[List[Sequence], Dict[str, List[str]]]:
    """
    Remove sequences with identical nucleotide sequences.

    This is the first step of deduplication - remove 100% identical sequences.

    Args:
        sequences: List of sequences

    Returns:
        Tuple of:
        - List of unique sequences (one per unique sequence)
        - Dictionary mapping unique sequence to list of duplicate IDs

    Example:
        >>> seqs = [
        ...     Sequence("copy1", "E. coli", "ATGC"),
        ...     Sequence("copy2", "E. coli", "ATGC"),  # exact duplicate
        ...     Sequence("copy3", "E. coli", "ATGC"),  # exact duplicate
        ...     Sequence("copy4", "P. aeruginosa", "GGCC")
        ... ]
        >>> unique, duplicates = remove_exact_duplicates(seqs)
        >>> len(unique)
        2  # One ATGC, one GGCC
        >>> duplicates["copy1"]
        ['copy1', 'copy2', 'copy3']
    """
    # Map sequence -> first occurrence
    seq_to_first = {}
    # Map first sequence ID -> list of all duplicate IDs
    duplicates_map = defaultdict(list)

    for seq in sequences:
        seq_str = seq.sequence.upper().replace('-', '')  # Remove gaps for comparison

        if seq_str not in seq_to_first:
            # First occurrence of this sequence
            seq_to_first[seq_str] = seq
            duplicates_map[seq.id].append(seq.id)
        else:
            # Duplicate found
            first_seq = seq_to_first[seq_str]
            duplicates_map[first_seq.id].append(seq.id)

    # Return unique sequences
    unique_seqs = list(seq_to_first.values())
    return unique_seqs, dict(duplicates_map)


def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate pairwise sequence identity (percentage of matching positions).

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Similarity as percentage (0-100)

    Example:
        >>> calculate_sequence_similarity("ATGC", "ATGC")
        100.0
        >>> calculate_sequence_similarity("ATGC", "ATGG")
        75.0
    """
    if len(seq1) != len(seq2):
        # For unequal lengths, use the shorter one as denominator
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        return (matches / min_len) * 100.0

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100.0


def cluster_similar_sequences(
    sequences: List[Sequence],
    similarity_threshold: float = 99.5,
    species_aware: bool = True
) -> List[List[Sequence]]:
    """
    Cluster sequences by similarity using simple greedy clustering.

    IMPORTANT: By default, only clusters sequences from the SAME species/genome.
    This prevents accidentally merging phylogenetically distinct sequences.

    Args:
        sequences: List of sequences to cluster
        similarity_threshold: Minimum similarity percentage (default 99.5%)
        species_aware: Only cluster sequences from same species (default True, RECOMMENDED)

    Returns:
        List of clusters, each cluster is a list of similar sequences

    Example:
        >>> seqs = [
        ...     Sequence("seq1", "E. coli", "ATGCATGC"),
        ...     Sequence("seq2", "E. coli", "ATGCATGC"),  # 100% similar to seq1, same species
        ...     Sequence("seq3", "Salmonella", "ATGCATGC"),  # 100% similar but DIFFERENT species
        ... ]
        >>> clusters = cluster_similar_sequences(seqs, similarity_threshold=99.0, species_aware=True)
        >>> len(clusters)
        2  # seq1+seq2 in one cluster, seq3 in another (different species)
    """
    if not sequences:
        return []

    # Clusters: list of lists
    clusters = []
    # Track which sequences are already clustered
    clustered = set()

    for i, seq in enumerate(sequences):
        if seq.id in clustered:
            continue

        # Start new cluster with this sequence
        cluster = [seq]
        clustered.add(seq.id)

        # Get species/genome identifier for this sequence
        if species_aware:
            # Use main accession as species identifier (e.g., U00096 for all E. coli K-12 sequences)
            species_id = seq.main_accession
        else:
            species_id = None

        # Find all sequences similar to this one
        for j, other_seq in enumerate(sequences):
            if i == j or other_seq.id in clustered:
                continue

            # If species-aware, only cluster sequences from same species
            if species_aware:
                if other_seq.main_accession != species_id:
                    continue  # Different species, skip

            similarity = calculate_sequence_similarity(
                seq.sequence.upper(),
                other_seq.sequence.upper()
            )

            if similarity >= similarity_threshold:
                cluster.append(other_seq)
                clustered.add(other_seq.id)

        clusters.append(cluster)

    return clusters


def smart_dereplicate(
    sequences: List[Sequence],
    remove_exact: bool = True,
    similarity_threshold: float = 99.5,
    selection_method: str = "longest",
    species_aware: bool = True,
    verbose: bool = False
) -> Tuple[List[Sequence], Dict[str, any]]:
    """
    Smart multi-step deduplication:
    1. Remove exact duplicates (100% identical sequences)
    2. Cluster highly similar sequences (>= similarity_threshold%) from SAME species
    3. Select one representative per cluster

    This is the recommended deduplication approach for phylogenetic analysis.

    Args:
        sequences: List of sequences to dereplicate
        remove_exact: Remove exact duplicates first (default True)
        similarity_threshold: Cluster sequences with >= this similarity (default 99.5%)
        selection_method: How to pick representative ("longest", "first", "median")
        species_aware: Only cluster sequences from same species (default True, RECOMMENDED)
        verbose: Print deduplication statistics

    Returns:
        Tuple of:
        - List of dereplicated sequences
        - Dictionary with statistics:
            - "original_count": Original number of sequences
            - "exact_duplicates_removed": Number of exact duplicates removed
            - "clusters_formed": Number of similarity clusters
            - "final_count": Final number of sequences
            - "reduction_percentage": Percentage of sequences removed

    Example:
        >>> seqs = [sequences from same species with multiple rRNA copies]
        >>> derep_seqs, stats = smart_dereplicate(seqs, verbose=True)
        >>> print(f"Reduced from {stats['original_count']} to {stats['final_count']} sequences")
    """
    stats = {
        "original_count": len(sequences),
        "exact_duplicates_removed": 0,
        "clusters_formed": 0,
        "final_count": 0,
        "reduction_percentage": 0.0
    }

    if not sequences:
        return [], stats

    working_seqs = sequences

    # Step 1: Remove exact duplicates
    if remove_exact:
        unique_seqs, dup_map = remove_exact_duplicates(working_seqs)
        num_exact_dups = len(working_seqs) - len(unique_seqs)
        stats["exact_duplicates_removed"] = num_exact_dups

        if verbose and num_exact_dups > 0:
            print(f"  Removed {num_exact_dups} exact duplicate(s)")
            print(f"  {len(unique_seqs)} unique sequences remain")

        working_seqs = unique_seqs

    # Step 2: Cluster by similarity (species-aware)
    clusters = cluster_similar_sequences(working_seqs, similarity_threshold, species_aware=species_aware)
    stats["clusters_formed"] = len(clusters)

    if verbose:
        print(f"  Formed {len(clusters)} cluster(s) at {similarity_threshold}% similarity")
        multi_seq_clusters = [c for c in clusters if len(c) > 1]
        if multi_seq_clusters:
            print(f"  {len(multi_seq_clusters)} cluster(s) contain multiple sequences:")
            for i, cluster in enumerate(multi_seq_clusters[:5]):  # Show first 5
                print(f"    Cluster {i+1}: {len(cluster)} sequences")

    # Step 3: Select representatives
    representatives = []
    for cluster in clusters:
        rep = select_representative(cluster, method=selection_method)
        representatives.append(rep)

    stats["final_count"] = len(representatives)
    stats["reduction_percentage"] = ((stats["original_count"] - stats["final_count"]) / stats["original_count"]) * 100

    if verbose:
        print(f"  Final: {stats['final_count']} representative sequences")
        print(f"  Reduced by {stats['reduction_percentage']:.1f}%")

    return representatives, stats


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

"""
Sampling strategies to handle database bias in phylogenetic analysis.

Database bias is a critical problem in phylogenetics:
- Billions of human sequences vs. a handful of rare species
- Model organisms (E. coli, mouse, yeast) massively overrepresented
- Geographic bias (Western countries vs. global diversity)
- Temporal bias (recent sequences vs. historical diversity)

This module provides strategies to create balanced, representative datasets.
"""

from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import random
from rrna_phylo.io.fasta_parser import Sequence


class DatabaseBias:
    """
    Types of database bias and their impact.
    """
    TAXONOMIC = "taxonomic"  # E.g., 1B human vs 10 rare species
    GEOGRAPHIC = "geographic"  # Western countries overrepresented
    TEMPORAL = "temporal"  # Recent sequences vs historical
    METHODOLOGICAL = "methodological"  # Certain sequencing methods overrepresented


def count_by_taxonomy(sequences: List[Sequence]) -> Dict[str, int]:
    """
    Count sequences by taxonomic group.

    Args:
        sequences: List of sequences

    Returns:
        Dictionary mapping taxonomy to count

    Example:
        >>> counts = count_by_taxonomy(sequences)
        >>> counts
        {'Homo sapiens': 1000000000, 'Rare species': 3}
    """
    counts = defaultdict(int)

    for seq in sequences:
        # Extract species from description
        # Format: "Genus species ..." or ";Genus species"
        if not seq.description:
            counts["Unknown"] += 1
            continue

        desc = seq.description
        # Try to extract species name (last element if semicolon-separated)
        if ';' in desc:
            species = desc.split(';')[-1].strip()
        else:
            # Take first two words as "Genus species"
            parts = desc.split()
            species = ' '.join(parts[:2]) if len(parts) >= 2 else desc

        counts[species] += 1

    return dict(counts)


def detect_bias(sequences: List[Sequence], threshold: float = 0.1) -> Dict[str, List[str]]:
    """
    Detect overrepresented taxa that may bias phylogenetic analysis.

    Args:
        sequences: List of sequences
        threshold: Fraction above which a taxon is considered overrepresented (default: 10%)

    Returns:
        Dictionary with bias analysis

    Example:
        >>> bias = detect_bias(sequences, threshold=0.1)
        >>> bias['overrepresented']
        ['Homo sapiens (90% of dataset)']
        >>> bias['underrepresented']
        ['Rare species A (0.0001%)', ...]
    """
    counts = count_by_taxonomy(sequences)
    total = len(sequences)

    overrepresented = []
    underrepresented = []
    balanced = []

    for species, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
        fraction = count / total

        if fraction > threshold:
            overrepresented.append(f"{species} ({fraction * 100:.1f}% of dataset, {count} seqs)")
        elif fraction < 0.01:  # Less than 1%
            underrepresented.append(f"{species} ({fraction * 100:.4f}%, {count} seqs)")
        else:
            balanced.append(f"{species} ({fraction * 100:.1f}%, {count} seqs)")

    return {
        'overrepresented': overrepresented,
        'underrepresented': underrepresented,
        'balanced': balanced,
        'total_species': len(counts),
        'total_sequences': total
    }


def stratified_sample(
    sequences: List[Sequence],
    max_per_species: int = 10,
    min_per_species: int = 1,
    random_seed: Optional[int] = None
) -> Tuple[List[Sequence], Dict[str, int]]:
    """
    Perform stratified sampling to balance overrepresented species.

    This is the key function for handling database bias!

    Strategy:
    1. Group sequences by species
    2. Cap overrepresented species at max_per_species
    3. Keep all sequences from rare species (>= min_per_species)
    4. Randomly sample from abundant species

    Args:
        sequences: List of all sequences
        max_per_species: Maximum sequences per species (default: 10)
        min_per_species: Minimum sequences per species to include (default: 1)
        random_seed: Random seed for reproducibility

    Returns:
        Tuple of (sampled_sequences, sampling_report)

    Example:
        >>> # Dataset: 1000 human, 3 rare species
        >>> sampled, report = stratified_sample(sequences, max_per_species=10)
        >>> # Result: 10 human, 3 rare = balanced!
        >>> report
        {'Homo sapiens': 10, 'Rare species': 3}
    """
    if random_seed is not None:
        random.seed(random_seed)

    # Group by species
    species_groups = defaultdict(list)
    for seq in sequences:
        # Extract species name
        if seq.description and ';' in seq.description:
            species = seq.description.split(';')[-1].strip()
        elif seq.description:
            parts = seq.description.split()
            species = ' '.join(parts[:2]) if len(parts) >= 2 else seq.description
        else:
            species = "Unknown"

        species_groups[species].append(seq)

    # Sample from each species
    sampled = []
    sampling_report = {}

    for species, species_seqs in species_groups.items():
        count = len(species_seqs)

        if count < min_per_species:
            # Skip species with too few sequences
            continue
        elif count <= max_per_species:
            # Keep all sequences from rare/balanced species
            sampled.extend(species_seqs)
            sampling_report[species] = count
        else:
            # Sample from overrepresented species
            sample = random.sample(species_seqs, max_per_species)
            sampled.extend(sample)
            sampling_report[species] = max_per_species

    return sampled, sampling_report


def taxonomic_coverage_sample(
    sequences: List[Sequence],
    target_diversity: int = 100,
    prioritize_rare: bool = True
) -> Tuple[List[Sequence], Dict]:
    """
    Sample to maximize taxonomic diversity while controlling dataset size.

    Strategy:
    1. Identify all unique species
    2. Select representatives to maximize phylogenetic diversity
    3. Prioritize rare species (if prioritize_rare=True)

    Args:
        sequences: List of all sequences
        target_diversity: Target number of species to include
        prioritize_rare: Whether to prioritize rare species

    Returns:
        Tuple of (sampled_sequences, diversity_report)
    """
    # Group by species
    species_groups = defaultdict(list)
    for seq in sequences:
        if seq.description and ';' in seq.description:
            species = seq.description.split(';')[-1].strip()
        elif seq.description:
            parts = seq.description.split()
            species = ' '.join(parts[:2]) if len(parts) >= 2 else seq.description
        else:
            species = "Unknown"

        species_groups[species].append(seq)

    # Sort species by abundance
    sorted_species = sorted(species_groups.items(), key=lambda x: len(x[1]))

    if prioritize_rare:
        # Take rarest species first
        selected_species = sorted_species[:target_diversity]
    else:
        # Take most abundant species first
        selected_species = sorted_species[-target_diversity:]

    # Select one representative per species (longest sequence)
    sampled = []
    for species, seqs in selected_species:
        # Select longest sequence
        representative = max(seqs, key=lambda s: len(s.sequence))
        sampled.append(representative)

    diversity_report = {
        'total_species_available': len(species_groups),
        'species_selected': len(selected_species),
        'total_sequences_available': len(sequences),
        'sequences_selected': len(sampled),
        'sampling_strategy': 'rare-first' if prioritize_rare else 'abundant-first'
    }

    return sampled, diversity_report


def geographic_balance_sample(
    sequences: List[Sequence],
    max_per_region: int = 50
) -> Tuple[List[Sequence], Dict[str, int]]:
    """
    Balance geographic representation in the dataset.

    Args:
        sequences: List of sequences
        max_per_region: Maximum sequences per geographic region

    Returns:
        Tuple of (sampled_sequences, region_report)

    Note:
        Requires geographic information in sequence descriptions.
        Format: "...;Country;..." or "...;Continent;..."
    """
    # This is a placeholder - would require geographic parsing
    # Most sequence databases don't have standardized geo info
    # Would need to integrate with external databases (INSDC, ENA)

    # For now, just return all sequences
    return sequences, {"note": "Geographic balancing not implemented"}


def get_sampling_recommendations(sequences: List[Sequence]) -> str:
    """
    Analyze dataset and recommend sampling strategy.

    Args:
        sequences: List of sequences to analyze

    Returns:
        Human-readable recommendations
    """
    bias = detect_bias(sequences, threshold=0.1)

    recommendations = ["=== SAMPLING RECOMMENDATIONS ===\n"]

    recommendations.append(f"Dataset: {bias['total_sequences']} sequences from {bias['total_species']} species\n")

    if bias['overrepresented']:
        recommendations.append("[!] OVERREPRESENTED SPECIES DETECTED:")
        for species in bias['overrepresented'][:5]:  # Top 5
            recommendations.append(f"  - {species}")
        recommendations.append("")
        recommendations.append("RECOMMENDATION: Use stratified sampling")
        recommendations.append("  Command: --stratify --max-per-species 10")
        recommendations.append("")

    if bias['underrepresented']:
        count = len(bias['underrepresented'])
        recommendations.append(f"[INFO] {count} rare species detected (each <1% of dataset)")
        recommendations.append("RECOMMENDATION: Ensure rare species are preserved")
        recommendations.append("  Command: --stratify --min-per-species 1 --max-per-species 10")
        recommendations.append("")

    if len(bias['balanced']) > 0:
        recommendations.append(f"[OK] {len(bias['balanced'])} species are reasonably balanced")
        recommendations.append("")

    # Overall recommendation
    if bias['overrepresented']:
        recommendations.append("=== RECOMMENDED WORKFLOW ===")
        recommendations.append("1. Dereplicate: --dereplicate (remove duplicate gene copies)")
        recommendations.append("2. Stratify: --stratify --max-per-species 10 (balance species)")
        recommendations.append("3. Add outgroup: --outgroup \"pattern\" (root the tree)")
        recommendations.append("4. Bootstrap: --bootstrap 100 (assess confidence)")
    else:
        recommendations.append("[OK] Dataset appears balanced - no stratification needed")

    return '\n'.join(recommendations)


if __name__ == "__main__":
    # Example with biased dataset
    from rrna_phylo.io.fasta_parser import Sequence

    # Simulate heavily biased dataset (like real NCBI)
    sequences = []

    # 1000 E. coli sequences (model organism)
    for i in range(1000):
        sequences.append(Sequence(f"ecoli_{i}", "Escherichia coli K-12", "ATGC" * 100))

    # 500 human sequences
    for i in range(500):
        sequences.append(Sequence(f"human_{i}", "Homo sapiens", "ATGC" * 100))

    # 3 rare species (just a few sequences each)
    sequences.append(Sequence("rare1", "Rare species A", "ATGC" * 100))
    sequences.append(Sequence("rare2", "Rare species B", "ATGC" * 100))
    sequences.append(Sequence("rare3", "Rare species C", "ATGC" * 100))

    print(get_sampling_recommendations(sequences))
    print("\n" + "=" * 70 + "\n")

    # Demonstrate stratified sampling
    sampled, report = stratified_sample(sequences, max_per_species=10)
    print(f"Original dataset: {len(sequences)} sequences")
    print(f"After stratified sampling: {len(sampled)} sequences")
    print("\nSampling report:")
    for species, count in sorted(report.items(), key=lambda x: x[1], reverse=True):
        print(f"  {species}: {count} sequences")

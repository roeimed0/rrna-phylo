#!/usr/bin/env python3
"""
Create a synthetic larger test dataset by replicating and mutating existing sequences.

This creates a realistic test dataset with:
- 100-150 sequences derived from real 16S rRNA
- Varying levels of sequence divergence (0.1% to 15%)
- Multiple copies per species (simulates rRNA operons)
- Species from different phyla

This is useful for testing:
1. Deduplication effectiveness at scale
2. Tree building performance with larger datasets
3. Bootstrap reliability
4. Memory usage and computational bottlenecks
"""

import random
from pathlib import Path
import sys

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rrna_phylo.io.fasta_parser import parse_fasta

def mutate_sequence(seq, mutation_rate=0.01):
    """
    Introduce random mutations into a sequence.

    Args:
        seq: DNA sequence string
        mutation_rate: Fraction of bases to mutate (0.01 = 1%)

    Returns:
        Mutated sequence string
    """
    bases = ['A', 'T', 'G', 'C']
    seq_list = list(seq.upper())
    n_mutations = int(len(seq_list) * mutation_rate)

    # Select random positions to mutate
    positions = random.sample(range(len(seq_list)), n_mutations)

    for pos in positions:
        current_base = seq_list[pos]
        # Choose different base
        new_base = random.choice([b for b in bases if b != current_base])
        seq_list[pos] = new_base

    return ''.join(seq_list)

def create_species_variants(original_seq, species_name, accession, n_variants=10,
                          divergence_range=(0.001, 0.05)):
    """
    Create multiple variants of a species with varying divergence.

    Simulates:
    - Multiple rRNA operons (very similar, <0.5% divergence)
    - Different strains (0.5-2% divergence)
    - Related species (2-5% divergence)
    """
    variants = []

    for i in range(n_variants):
        # Divergence increases with variant number
        # First few are very similar (operons), later ones more divergent (strains/species)
        if i < 3:
            # Simulate rRNA operons within same genome (0.1-0.5%)
            mutation_rate = random.uniform(0.001, 0.005)
            variant_name = f"{species_name} (operon {i+1})"
        elif i < 6:
            # Simulate different strains (0.5-2%)
            mutation_rate = random.uniform(0.005, 0.02)
            variant_name = f"{species_name} strain {i-2}"
        else:
            # Simulate related species/subspecies (2-5%)
            mutation_rate = random.uniform(0.02, divergence_range[1])
            variant_name = f"{species_name} variant {i-5}"

        mutated_seq = mutate_sequence(original_seq, mutation_rate)

        # Create unique accession
        variant_accession = f"{accession}.VAR{i+1:03d}"

        variants.append({
            'name': variant_name,
            'accession': variant_accession,
            'sequence': mutated_seq,
            'divergence': mutation_rate * 100  # Convert to percentage
        })

    return variants

def main():
    """Create synthetic larger dataset."""

    # Set random seed for reproducibility
    random.seed(42)

    # Load existing test dataset
    input_file = Path(__file__).parent.parent / "test_real_rrana.fasta"
    output_file = Path(__file__).parent.parent / "test_large_synthetic.fasta"

    print("=" * 80)
    print("CREATING SYNTHETIC LARGER DATASET")
    print("=" * 80)
    print()

    print(f"Loading sequences from: {input_file}")
    sequences = list(parse_fasta(input_file))
    print(f"Loaded {len(sequences)} original sequences")
    print()

    # Create variants for each sequence
    all_variants = []

    print("Creating variants...")
    print("-" * 80)

    for seq in sequences:
        # Extract species name (simplified)
        desc_parts = seq.description.split()
        if len(desc_parts) >= 2:
            species_name = ' '.join(desc_parts[:2])  # First two words
        else:
            species_name = desc_parts[0] if desc_parts else "Unknown"

        # Determine number of variants based on sequence
        # Create more variants for some species to simulate real-world bias
        if "Escherichia" in species_name or "Salmonella" in species_name:
            n_variants = 12  # Well-studied organisms have more sequences
        elif "Pseudomonas" in species_name or "Staphylococcus" in species_name:
            n_variants = 10
        else:
            n_variants = 8

        variants = create_species_variants(
            str(seq.sequence),
            species_name,
            seq.main_accession,  # Use main_accession property instead of accession
            n_variants=n_variants,
            divergence_range=(0.001, 0.05)
        )

        all_variants.extend(variants)
        print(f"  {species_name}: {len(variants)} variants")

    print()
    print(f"Total variants created: {len(all_variants)}")
    print()

    # Write to FASTA
    print(f"Writing to: {output_file}")
    with open(output_file, 'w') as f:
        for i, variant in enumerate(all_variants, 1):
            # FASTA format: >accession description
            header = f">{variant['accession']} {variant['name']}"
            f.write(header + "\n")

            # Write sequence in 80-character lines
            seq = variant['sequence']
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + "\n")

    print()
    print("=" * 80)
    print("SYNTHETIC DATASET CREATED!")
    print("=" * 80)
    print()

    # Statistics
    lengths = [len(v['sequence']) for v in all_variants]
    divergences = [v['divergence'] for v in all_variants]

    print("Dataset Statistics:")
    print(f"  Total sequences: {len(all_variants)}")
    print()
    print("  Sequence lengths:")
    print(f"    Min: {min(lengths)} bp")
    print(f"    Max: {max(lengths)} bp")
    print(f"    Mean: {sum(lengths) / len(lengths):.0f} bp")
    print()
    print("  Divergence from originals:")
    print(f"    Min: {min(divergences):.3f}%")
    print(f"    Max: {max(divergences):.3f}%")
    print(f"    Mean: {sum(divergences) / len(divergences):.3f}%")
    print()

    # Deduplication predictions
    print("Expected deduplication results:")

    # Count sequences by similarity level
    very_similar = sum(1 for d in divergences if d < 0.5)
    similar = sum(1 for d in divergences if 0.5 <= d < 2.0)
    divergent = sum(1 for d in divergences if d >= 2.0)

    print(f"  Very similar (<0.5% divergence): {very_similar} sequences")
    print(f"  Similar (0.5-2% divergence): {similar} sequences")
    print(f"  Divergent (>2% divergence): {divergent} sequences")
    print()

    # Estimate deduplication effectiveness
    # With exact duplicate removal: should remove very few (our mutations prevent exact matches)
    # With 99.5% similarity clustering: should cluster most of "very similar" group
    expected_after_exact = len(all_variants)  # No exact duplicates due to mutations
    expected_after_similarity = divergent + similar + (very_similar // 3)  # Rough estimate

    print(f"  Expected after exact deduplication: ~{expected_after_exact} sequences")
    print(f"  Expected after similarity clustering (99.5%): ~{expected_after_similarity} sequences")
    print(f"  Expected reduction: ~{100 * (1 - expected_after_similarity / expected_after_exact):.1f}%")
    print()

    print(f"Output file: {output_file}")
    print()

if __name__ == "__main__":
    main()

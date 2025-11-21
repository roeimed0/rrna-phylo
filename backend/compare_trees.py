"""
Compare UPGMA and BioNJ trees on the same input data.

This demonstrates the difference between methods that assume
a molecular clock (UPGMA) vs those that handle variable rates (BioNJ).
"""

import numpy as np
from fasta_parser import Sequence
from distance import calculate_distance_matrix
from upgma import build_upgma_tree, UPGMABuilder
from bionj import build_bionj_tree, BioNJBuilder


def compare_simple_example():
    """Compare with simple synthetic data."""
    print("=" * 80)
    print("COMPARISON 1: Simple Synthetic Data")
    print("=" * 80)

    # Distance matrix where one lineage evolved much faster
    # A and B are close (0.1), C is medium (0.5), D evolved very fast (0.9)
    distance_matrix = np.array([
        [0.0, 0.1, 0.5, 0.9],
        [0.1, 0.0, 0.5, 0.9],
        [0.5, 0.5, 0.0, 0.8],
        [0.9, 0.9, 0.8, 0.0]
    ])

    labels = ["A_recent", "B_recent", "C_medium", "D_fast_evolving"]

    print("\nDistance Matrix (D evolved much faster):")
    print(f"{'':20}", end='')
    for label in labels:
        print(f"{label:20}", end='')
    print()
    for i, label in enumerate(labels):
        print(f"{label:20}", end='')
        for j in range(len(labels)):
            print(f"{distance_matrix[i][j]:20.3f}", end='')
        print()

    # Build both trees
    print("\n" + "=" * 80)
    print("UPGMA Tree (assumes constant evolutionary rate - MOLECULAR CLOCK)")
    print("=" * 80)
    upgma_tree = build_upgma_tree(distance_matrix, labels)
    upgma_builder = UPGMABuilder()
    upgma_builder.print_tree(upgma_tree)
    print(f"\nNewick: {upgma_tree.to_newick()};")

    print("\n" + "=" * 80)
    print("BioNJ Tree (handles variable rates - NO MOLECULAR CLOCK)")
    print("=" * 80)
    bionj_tree = build_bionj_tree(distance_matrix, labels)
    bionj_builder = BioNJBuilder()
    bionj_builder.print_tree(bionj_tree)
    print(f"\nNewick: {bionj_tree.to_newick()};")

    print("\n" + "=" * 80)
    print("ANALYSIS:")
    print("=" * 80)
    print("- Notice branch lengths differ between methods")
    print("- UPGMA forces equal distances to tips (molecular clock)")
    print("- BioNJ allows D to have longer branch (realistic!)")
    print("- BioNJ is more accurate when evolution rates vary")


def compare_real_sequences():
    """Compare with real sequences."""
    print("\n\n" + "=" * 80)
    print("COMPARISON 2: Real Sequence Data")
    print("=" * 80)

    # Create test sequences (simulating different species)
    sequences = [
        Sequence("E_coli", "Escherichia coli", "ATGCATGCATGCATGC"),
        Sequence("Salmonella", "Salmonella enterica", "ATGCATGCATGCATCC"),  # 1 diff
        Sequence("B_subtilis", "Bacillus subtilis", "ATCCATGCATGCATGC"),  # 2 diff
        Sequence("S_aureus", "Staphylococcus aureus", "ATGCATCCATGCATGC"),  # 2 diff
        Sequence("P_aeruginosa", "Pseudomonas aeruginosa", "TTGCATGCATGCATGC"),  # 2 diff
    ]

    print("\nSequences:")
    for seq in sequences:
        print(f"  {seq.id:20} {seq.sequence}")

    # Calculate distance matrix
    print("\nCalculating Jukes-Cantor distances...")
    dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

    print("\nDistance Matrix:")
    print(f"{'':20}", end='')
    for id in ids:
        print(f"{id:20}", end='')
    print()
    for i, id in enumerate(ids):
        print(f"{id:20}", end='')
        for j in range(len(ids)):
            print(f"{dist_matrix[i][j]:20.6f}", end='')
        print()

    # Build both trees
    print("\n" + "=" * 80)
    print("UPGMA Tree")
    print("=" * 80)
    upgma_tree = build_upgma_tree(dist_matrix, ids)
    upgma_builder = UPGMABuilder()
    upgma_builder.print_tree(upgma_tree)
    print(f"\nNewick: {upgma_tree.to_newick()};")

    print("\n" + "=" * 80)
    print("BioNJ Tree")
    print("=" * 80)
    bionj_tree = build_bionj_tree(dist_matrix, ids)
    bionj_builder = BioNJBuilder()
    bionj_builder.print_tree(bionj_tree)
    print(f"\nNewick: {bionj_tree.to_newick()};")

    print("\n" + "=" * 80)
    print("ANALYSIS:")
    print("=" * 80)
    print("- Both methods cluster closely related species together")
    print("- Branch lengths may differ due to different assumptions")
    print("- BioNJ generally more accurate for phylogenetics")


def compare_extreme_case():
    """Compare with extreme rate variation."""
    print("\n\n" + "=" * 80)
    print("COMPARISON 3: Extreme Rate Variation")
    print("=" * 80)

    # One lineage evolved EXTREMELY fast
    distance_matrix = np.array([
        [0.0, 0.05, 0.10, 0.95],
        [0.05, 0.0, 0.10, 0.95],
        [0.10, 0.10, 0.0, 0.95],
        [0.95, 0.95, 0.95, 0.0]
    ])

    labels = ["Human", "Chimp", "Gorilla", "FastBug"]

    print("\nScenario: Three primates + one rapidly evolving organism")
    print("\nDistance Matrix:")
    print(f"{'':15}", end='')
    for label in labels:
        print(f"{label:15}", end='')
    print()
    for i, label in enumerate(labels):
        print(f"{label:15}", end='')
        for j in range(len(labels)):
            print(f"{distance_matrix[i][j]:15.3f}", end='')
        print()

    # Build both trees
    print("\n" + "=" * 80)
    print("UPGMA Tree (WRONG - forced molecular clock)")
    print("=" * 80)
    upgma_tree = build_upgma_tree(distance_matrix, labels)
    upgma_builder = UPGMABuilder()
    upgma_builder.print_tree(upgma_tree)

    print("\n" + "=" * 80)
    print("BioNJ Tree (CORRECT - handles rate variation)")
    print("=" * 80)
    bionj_tree = build_bionj_tree(distance_matrix, labels)
    bionj_builder = BioNJBuilder()
    bionj_builder.print_tree(bionj_tree)

    print("\n" + "=" * 80)
    print("ANALYSIS:")
    print("=" * 80)
    print("- UPGMA struggles with extreme rate variation")
    print("- BioNJ correctly identifies FastBug as long branch")
    print("- This is why BioNJ is preferred for phylogenetics!")
    print("- For forensics: USE BIONJ, not UPGMA")


def main():
    """Run all comparisons."""
    print("\n" + "=" * 80)
    print("TREE COMPARISON: UPGMA vs BioNJ")
    print("=" * 80)
    print("\nThis demonstrates why BioNJ is better than UPGMA for phylogenetics.")
    print("UPGMA assumes molecular clock (constant rates) - rarely true!")
    print("BioNJ handles variable rates - more realistic and accurate.")

    compare_simple_example()
    compare_real_sequences()
    compare_extreme_case()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("✓ UPGMA: Simple, fast, but assumes molecular clock")
    print("✓ BioNJ: More accurate, handles variable evolutionary rates")
    print("✓ For forensics/publication: Use BioNJ (or better: ML/Bayesian)")
    print("\nNext steps: Add Maximum Likelihood (RAxML) and Bayesian (MrBayes)")
    print("for even better trees, then build consensus!")
    print()


if __name__ == "__main__":
    main()

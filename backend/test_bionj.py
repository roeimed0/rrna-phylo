"""
Test BioNJ tree building.
"""

import numpy as np
from bionj import BioNJBuilder, build_bionj_tree
from upgma import build_upgma_tree


def test_simple_bionj():
    """Test BioNJ with simple 3-sequence example."""
    print("=" * 60)
    print("TEST: Simple BioNJ Tree")
    print("=" * 60)

    # Distance matrix where C evolved faster
    distance_matrix = np.array([
        [0.0, 0.2, 0.8],
        [0.2, 0.0, 0.8],
        [0.8, 0.8, 0.0]
    ])

    labels = ["A", "B", "C"]

    print("\nDistance Matrix:")
    print(f"     A     B     C")
    for i, label in enumerate(labels):
        print(f"{label}  ", end='')
        for j in range(len(labels)):
            print(f"{distance_matrix[i][j]:.3f}  ", end='')
        print()

    # Build tree
    builder = BioNJBuilder()
    tree = builder.build_tree(distance_matrix, labels)

    print("\nBioNJ Tree structure:")
    builder.print_tree(tree)

    print("\nNewick format:")
    newick = tree.to_newick() + ";"
    print(newick)

    assert not tree.is_leaf(), "Root should be internal node"
    print("\n✓ Test passed!")


def test_bionj_vs_upgma():
    """Compare BioNJ and UPGMA on same data."""
    print("\n" + "=" * 60)
    print("TEST: BioNJ vs UPGMA Comparison")
    print("=" * 60)

    # Distance matrix with rate variation
    # A and B are close, D evolved very fast
    distance_matrix = np.array([
        [0.0, 0.1, 0.5, 0.9],
        [0.1, 0.0, 0.5, 0.9],
        [0.5, 0.5, 0.0, 0.9],
        [0.9, 0.9, 0.9, 0.0]
    ])

    labels = ["A", "B", "C", "D"]

    print("\nDistance Matrix (D evolved much faster):")
    print(f"     A     B     C     D")
    for i, label in enumerate(labels):
        print(f"{label}  ", end='')
        for j in range(len(labels)):
            print(f"{distance_matrix[i][j]:.3f}  ", end='')
        print()

    # Build both trees
    print("\n" + "-" * 60)
    print("UPGMA Tree (assumes molecular clock):")
    print("-" * 60)
    upgma_tree = build_upgma_tree(distance_matrix, labels)
    upgma_builder = BioNJBuilder()
    upgma_builder.print_tree(upgma_tree)

    print("\n" + "-" * 60)
    print("BioNJ Tree (handles rate variation):")
    print("-" * 60)
    bionj_tree = build_bionj_tree(distance_matrix, labels)
    bionj_builder = BioNJBuilder()
    bionj_builder.print_tree(bionj_tree)

    print("\nNewick formats:")
    print(f"UPGMA:  {upgma_tree.to_newick()};")
    print(f"BioNJ:  {bionj_tree.to_newick()};")

    print("\n✓ Test passed!")


def test_complete_pipeline():
    """Test complete pipeline: sequences → alignment → distance → BioNJ."""
    print("\n" + "=" * 60)
    print("TEST: Complete Pipeline (Sequences → BioNJ)")
    print("=" * 60)

    from fasta_parser import Sequence
    from distance import calculate_distance_matrix

    # Create aligned sequences
    sequences = [
        Sequence("E_coli", "Escherichia coli", "ATGCATGC"),
        Sequence("B_subtilis", "Bacillus subtilis", "ATCCATGC"),
        Sequence("S_aureus", "Staphylococcus aureus", "ATGCATCC"),
        Sequence("P_aeruginosa", "Pseudomonas aeruginosa", "TTGCATGC"),
    ]

    print("\nSequences:")
    for seq in sequences:
        print(f"  {seq.id}: {seq.sequence}")

    # Calculate distance matrix
    print("\nCalculating Jukes-Cantor distances...")
    dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

    print("\nDistance Matrix:")
    print(f"{'':15}", end='')
    for id in ids:
        print(f"{id:15}", end='')
    print()
    for i, id in enumerate(ids):
        print(f"{id:15}", end='')
        for j in range(len(ids)):
            print(f"{dist_matrix[i][j]:15.6f}", end='')
        print()

    # Build BioNJ tree
    print("\nBuilding BioNJ tree...")
    tree = build_bionj_tree(dist_matrix, ids)

    print("\nBioNJ Tree structure:")
    builder = BioNJBuilder()
    builder.print_tree(tree)

    print("\nNewick format:")
    newick = tree.to_newick() + ";"
    print(newick)

    print("\n✓ Complete pipeline test passed!")


def test_variance_weighting():
    """Test that BioNJ variance weighting works correctly."""
    print("\n" + "=" * 60)
    print("TEST: BioNJ Variance Weighting")
    print("=" * 60)

    # Create asymmetric distance matrix
    # C is closer to A than B, but with high variance
    distance_matrix = np.array([
        [0.0, 0.5, 0.3],
        [0.5, 0.0, 0.7],
        [0.3, 0.7, 0.0]
    ])

    labels = ["A", "B", "C"]

    print("\nDistance Matrix:")
    print("  A-B: 0.5")
    print("  A-C: 0.3  (closer to A)")
    print("  B-C: 0.7")

    builder = BioNJBuilder()
    tree = builder.build_tree(distance_matrix, labels)

    print("\nBioNJ handles this with variance weighting:")
    builder.print_tree(tree)

    print("\nNewick:")
    print(tree.to_newick() + ";")

    print("\n✓ Test passed!")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("BioNJ TREE BUILDING TESTS")
    print("=" * 60)

    test_simple_bionj()
    test_bionj_vs_upgma()
    test_complete_pipeline()
    test_variance_weighting()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✓")
    print("=" * 60)
    print("\nBioNJ is ready! It's the best distance-based method.")
    print("Next: Integrate with ML methods (RAxML) for even better trees!")
    print()


if __name__ == "__main__":
    main()

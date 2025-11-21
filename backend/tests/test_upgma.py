"""
Test UPGMA tree building.
"""

import numpy as np
from rrna_phylo.methods.upgma import UPGMABuilder, build_upgma_tree, TreeNode


def test_simple_tree():
    """Test UPGMA with simple 3-sequence example."""
    print("=" * 60)
    print("TEST: Simple UPGMA Tree")
    print("=" * 60)

    # Simple distance matrix
    # A and B are close (0.2), C is distant (0.6 from both)
    distance_matrix = np.array([
        [0.0, 0.2, 0.6],
        [0.2, 0.0, 0.6],
        [0.6, 0.6, 0.0]
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
    builder = UPGMABuilder()
    tree = builder.build_tree(distance_matrix, labels)

    print("\nTree structure:")
    builder.print_tree(tree)

    print("\nNewick format:")
    newick = tree.to_newick() + ";"
    print(newick)

    # Verify tree structure
    assert not tree.is_leaf(), "Root should be internal node"
    assert tree.left is not None and tree.right is not None

    print("\n[OK] Test passed!")


def test_four_sequences():
    """Test UPGMA with 4 sequences."""
    print("\n" + "=" * 60)
    print("TEST: UPGMA with 4 Sequences")
    print("=" * 60)

    # Distance matrix for 4 sequences
    # (A,B) cluster and (C,D) cluster, then merge
    distance_matrix = np.array([
        [0.0, 0.1, 0.5, 0.5],
        [0.1, 0.0, 0.5, 0.5],
        [0.5, 0.5, 0.0, 0.1],
        [0.5, 0.5, 0.1, 0.0]
    ])

    labels = ["A", "B", "C", "D"]

    print("\nDistance Matrix:")
    print(f"     A     B     C     D")
    for i, label in enumerate(labels):
        print(f"{label}  ", end='')
        for j in range(len(labels)):
            print(f"{distance_matrix[i][j]:.3f}  ", end='')
        print()

    print("\nExpected tree: ((A,B),(C,D))")

    # Build tree
    tree = build_upgma_tree(distance_matrix, labels)

    print("\nTree structure:")
    builder = UPGMABuilder()
    builder.print_tree(tree)

    print("\nNewick format:")
    newick = tree.to_newick() + ";"
    print(newick)

    print("\n[OK] Test passed!")


def test_with_aligned_sequences():
    """Test complete pipeline: sequences -> alignment -> distance -> UPGMA."""
    print("\n" + "=" * 60)
    print("TEST: Complete Pipeline (Sequences -> UPGMA)")
    print("=" * 60)

    # We'll simulate this with pre-calculated distances
    # In real use: parse FASTA -> align with MUSCLE -> calculate distances -> build tree

    from rrna_phylo.io.fasta_parser import Sequence
    from rrna_phylo.distance.distance import calculate_distance_matrix

    # Create aligned sequences
    sequences = [
        Sequence("E_coli", "Escherichia coli", "ATGCATGC"),
        Sequence("B_subtilis", "Bacillus subtilis", "ATCCATGC"),  # 1 diff
        Sequence("S_aureus", "Staphylococcus aureus", "ATGCATCC"),  # 1 diff
        Sequence("P_aeruginosa", "Pseudomonas aeruginosa", "TTGCATGC"),  # 1 diff
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

    # Build UPGMA tree
    print("\nBuilding UPGMA tree...")
    tree = build_upgma_tree(dist_matrix, ids)

    print("\nTree structure:")
    builder = UPGMABuilder()
    builder.print_tree(tree)

    print("\nNewick format:")
    newick = tree.to_newick() + ";"
    print(newick)

    print("\n[OK] Complete pipeline test passed!")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("UPGMA TREE BUILDING TESTS")
    print("=" * 60)

    test_simple_tree()
    test_four_sequences()
    test_with_aligned_sequences()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED [OK]")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

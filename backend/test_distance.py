"""
Test distance calculations.
"""

from rrna_phylo.distance.distance import DistanceCalculator, calculate_distance_matrix
from rrna_phylo.io.fasta_parser import Sequence


def test_p_distance():
    """Test simple p-distance calculation."""
    print("=" * 60)
    print("TEST: P-Distance Calculation")
    print("=" * 60)

    # Create two sequences with known differences
    seq1 = Sequence("seq1", "Test 1", "ATGCATGC")
    seq2 = Sequence("seq2", "Test 2", "ATCCATGC")
    #                                    ^  (1 difference out of 8)

    calc = DistanceCalculator(model="p-distance")
    dist = calc.pairwise_distance(seq1, seq2)

    print(f"\nSeq1: {seq1.sequence}")
    print(f"Seq2: {seq2.sequence}")
    print(f"P-distance: {dist:.4f}")
    print(f"Expected: 0.1250 (1/8)")

    assert abs(dist - 0.125) < 0.001

    print("\n✓ Test passed!")


def test_jukes_cantor():
    """Test Jukes-Cantor corrected distance."""
    print("\n" + "=" * 60)
    print("TEST: Jukes-Cantor Distance")
    print("=" * 60)

    seq1 = Sequence("seq1", "Test 1", "ATGCATGC")
    seq2 = Sequence("seq2", "Test 2", "ATCCATGC")

    calc = DistanceCalculator(model="jukes-cantor")
    dist = calc.pairwise_distance(seq1, seq2)

    print(f"\nSeq1: {seq1.sequence}")
    print(f"Seq2: {seq2.sequence}")
    print(f"P-distance: 0.1250")
    print(f"Jukes-Cantor: {dist:.4f}")
    print(f"(JC > p-distance because it corrects for hidden substitutions)")

    assert dist > 0.125  # JC should be larger than p-distance

    print("\n✓ Test passed!")


def test_distance_with_gaps():
    """Test distance calculation with gaps."""
    print("\n" + "=" * 60)
    print("TEST: Distance with Gaps")
    print("=" * 60)

    seq1 = Sequence("seq1", "Test 1", "ATGC-ATGC")
    seq2 = Sequence("seq2", "Test 2", "ATCC-ATGC")
    #                                    ^  (1 diff, 1 gap ignored)

    calc = DistanceCalculator(model="p-distance")
    dist = calc.pairwise_distance(seq1, seq2)

    print(f"\nSeq1: {seq1.sequence}")
    print(f"Seq2: {seq2.sequence}")
    print(f"Valid positions: 8 (gap ignored)")
    print(f"Differences: 1")
    print(f"P-distance: {dist:.4f}")

    assert abs(dist - 0.125) < 0.001

    print("\n✓ Test passed!")


def test_distance_matrix():
    """Test full distance matrix calculation."""
    print("\n" + "=" * 60)
    print("TEST: Distance Matrix")
    print("=" * 60)

    # Create 4 aligned sequences
    sequences = [
        Sequence("seq1", "Species A", "ATGCATGC"),
        Sequence("seq2", "Species B", "ATCCATGC"),  # 1 diff from seq1
        Sequence("seq3", "Species C", "ATGCATCC"),  # 1 diff from seq1
        Sequence("seq4", "Species D", "TTGCATGC"),  # 1 diff from seq1
    ]

    print("\nSequences:")
    for seq in sequences:
        print(f"  {seq.id}: {seq.sequence}")

    calc = DistanceCalculator(model="jukes-cantor")
    matrix, ids = calc.distance_matrix(sequences)

    calc.print_matrix(matrix, ids)

    # Check matrix properties
    print("\nMatrix properties:")
    print(f"  Shape: {matrix.shape}")
    print(f"  Diagonal (should be all zeros): {matrix.diagonal()}")
    print(f"  Symmetric: {(matrix == matrix.T).all()}")

    assert matrix.shape == (4, 4)
    assert (matrix.diagonal() == 0).all()
    assert (matrix == matrix.T).all()

    print("\n✓ Test passed!")


def test_convenience_function():
    """Test convenience function."""
    print("\n" + "=" * 60)
    print("TEST: Convenience Function")
    print("=" * 60)

    sequences = [
        Sequence("A", "Test A", "ATGCAT"),
        Sequence("B", "Test B", "ATCCAT"),
        Sequence("C", "Test C", "TTGCAT"),
    ]

    matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

    print(f"\n✓ Calculated {len(ids)}x{len(ids)} distance matrix")
    print(f"IDs: {ids}")

    assert len(ids) == 3
    assert matrix.shape == (3, 3)

    print("\n✓ Test passed!")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("DISTANCE CALCULATION TESTS")
    print("=" * 60)

    test_p_distance()
    test_jukes_cantor()
    test_distance_with_gaps()
    test_distance_matrix()
    test_convenience_function()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✓")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

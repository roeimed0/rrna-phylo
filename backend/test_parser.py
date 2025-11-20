"""
Test the FASTA parser.
"""

import os
from fasta_parser import FastaParser, parse_fasta, Sequence


def create_test_fasta(filepath: str, aligned: bool = False):
    """Create a test FASTA file."""
    if aligned:
        # Aligned sequences (same length with gaps)
        content = """>seq1 Escherichia coli
ATGCATGC-ATGC
>seq2 Bacillus subtilis
ATGC-TGCGATGC
>seq3 Staphylococcus aureus
ATGCATGCGATGC
>seq4 Pseudomonas aeruginosa
ATG-ATGC-ATGC
"""
    else:
        # Unaligned sequences (different lengths)
        content = """>seq1 Escherichia coli
ATGCATGCATGC
>seq2 Bacillus subtilis
ATGCTGCGATGC
>seq3 Staphylococcus aureus
ATGCATGCGATGC
>seq4 Pseudomonas aeruginosa
ATGATGCATGC
"""

    with open(filepath, "w", newline="\n") as f:
        f.write(content)


def test_parser():
    """Test basic parsing functionality."""
    print("=" * 60)
    print("TEST: Basic FASTA Parsing")
    print("=" * 60)

    test_file = "test_sequences.fasta"
    create_test_fasta(test_file, aligned=False)

    parser = FastaParser()
    sequences = parser.parse(test_file)

    print(f"\n✓ Parsed {len(sequences)} sequences")

    for i, seq in enumerate(sequences, 1):
        print(f"\n  Sequence {i}:")
        print(f"    ID: {seq.id}")
        print(f"    Description: {seq.description}")
        print(f"    Length: {seq.length} bp")
        print(f"    Sequence: {seq.sequence}")
        print(f"    Is nucleotide: {seq.is_nucleotide()}")
        print(f"    Is protein: {seq.is_protein()}")

    assert len(sequences) == 4

    os.remove(test_file)
    print("\n✓ Test passed!")


def test_validation_unaligned():
    """Test validation of unaligned sequences."""
    print("\n" + "=" * 60)
    print("TEST: Validation - Unaligned Sequences")
    print("=" * 60)

    test_file = "test_unaligned.fasta"
    create_test_fasta(test_file, aligned=False)

    parser = FastaParser()
    sequences = parser.parse(test_file)
    validation = parser.validate_sequences(sequences)

    print("\nValidation Results:")
    for key, value in validation.items():
        print(f"  {key}: {value}")

    assert validation["valid"] is True
    assert validation["count"] == 4
    assert validation["type"] == "nucleotide"
    assert validation["aligned"] is False
    assert validation["needs_alignment"] is True

    os.remove(test_file)
    print("\n✓ Test passed!")


def test_validation_aligned():
    """Test validation of aligned sequences."""
    print("\n" + "=" * 60)
    print("TEST: Validation - Aligned Sequences")
    print("=" * 60)

    test_file = "test_aligned.fasta"
    create_test_fasta(test_file, aligned=True)

    parser = FastaParser()
    sequences = parser.parse(test_file)
    validation = parser.validate_sequences(sequences)

    print("\nValidation Results:")
    for key, value in validation.items():
        print(f"  {key}: {value}")

    assert validation["valid"] is True
    assert validation["count"] == 4
    assert validation["type"] == "nucleotide"
    assert validation["aligned"] is True
    assert validation["needs_alignment"] is False

    os.remove(test_file)
    print("\n✓ Test passed!")


def test_sequence_properties():
    """Test Sequence class properties."""
    print("\n" + "=" * 60)
    print("TEST: Sequence Properties")
    print("=" * 60)

    seq = Sequence(
        id="test1",
        description="Test sequence with gaps",
        sequence="ATGC-ATGC--ATGC"
    )

    print(f"\nSequence: {seq.sequence}")
    print(f"  Aligned length: {seq.aligned_length}")
    print(f"  Actual length (no gaps): {seq.length}")
    print(f"  Is nucleotide: {seq.is_nucleotide()}")
    print(f"  Is protein: {seq.is_protein()}")

    assert seq.aligned_length == 15
    assert seq.length == 12
    assert seq.is_nucleotide() is True
    assert seq.is_protein() is False

    print("\n✓ Test passed!")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("FASTA PARSER TESTS")
    print("=" * 60)

    test_parser()
    test_validation_unaligned()
    test_validation_aligned()
    test_sequence_properties()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✓")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

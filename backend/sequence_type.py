"""
Sequence type detection and model selection.

This module handles the three different molecular sequence types:
1. DNA (deoxyribonucleic acid) - uses A, C, G, T
2. RNA (ribonucleic acid) - uses A, C, G, U
3. Protein (amino acid sequences) - uses 20 amino acids

Each type requires different substitution models and matrices.
"""

from enum import Enum
from typing import List
from fasta_parser import Sequence


class SequenceType(Enum):
    """Enumeration of supported sequence types."""
    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    UNKNOWN = "unknown"


class SequenceTypeDetector:
    """
    Detect sequence type and validate consistency.

    Critical for phylogenetics:
    - DNA uses 4x4 GTR model
    - RNA uses 4x4 GTR model (U treated as T internally)
    - Protein uses 20x20 empirical models (WAG, LG, JTT)

    Cannot mix types in one analysis!
    """

    # Character sets for each type
    DNA_CHARS = set('ACGT')
    RNA_CHARS = set('ACGU')
    PROTEIN_CHARS = set('ACDEFGHIKLMNPQRSTVWY')

    # Ambiguity codes (IUPAC)
    DNA_AMBIG = set('NRYWSMKHBVD')
    RNA_AMBIG = set('NRYWSMKHBVD')
    PROTEIN_AMBIG = set('XBZJ')

    # Gap characters
    GAP_CHARS = set('-.')

    def __init__(self):
        """Initialize detector."""
        pass

    def detect_single(self, sequence: Sequence) -> SequenceType:
        """
        Detect type of a single sequence.

        Args:
            sequence: Sequence to analyze

        Returns:
            SequenceType enum value
        """
        # Get unique characters (uppercase, no gaps)
        seq_chars = set(sequence.sequence.upper()) - self.GAP_CHARS

        if not seq_chars:
            return SequenceType.UNKNOWN

        # Check for protein-specific characters
        protein_specific = seq_chars & (self.PROTEIN_CHARS - self.DNA_CHARS - self.RNA_CHARS)
        if protein_specific:
            return SequenceType.PROTEIN

        # Check for RNA (has U but no T)
        has_u = 'U' in seq_chars
        has_t = 'T' in seq_chars

        if has_u and not has_t:
            return SequenceType.RNA
        elif has_t and not has_u:
            return SequenceType.DNA
        elif has_u and has_t:
            # Mixed U and T - probably an error, but treat as DNA
            return SequenceType.DNA

        # Only A, C, G (no T or U)
        # Could be either DNA or RNA - check length heuristic
        # DNA/RNA sequences are typically longer
        if sequence.length > 50:
            return SequenceType.DNA  # Default to DNA
        else:
            # Very short, might be peptide
            return SequenceType.PROTEIN if sequence.length < 20 else SequenceType.DNA

    def detect_alignment(self, sequences: List[Sequence]) -> SequenceType:
        """
        Detect type of an alignment (must be consistent).

        Args:
            sequences: List of aligned sequences

        Returns:
            SequenceType (or UNKNOWN if inconsistent)

        Raises:
            ValueError: If sequences have mixed types
        """
        if not sequences:
            return SequenceType.UNKNOWN

        # Detect type of each sequence
        types = [self.detect_single(seq) for seq in sequences]

        # Remove unknowns
        types = [t for t in types if t != SequenceType.UNKNOWN]

        if not types:
            return SequenceType.UNKNOWN

        # Check consistency
        unique_types = set(types)

        if len(unique_types) > 1:
            raise ValueError(
                f"Mixed sequence types detected: {unique_types}. "
                f"All sequences in an alignment must be the same type."
            )

        return types[0]

    def validate_for_phylogenetics(self, sequences: List[Sequence]) -> tuple[SequenceType, str]:
        """
        Validate sequences are suitable for phylogenetic analysis.

        Args:
            sequences: Sequences to validate

        Returns:
            (sequence_type, recommendation)

        Raises:
            ValueError: If sequences are invalid
        """
        seq_type = self.detect_alignment(sequences)

        if seq_type == SequenceType.UNKNOWN:
            raise ValueError("Cannot determine sequence type. Check your input.")

        # Generate recommendation
        if seq_type == SequenceType.DNA:
            recommendation = (
                "DNA sequences detected. Will use GTR+Gamma model.\n"
                "Tip: For rRNA phylogenetics, use rRNA sequences directly (RNA type)."
            )
        elif seq_type == SequenceType.RNA:
            recommendation = (
                "RNA sequences detected. Will use GTR+Gamma model (U→T mapping).\n"
                "Perfect for rRNA phylogenetics (16S, 23S, 18S, etc.)."
            )
        elif seq_type == SequenceType.PROTEIN:
            recommendation = (
                "Protein sequences detected. Will use WAG or LG model.\n"
                "Note: Protein phylogenetics requires different models than DNA/RNA."
            )
        else:
            recommendation = "Unknown sequence type."

        return seq_type, recommendation


def get_appropriate_model(seq_type: SequenceType) -> str:
    """
    Get the appropriate substitution model for a sequence type.

    Args:
        seq_type: Type of sequences

    Returns:
        Model name as string
    """
    model_map = {
        SequenceType.DNA: "GTR+G",
        SequenceType.RNA: "GTR+G",  # Same as DNA, U→T internally
        SequenceType.PROTEIN: "WAG",  # Default protein model
    }

    return model_map.get(seq_type, "GTR+G")


# Example usage
if __name__ == "__main__":
    detector = SequenceTypeDetector()

    # Test cases
    test_sequences = [
        Sequence("dna1", "DNA test", "ATGCATGC"),
        Sequence("rna1", "RNA test", "AUGCAUGC"),
        Sequence("protein1", "Protein test", "MKTAYIAKQRQISFVKSHFSRQ"),
    ]

    print("Sequence Type Detection Tests")
    print("=" * 60)

    for seq in test_sequences:
        seq_type = detector.detect_single(seq)
        print(f"{seq.id:15} {seq.sequence[:20]:20} → {seq_type.value}")

    print("\nTesting alignment validation:")
    print("-" * 60)

    # DNA alignment (should work)
    dna_seqs = [
        Sequence("seq1", "Test 1", "ATGCAT"),
        Sequence("seq2", "Test 2", "ATGCAC"),
    ]
    seq_type, rec = detector.validate_for_phylogenetics(dna_seqs)
    print(f"DNA alignment: {seq_type.value}")
    print(f"Recommendation: {rec}\n")

    # RNA alignment (should work)
    rna_seqs = [
        Sequence("rna1", "16S rRNA", "AUGCAUGC"),
        Sequence("rna2", "16S rRNA", "AUGCAUGC"),
    ]
    seq_type, rec = detector.validate_for_phylogenetics(rna_seqs)
    print(f"RNA alignment: {seq_type.value}")
    print(f"Recommendation: {rec}\n")

    # Mixed alignment (should fail)
    mixed_seqs = [
        Sequence("seq1", "DNA", "ATGCAT"),
        Sequence("seq2", "RNA", "AUGCAU"),
    ]
    try:
        seq_type, rec = detector.validate_for_phylogenetics(mixed_seqs)
    except ValueError as e:
        print(f"Mixed alignment error (expected): {e}")

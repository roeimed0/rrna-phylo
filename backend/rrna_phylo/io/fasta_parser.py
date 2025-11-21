"""
Simple FASTA parser for phylogenetic analysis.

Reads FASTA files, validates sequences, returns structured data.
"""

from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class Sequence:
    """Represents a single biological sequence."""
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        """Return sequence length (without gaps)."""
        return len(self.sequence.replace('-', '').replace('.', ''))

    @property
    def aligned_length(self) -> int:
        """Return full length including gaps."""
        return len(self.sequence)

    def is_nucleotide(self) -> bool:
        """Check if sequence is nucleotide (DNA/RNA)."""
        nucleotide_chars = set('ACGTUN-.')
        seq_chars = set(self.sequence.upper())
        return seq_chars.issubset(nucleotide_chars)

    def is_protein(self) -> bool:
        """
        Check if sequence is protein (amino acid).

        Returns True only if it contains protein-specific characters
        not found in nucleotides (E, F, I, L, P, Q, etc.)
        """
        nucleotide_only = set('ACGTUN-.')
        protein_specific = set('EFHIKLMPQRSVWY*')
        seq_chars = set(self.sequence.upper())

        # If it has protein-specific chars, it's definitely protein
        if seq_chars.intersection(protein_specific):
            return True

        # If it only has nucleotide chars, it's nucleotide not protein
        if seq_chars.issubset(nucleotide_only):
            return False

        # Has D or other ambiguous chars - could be either
        # Default to protein if it passes protein validation
        protein_chars = set('ACDEFGHIKLMNPQRSTVWY*-.')
        return seq_chars.issubset(protein_chars)


class FastaParser:
    """Parse FASTA format files."""

    def parse(self, filepath: str) -> List[Sequence]:
        """
        Parse FASTA file and return list of sequences.

        Args:
            filepath: Path to FASTA file

        Returns:
            List of Sequence objects

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file is not valid FASTA format
        """
        sequences = []
        current_id = None
        current_desc = None
        current_seq = []

        try:
            with open(filepath, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines
                    if not line:
                        continue

                    # Header line
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_id is not None:
                            seq = ''.join(current_seq)
                            sequences.append(Sequence(
                                id=current_id,
                                description=current_desc,
                                sequence=seq
                            ))

                        # Parse new header
                        header = line[1:].strip()
                        parts = header.split(None, 1)  # Split on first whitespace
                        current_id = parts[0]
                        current_desc = parts[1] if len(parts) > 1 else ""
                        current_seq = []

                    # Sequence line
                    else:
                        if current_id is None:
                            raise ValueError(
                                f"Line {line_num}: Sequence data before header"
                            )
                        current_seq.append(line)

                # Don't forget last sequence
                if current_id is not None:
                    seq = ''.join(current_seq)
                    sequences.append(Sequence(
                        id=current_id,
                        description=current_desc,
                        sequence=seq
                    ))

        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filepath}")

        if not sequences:
            raise ValueError("No sequences found in file")

        return sequences

    def validate_sequences(self, sequences: List[Sequence]) -> Dict[str, any]:
        """
        Validate a list of sequences for phylogenetic analysis.

        Args:
            sequences: List of Sequence objects

        Returns:
            Dictionary with validation results
        """
        if not sequences:
            return {
                'valid': False,
                'error': 'No sequences provided'
            }

        if len(sequences) < 3:
            return {
                'valid': False,
                'error': 'Need at least 3 sequences for phylogenetic analysis'
            }

        # Check if all sequences are same type (nucleotide or protein)
        first_seq = sequences[0]
        is_nucleotide = first_seq.is_nucleotide()
        is_protein = first_seq.is_protein()

        mixed_types = False
        for seq in sequences[1:]:
            if is_nucleotide and not seq.is_nucleotide():
                mixed_types = True
                break
            if is_protein and not seq.is_protein():
                mixed_types = True
                break

        if mixed_types:
            return {
                'valid': False,
                'error': 'Mixed sequence types (nucleotide and protein)'
            }

        # Check for duplicate IDs
        ids = [seq.id for seq in sequences]
        if len(ids) != len(set(ids)):
            return {
                'valid': False,
                'error': 'Duplicate sequence IDs found'
            }

        # Check if sequences are aligned (same length)
        lengths = [seq.aligned_length for seq in sequences]
        all_same_length = len(set(lengths)) == 1

        seq_type = 'nucleotide' if is_nucleotide else 'protein'

        return {
            'valid': True,
            'count': len(sequences),
            'type': seq_type,
            'aligned': all_same_length,
            'lengths': lengths,
            'min_length': min(lengths),
            'max_length': max(lengths),
            'needs_alignment': not all_same_length
        }


def parse_fasta(filepath: str) -> List[Sequence]:
    """
    Convenience function to parse FASTA file.

    Args:
        filepath: Path to FASTA file

    Returns:
        List of Sequence objects
    """
    parser = FastaParser()
    return parser.parse(filepath)

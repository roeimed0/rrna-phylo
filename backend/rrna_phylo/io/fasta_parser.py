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
    _unique_display_name: str = None  # Set by assign_unique_display_names()

    @property
    def length(self) -> int:
        """Return sequence length (without gaps)."""
        return len(self.sequence.replace('-', '').replace('.', ''))

    @property
    def aligned_length(self) -> int:
        """Return full length including gaps."""
        return len(self.sequence)

    @property
    def species_name(self) -> str:
        """
        Extract species name from description.

        For taxonomy descriptions like "Bacteria;...;Genus;Species strain",
        returns the last part (species + strain).

        Returns:
            Species name or empty string if not found
        """
        if not self.description:
            return ""

        # Handle taxonomy format (semicolon-separated)
        if ';' in self.description:
            parts = self.description.split(';')
            # Last part is usually the species/strain
            return parts[-1].strip()

        # Handle space-separated format
        # Take first two words as genus + species
        words = self.description.split()
        if len(words) >= 2:
            return ' '.join(words[:2])
        elif len(words) == 1:
            return words[0]

        return ""

    @property
    def main_accession(self) -> str:
        """
        Extract main accession number from ID.

        For IDs like "AE006468.4100145.4101688", returns "AE006468"
        For IDs like "U00096.223771.225312", returns "U00096"

        Returns:
            Main accession number (first part before dot)
        """
        if '.' in self.id:
            return self.id.split('.')[0]
        return self.id

    @property
    def display_name(self) -> str:
        """
        Get display name for tree visualization.

        Returns species name with main accession in parentheses,
        e.g., "Bacillus subtilis (AJ276351)"

        If assign_unique_display_names() was called, returns the unique name
        with counter (e.g., "Escherichia coli (U00096) #1")

        Falls back to ID if no species name available.
        """
        # Return unique name if it was assigned
        if self._unique_display_name:
            return self._unique_display_name

        # Use species name with main accession in parentheses
        # NOTE: Do NOT add quotes here - they will be added by to_newick() if needed
        species = self.species_name
        accession = self.main_accession

        if species:
            return f"{species} ({accession})"
        else:
            # No species name, just use the accession
            return accession

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


def assign_unique_display_names(sequences: List[Sequence]) -> List[Sequence]:
    """
    Assign unique display names to sequences, adding counters for duplicates.

    When multiple sequences have the same species name, adds #1, #2, etc.
    to differentiate them in tree visualizations.

    Args:
        sequences: List of Sequence objects

    Returns:
        Same list of sequences (modified in place)

    Example:
        >>> seqs = parse_fasta("sequences.fasta")
        >>> assign_unique_display_names(seqs)
        # Now seqs will have unique display names like:
        # "Escherichia coli (U00096) #1"
        # "Escherichia coli (U00096) #2"
    """
    # Count occurrences of each species name
    species_counts = {}
    species_current = {}

    for seq in sequences:
        base_name = seq.display_name
        species_counts[base_name] = species_counts.get(base_name, 0) + 1

    # Add counters only if there are duplicates
    for seq in sequences:
        base_name = seq.display_name

        if species_counts[base_name] > 1:
            # Multiple sequences with same species name
            counter = species_current.get(base_name, 0) + 1
            species_current[base_name] = counter
            # Override the display name with a counter
            seq._unique_display_name = f"{base_name} #{counter}"
        else:
            # Unique species name, no counter needed
            seq._unique_display_name = base_name

    return sequences

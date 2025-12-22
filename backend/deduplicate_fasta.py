#!/usr/bin/env python3
"""
Deduplicate FASTA file by keeping only the best sequence per species.

Strategy:
1. Parse all sequences
2. Group by species name (last field in header)
3. For each species, keep the longest sequence
4. Write deduplicated output

Usage:
    python deduplicate_fasta.py input.fasta output.fasta
"""

import sys
from pathlib import Path
from collections import defaultdict


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return dict: {species: [(header, sequence), ...]}
    """
    sequences_by_species = defaultdict(list)

    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    seq = ''.join(current_seq)
                    # Extract species name (last field after |)
                    species = current_header.split('|')[-1]
                    sequences_by_species[species].append((current_header, seq))

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header is not None:
            seq = ''.join(current_seq)
            species = current_header.split('|')[-1]
            sequences_by_species[species].append((current_header, seq))

    return sequences_by_species


def select_best_sequence(sequences):
    """
    Select the best sequence from a list.

    Strategy: Choose the longest sequence.
    """
    if len(sequences) == 1:
        return sequences[0]

    # Find longest sequence
    best = max(sequences, key=lambda x: len(x[1]))
    return best


def deduplicate_fasta(input_file, output_file):
    """
    Deduplicate FASTA file by keeping only best sequence per species.
    """
    print(f"Reading: {input_file}")
    sequences_by_species = parse_fasta(input_file)

    print(f"Found {len(sequences_by_species)} unique species")

    # Count sequences before deduplication
    total_before = sum(len(seqs) for seqs in sequences_by_species.values())
    print(f"Total sequences before: {total_before}")

    # Select best sequence for each species
    deduplicated = {}
    for species, sequences in sequences_by_species.items():
        if len(sequences) > 1:
            print(f"  {species}: {len(sequences)} sequences -> keeping longest ({len(select_best_sequence(sequences)[1])} bp)")
        deduplicated[species] = select_best_sequence(sequences)

    # Write output
    print(f"\nWriting deduplicated file: {output_file}")
    with open(output_file, 'w') as f:
        for species in sorted(deduplicated.keys()):
            header, seq = deduplicated[species]
            f.write(f'>{header}\n')
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f"\nDone!")
    print(f"  Input:  {total_before} sequences")
    print(f"  Output: {len(deduplicated)} sequences (one per species)")
    print(f"  Removed: {total_before - len(deduplicated)} duplicate sequences")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python deduplicate_fasta.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not Path(input_file).exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    deduplicate_fasta(input_file, output_file)

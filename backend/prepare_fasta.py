#!/usr/bin/env python3
"""
Prepare FASTA file for phylogenetic analysis.

This script combines two preprocessing steps:
1. Deduplication - Keep only the longest sequence per species
2. Header cleaning - Convert to standard ID|Species_name format

Usage:
    python prepare_fasta.py input.fasta output.fasta

Example:
    python prepare_fasta.py data/raw_data.fasta data/clean_data.fasta
"""

import sys
from pathlib import Path
from collections import defaultdict
import re


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return list of (header, sequence) tuples.
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def extract_species_name(header):
    """
    Extract species name from FASTA header.

    Tries multiple strategies:
    1. PR2 format: last field after |
    2. NCBI format: second and third words
    3. Use header as-is if simple format

    Args:
        header: FASTA header line (without >)

    Returns:
        Species name with underscores
    """
    # Strategy 1: PR2 format (pipe-separated, species at end)
    if '|' in header:
        parts = header.split('|')
        species = parts[-1].strip()

        # Replace spaces with underscores
        species = species.replace(' ', '_')

        # Remove any trailing description after first space
        species = species.split()[0] if ' ' in species else species

        return species

    # Strategy 2: NCBI format (space-separated)
    # >AB571241.1 Homo sapiens 18S ribosomal RNA
    parts = header.split()
    if len(parts) >= 2:
        # Take second and third words (genus + species)
        if len(parts) >= 3:
            species = f"{parts[1]}_{parts[2]}"
        else:
            species = parts[1]

        # Remove special characters
        species = re.sub(r'[^\w_]', '', species)

        return species

    # Strategy 3: Use as-is (already simple format)
    species = header.strip().replace(' ', '_')
    species = re.sub(r'[^\w_]', '', species)

    return species


def extract_accession(header):
    """
    Extract accession ID from FASTA header.

    Args:
        header: FASTA header line (without >)

    Returns:
        Cleaned accession ID
    """
    if '|' in header:
        # PR2 format: first field is ID
        accession = header.split('|')[0]
    elif ' ' in header:
        # NCBI format: first field is ID
        accession = header.split()[0]
    else:
        # Simple format: use whole header as ID
        accession = header

    # Clean accession (remove special characters after dot)
    accession = accession.split('.')[0] if '.' in accession else accession
    accession = accession.split('_')[0] if '_U' in accession else accession

    return accession


def deduplicate_sequences(sequences):
    """
    Deduplicate sequences by keeping longest sequence per species.

    Args:
        sequences: List of (header, sequence) tuples

    Returns:
        Tuple of (deduplicated_sequences, duplicates_found)
    """
    sequences_by_species = defaultdict(list)

    # Group by species name
    for header, seq in sequences:
        species = extract_species_name(header)
        sequences_by_species[species].append((header, seq))

    # Select best (longest) sequence for each species
    deduplicated = []
    duplicates_found = []

    for species in sorted(sequences_by_species.keys()):
        seqs = sequences_by_species[species]

        if len(seqs) > 1:
            duplicates_found.append((species, len(seqs)))

        # Keep longest sequence
        best = max(seqs, key=lambda x: len(x[1]))
        deduplicated.append((species, best[0], best[1]))

    return deduplicated, duplicates_found


def clean_headers(deduplicated_sequences):
    """
    Convert headers to standard ID|Species_name format.

    Args:
        deduplicated_sequences: List of (species, original_header, sequence) tuples

    Returns:
        List of (cleaned_header, sequence) tuples
    """
    cleaned = []

    for species, original_header, seq in deduplicated_sequences:
        accession = extract_accession(original_header)
        new_header = f"{accession}|{species}"
        cleaned.append((new_header, seq))

    return cleaned


def write_fasta(sequences, output_file):
    """
    Write sequences to FASTA file.

    Args:
        sequences: List of (header, sequence) tuples
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f'>{header}\n')
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def prepare_fasta(input_file, output_file):
    """
    Complete FASTA preparation workflow:
    1. Parse sequences
    2. Deduplicate (keep longest per species)
    3. Clean headers to ID|Species_name format
    4. Write output

    Args:
        input_file: Input FASTA file (raw data)
        output_file: Output FASTA file (prepared data)
    """
    print("=" * 70)
    print("FASTA PREPARATION WORKFLOW")
    print("=" * 70)
    print()

    # Step 1: Parse
    print(f"[1/4] Reading: {input_file}")
    sequences = parse_fasta(input_file)
    print(f"      Loaded {len(sequences)} sequences")
    print()

    # Step 2: Deduplicate
    print("[2/4] Deduplicating sequences (keeping longest per species)...")
    deduplicated, duplicates_found = deduplicate_sequences(sequences)

    if duplicates_found:
        print(f"      Found duplicates for {len(duplicates_found)} species:")
        for species, count in duplicates_found[:10]:  # Show first 10
            best_seq = next(seq for sp, _, seq in deduplicated if sp == species)
            print(f"        {species}: {count} sequences -> kept longest ({len(best_seq)} bp)")
        if len(duplicates_found) > 10:
            print(f"        ... and {len(duplicates_found) - 10} more")
    else:
        print("      No duplicates found (all species are unique)")

    print(f"      Result: {len(sequences)} -> {len(deduplicated)} sequences")
    print()

    # Step 3: Clean headers
    print("[3/4] Cleaning headers to standard format (ID|Species_name)...")
    cleaned = clean_headers(deduplicated)
    print(f"      Converted {len(cleaned)} headers")
    print()
    print("      Example transformations:")
    for i in range(min(3, len(cleaned))):
        original = deduplicated[i][1]
        new = cleaned[i][0]
        # Show shortened version if too long
        if len(original) > 50:
            original = original[:47] + "..."
        print(f"        Before: >{original}")
        print(f"        After:  >{new}")
        print()

    # Step 4: Write output
    print(f"[4/4] Writing prepared data: {output_file}")
    write_fasta(cleaned, output_file)

    output_size = Path(output_file).stat().st_size / 1024
    print(f"      File size: {output_size:.1f} KB")
    print()

    # Summary
    print("=" * 70)
    print("PREPARATION COMPLETE!")
    print("=" * 70)
    print()
    print(f"Input:  {len(sequences)} sequences")
    print(f"Output: {len(cleaned)} unique species")
    if duplicates_found:
        print(f"Removed: {len(sequences) - len(cleaned)} duplicate sequences")
    print()
    print("Standard format: >ID|Species_name")
    print("  - ID: Traceable accession number")
    print("  - Species_name: Readable name for phylogenetic trees")
    print()
    print(f"Ready for tree building: {output_file}")
    print("=" * 70)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("=" * 70)
        print("FASTA Preparation Tool")
        print("=" * 70)
        print()
        print("Prepares FASTA files for phylogenetic analysis by:")
        print("  1. Deduplicating sequences (keeps longest per species)")
        print("  2. Cleaning headers to standard ID|Species_name format")
        print()
        print("USAGE:")
        print("  python prepare_fasta.py input.fasta output.fasta")
        print()
        print("EXAMPLE:")
        print("  python prepare_fasta.py data/raw_pr2.fasta data/clean_pr2.fasta")
        print()
        print("STANDARD FORMAT:")
        print("  >AB571241|Homo_sapiens")
        print("  >NC_012920|Callithrix_jacchus")
        print()
        print("Benefits:")
        print("  - Removes duplicate sequences (multiple per species)")
        print("  - Keeps longest/best sequence for each species")
        print("  - Standardizes headers for readable phylogenetic trees")
        print("  - One command instead of two separate tools")
        print("=" * 70)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not Path(input_file).exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    prepare_fasta(input_file, output_file)

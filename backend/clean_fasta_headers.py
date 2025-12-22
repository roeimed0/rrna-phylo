#!/usr/bin/env python3
"""
Clean FASTA headers to standard format: ID|Species_name

STANDARD FORMAT for rRNA-Phylo project:
    >ID|Species_name

Example:
    >AB571241|Homo_sapiens
    >NC_012920|Callithrix_jacchus

This format provides:
- Traceability: Keep original ID (accession, genome ID, etc.)
- Readability: Species name shown in phylogenetic trees
- Consistency: All trees use the same labeling scheme

Supports conversion from multiple input formats:
- PR2 format: >accession|taxonomy|...|Species_name
- NCBI format: >accession Species name description
- Simple format: >Species_name (ID assigned as species name)

Usage:
    python clean_fasta_headers.py input.fasta output.fasta
"""

import sys
from pathlib import Path
import re


def extract_species_name(header):
    """
    Extract species name from FASTA header.

    Tries multiple strategies:
    1. PR2 format: last field after |
    2. NCBI format: second field
    3. Use header as-is if simple format

    Args:
        header: FASTA header line (without >)

    Returns:
        Species name
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


def clean_fasta_headers(input_file, output_file, standard_format=True):
    """
    Clean FASTA headers to standard ID|Species_name format.

    Args:
        input_file: Input FASTA file
        output_file: Output FASTA file with cleaned headers
        standard_format: If True, use ID|Species_name format (default, recommended)
                        If False, use Species_name only
    """
    print(f"Reading: {input_file}")

    sequences = []
    current_header = None
    current_seq = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))

                # Process new header
                original_header = line[1:]  # Remove >
                species_name = extract_species_name(original_header)

                # Extract ID (accession number)
                if '|' in original_header:
                    # PR2 format: first field is ID
                    accession = original_header.split('|')[0]
                elif ' ' in original_header:
                    # NCBI format: first field is ID
                    accession = original_header.split()[0]
                else:
                    # Simple format: use species name as ID
                    accession = species_name

                # Clean accession (remove special characters after dot)
                accession = accession.split('.')[0] if '.' in accession else accession
                accession = accession.split('_')[0] if '_U' in accession else accession

                if standard_format:
                    # Standard format: ID|Species_name
                    new_header = f"{accession}|{species_name}"
                else:
                    # Species name only
                    new_header = species_name

                current_header = new_header
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    print(f"Processed {len(sequences)} sequences")

    # Check for duplicates
    species_counts = {}
    for header, seq in sequences:
        name = header.split('|')[-1] if '|' in header else header
        species_counts[name] = species_counts.get(name, 0) + 1

    duplicates = {k: v for k, v in species_counts.items() if v > 1}
    if duplicates:
        print(f"\nWarning: Found duplicate species names:")
        for name, count in duplicates.items():
            print(f"  {name}: {count} sequences")
        print("\nConsider running deduplicate_fasta.py first!")

    # Write output
    print(f"\nWriting: {output_file}")
    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f'>{header}\n')
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f"\nDone! Headers cleaned to standard format.")
    print(f"\nExample transformation:")
    print(f"  Before: >AB571241.1.1364_U|Obazoa|...|Callithrix_jacchus")
    print(f"  After:  >{sequences[0][0]}")
    print(f"\nStandard format: ID|Species_name")
    print(f"  - ID: Traceable accession/identifier")
    print(f"  - Species_name: Readable name shown in phylogenetic trees")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("=" * 70)
        print("FASTA Header Cleaning Tool")
        print("=" * 70)
        print()
        print("Converts FASTA headers to standard format: ID|Species_name")
        print()
        print("USAGE:")
        print("  python clean_fasta_headers.py input.fasta output.fasta [OPTIONS]")
        print()
        print("OPTIONS:")
        print("  --species-only      Use species name only (no ID)")
        print("                      Default: ID|Species_name (recommended)")
        print()
        print("EXAMPLES:")
        print("  # Standard format (ID|Species_name)")
        print("  python clean_fasta_headers.py pr2_data.fasta cleaned.fasta")
        print()
        print("  # Species name only")
        print("  python clean_fasta_headers.py pr2_data.fasta cleaned.fasta --species-only")
        print()
        print("STANDARD FORMAT:")
        print("  >AB571241|Homo_sapiens")
        print("  >NC_012920|Callithrix_jacchus")
        print()
        print("  Benefits:")
        print("    - Traceability: Keep original ID for reference")
        print("    - Readability: Species name shown in phylogenetic trees")
        print("    - Consistency: All FASTA files use same format")
        print("=" * 70)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    species_only = '--species-only' in sys.argv
    standard_format = not species_only  # Default to standard format

    if not Path(input_file).exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    clean_fasta_headers(input_file, output_file, standard_format)

#!/usr/bin/env python3
"""
Test script to verify output file organization.

This creates a minimal test dataset and runs build_trees.py to verify
that output files are organized in subfolders correctly.
"""

import sys
from pathlib import Path


def create_test_fasta():
    """Create a test FASTA file with 50 sequences of 1000bp each."""
    test_file = Path('test_sequences.fasta')

    # Base sequence: 1000bp repeating pattern with variations
    base_seq = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT' * 25  # 1000bp

    sequences = []
    for i in range(50):
        # Create variations by substituting bases at different positions
        seq = list(base_seq)

        # Add unique mutations for each species
        for j in range(i * 3, min(i * 3 + 10, 1000)):
            # Rotate bases to create diversity
            bases = ['A', 'C', 'G', 'T']
            current_base = seq[j]
            if current_base in bases:
                seq[j] = bases[(bases.index(current_base) + i) % 4]

        seq_str = ''.join(seq)
        species_name = f'Species_{i+1:03d}'
        sequences.append((species_name, seq_str))

    with open(test_file, 'w') as f:
        for name, seq in sequences:
            f.write(f'>{name}\n{seq}\n')

    print(f'Created test file: {test_file}')
    print(f'  Sequences: 50')
    print(f'  Length: 1000 bp each')
    return test_file


def verify_output_structure(input_basename):
    """Verify that output files are in the correct location."""
    expected_dir = Path('results') / input_basename

    expected_files = [
        'upgma_tree.nwk',
        'bionj_tree.nwk',
        'ml_tree.nwk',
        'upgma_tree.txt',
        'bionj_tree.txt',
        'ml_tree.txt',
        'summary.txt'
    ]

    print(f'\nVerifying output structure in: {expected_dir}')
    print('=' * 70)

    if not expected_dir.exists():
        print(f'[FAIL] ERROR: Output directory not found: {expected_dir}')
        return False

    print(f'[OK] Output directory exists: {expected_dir}')
    print()

    all_good = True
    for filename in expected_files:
        filepath = expected_dir / filename
        if filepath.exists():
            size = filepath.stat().st_size
            print(f'[OK] {filename:<20} ({size:>6} bytes)')
        else:
            print(f'[FAIL] {filename:<20} (MISSING)')
            all_good = False

    print()
    if all_good:
        print('[OK] All expected files created successfully!')
        print(f'[OK] Output organized in: {expected_dir.absolute()}')
    else:
        print('[FAIL] Some files are missing!')

    return all_good


def main():
    """Run the test."""
    print('=' * 70)
    print('TESTING OUTPUT FILE ORGANIZATION')
    print('=' * 70)
    print()

    # Create test data
    test_file = create_test_fasta()
    input_basename = test_file.stem  # 'test_sequences'

    print()
    print('Now run:')
    print(f'  python build_trees.py {test_file}')
    print()
    print('After the trees are built, run:')
    print('  python test_output_organization.py --verify')
    print()
    print('To clean up test files, run:')
    print('  bash cleanup_test.sh')
    print()

    # If --verify flag is passed, verify the output
    if '--verify' in sys.argv:
        print()
        success = verify_output_structure(input_basename)
        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()

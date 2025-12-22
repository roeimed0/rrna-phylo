#!/bin/bash
#
# Cleanup script for test output organization
#
# This removes all test files created by test_output_organization.py
#

echo "=========================================="
echo "CLEANING UP TEST FILES"
echo "=========================================="
echo

# Remove test FASTA file
if [ -f "test_sequences.fasta" ]; then
    echo "Removing: test_sequences.fasta"
    rm test_sequences.fasta
fi

# Remove test output directory
if [ -d "results/test_sequences" ]; then
    echo "Removing: results/test_sequences/"
    rm -rf results/test_sequences
fi

# Check if results/ directory is empty, if so remove it
if [ -d "results" ]; then
    if [ -z "$(ls -A results)" ]; then
        echo "Removing: results/ (empty)"
        rmdir results
    else
        echo "Keeping: results/ (contains other files)"
    fi
fi

echo
echo "✓ Cleanup complete!"
echo

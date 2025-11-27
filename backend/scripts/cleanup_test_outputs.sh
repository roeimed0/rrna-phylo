#!/bin/bash
# Clean up test output files and directories
# Run this before starting new tests to avoid confusion

echo "Cleaning up test output files..."

cd "$(dirname "$0")/.." || exit 1

# Remove test output text files
rm -f test_*performance*.txt
rm -f test_*output.txt
rm -f test_smart*.txt

# Remove test output directories
rm -rf test_*_output/

# Remove temporary aligned FASTA files (keep the original datasets)
rm -f test_*aligned.fasta

echo "Cleanup complete!"
echo ""
echo "Files kept:"
ls -lh test_real_*.fasta 2>/dev/null || echo "  (no dataset files found)"

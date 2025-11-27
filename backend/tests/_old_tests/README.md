# Old/Archived Tests

This folder contains older test files that have been superseded by more comprehensive tests.

## Archived Tests

### Individual Method Tests (Superseded by `test_phylo_builder.py` and `test_full_pipeline.py`)
- `test_upgma.py` - UPGMA tree building (now tested in `test_phylo_builder.py`)
- `test_bionj.py` - BioNJ tree building (now tested in `test_phylo_builder.py`)
- `test_distance.py` - Distance matrix calculations (now tested in `test_phylo_builder.py`)
- `test_ml_tree.py` - Basic ML tree (superseded by `test_ml_level3.py`)
- `test_ml_level2.py` - Intermediate ML (superseded by `test_ml_level3.py`)

## Why Archived?

These tests were created during early development to test individual methods. They have been superseded by:

1. **test_phylo_builder.py** - Comprehensive tests for all tree building methods (UPGMA, BioNJ, ML)
2. **test_full_pipeline.py** - End-to-end pipeline testing
3. **test_ml_level3.py** - Advanced ML testing with full features
4. **test_large_dataset_performance.py** - Large-scale real-world testing

## Active Tests

The following tests in `backend/tests/` are actively maintained:

- `test_aligner.py` - MUSCLE alignment
- `test_parser.py` - FASTA parsing
- `test_sequence_type.py` - Sequence type detection
- `test_strain_handler.py` - Deduplication
- `test_phylo_builder.py` - All tree building methods
- `test_ml_level3.py` - Advanced ML features
- `test_protein_phylo.py` - Protein phylogenetics
- `test_full_pipeline.py` - Complete pipeline
- `test_large_dataset_performance.py` - Large-scale performance testing

## Notes

These old tests are kept for reference but are not run in the main test suite. They may be deleted in the future or used for reference when debugging specific features.

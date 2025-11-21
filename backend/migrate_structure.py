"""
Script to reorganize backend codebase into proper package structure.

This migrates from flat structure to organized modules:
backend/
├── rrna_phylo/     # Main package
└── tests/          # Separate test directory
"""

import os
import shutil
from pathlib import Path

# Get backend directory
BACKEND_DIR = Path(__file__).parent
PACKAGE_DIR = BACKEND_DIR / "rrna_phylo"

# File mappings: source -> destination
FILE_MOVES = {
    # Already done:
    # "fasta_parser.py": "rrna_phylo/io/fasta_parser.py",
    # "aligner.py": "rrna_phylo/io/aligner.py",
    # "sequence_type.py": "rrna_phylo/core/sequence_type.py",

    # Distance modules
    "distance.py": "rrna_phylo/distance/nucleotide.py",
    "protein_distance.py": "rrna_phylo/distance/protein.py",

    # Model modules
    "ml_tree.py": "rrna_phylo/models/dna_models.py",  # GTR model
    "protein_models.py": "rrna_phylo/models/protein_models.py",

    # Method modules (tree building)
    "upgma.py": "rrna_phylo/methods/upgma.py",
    "bionj.py": "rrna_phylo/methods/bionj.py",
    "ml_tree_level2.py": "rrna_phylo/methods/_ml_level2.py",  # Temporary
    "ml_tree_level3.py": "rrna_phylo/methods/ml_nucleotide.py",  # Main ML
    "protein_ml.py": "rrna_phylo/methods/ml_protein.py",

    # Core builder
    "phylo_builder.py": "rrna_phylo/core/builder.py",

    # Utils
    "visualize_trees.py": "rrna_phylo/utils/visualization.py",
    "compare_trees.py": "rrna_phylo/utils/comparison.py",

    # Tests
    "test_parser.py": "tests/test_io.py",
    "test_aligner.py": "tests/test_io_aligner.py",
    "test_distance.py": "tests/test_distance.py",
    "test_upgma.py": "tests/test_upgma.py",
    "test_bionj.py": "tests/test_bionj.py",
    "test_ml_tree.py": "tests/_test_ml_level1.py",  # Archive
    "test_ml_level2.py": "tests/_test_ml_level2.py",  # Archive
    "test_ml_level3.py": "tests/test_ml_nucleotide.py",
    "test_protein_phylo.py": "tests/test_protein.py",
    "test_sequence_type.py": "tests/test_sequence_type.py",
    "test_phylo_builder.py": "tests/test_builder.py",
}

def main():
    print("=" * 70)
    print("REORGANIZING BACKEND STRUCTURE")
    print("=" * 70)

    # Show plan
    print("\nPlan:")
    print("-" * 70)
    for src, dst in FILE_MOVES.items():
        if (BACKEND_DIR / src).exists():
            print(f"  {src:30} → {dst}")
        else:
            print(f"  {src:30} → {dst} (SKIP - not found)")

    # Ask confirmation
    response = input("\nProceed with migration? (yes/no): ")
    if response.lower() != "yes":
        print("Aborted.")
        return

    # Do migration
    print("\nMigrating files...")
    for src, dst in FILE_MOVES.items():
        src_path = BACKEND_DIR / src
        dst_path = BACKEND_DIR / dst

        if not src_path.exists():
            print(f"  ⊗ SKIP {src} (not found)")
            continue

        # Create destination directory if needed
        dst_path.parent.mkdir(parents=True, exist_ok=True)

        # Copy file
        shutil.copy2(src_path, dst_path)
        print(f"  ✓ {src} → {dst}")

    print("\n" + "=" * 70)
    print("MIGRATION COMPLETE!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Update imports in moved files")
    print("  2. Create __init__.py files")
    print("  3. Test everything")
    print("  4. Remove old files if all tests pass")

if __name__ == "__main__":
    main()

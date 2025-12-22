#!/usr/bin/env python3
"""
rRNA-Phylo Application Launcher

Main entry point for the rRNA-Phylo phylogenetic tree builder.
Provides both CLI and interactive menu interfaces.

Usage:
    # Interactive menu (recommended for first-time users)
    python app.py

    # Direct CLI (for advanced users and scripting)
    python app.py build sequences.fasta
    python app.py build sequences.fasta --method all --visualize pdf
    python app.py test
    python app.py clean

For help:
    python app.py --help
"""

import os
import sys
from pathlib import Path

# Prevent Intel OpenMP library conflicts
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'


def main():
    """Main application entry point."""

    # If no arguments, launch interactive menu
    if len(sys.argv) == 1:
        from rrna_phylo_app import main_menu
        main_menu()
        return

    # Otherwise, handle CLI commands
    command = sys.argv[1].lower()

    if command in ['--help', '-h', 'help']:
        show_help()
    elif command == 'build':
        run_build()
    elif command == 'test':
        run_test()
    elif command == 'clean':
        run_clean()
    elif command == 'menu':
        from rrna_phylo_app import main_menu
        main_menu()
    else:
        print(f"Unknown command: {command}")
        print("Run 'python app.py --help' for usage information")
        sys.exit(1)


def show_help():
    """Show help information."""
    print("""
rRNA-Phylo: Multi-Method Phylogenetic Tree Builder

USAGE:
    python app.py                           Launch interactive menu
    python app.py build <file> [options]    Build phylogenetic trees
    python app.py test                      Create and test with sample data
    python app.py clean                     Clean up test files
    python app.py menu                      Launch interactive menu
    python app.py --help                    Show this help

INTERACTIVE MENU:
    python app.py

    Provides guided interface with:
    • Quick build (all 3 methods)
    • Custom build (choose options)
    • Advanced ML (model selection)
    • Test data creation
    • Results viewer
    • Cleanup tool
    • ETE3 installer
    • Built-in help

BUILD COMMAND:
    python app.py build <fasta_file> [options]

    Options:
        --method {all,upgma,bionj,ml}       Tree method (default: all)
        --bootstrap N                       Bootstrap replicates (default: 0)
        --visualize {pdf,png,svg}          Visualization format
        --dpi N                            DPI for PNG (default: 300)
        --output DIR                        Output directory (default: results)

    Examples:
        python app.py build sequences.fasta
        python app.py build sequences.fasta --method ml --bootstrap 100
        python app.py build sequences.fasta --visualize pdf

TEST COMMAND:
    python app.py test

    Creates test data (50 sequences, 1000bp), builds all 3 trees,
    verifies output structure, and reports results.

CLEAN COMMAND:
    python app.py clean

    Removes test_sequences.fasta and results/ directory.
    Use with caution!

DOCUMENTATION:
    See README.md for comprehensive documentation.

WEBSITE:
    https://github.com/yourusername/rrna-phylo
""")


def run_build():
    """Run tree building with CLI arguments."""
    # Remove 'build' command from args and pass rest to build_trees.py
    import subprocess

    if len(sys.argv) < 3:
        print("Error: No input file specified")
        print("Usage: python app.py build <fasta_file> [options]")
        sys.exit(1)

    # Get input file
    input_file = sys.argv[2]

    if not Path(input_file).exists():
        print(f"Error: File not found: {input_file}")
        sys.exit(1)

    # Pass all remaining arguments to rrna_phylo_cli.py
    args = ' '.join([f'"{arg}"' if ' ' in arg else arg for arg in sys.argv[2:]])
    cmd = f'python rrna_phylo_cli.py {args}'

    print(f"Building trees from: {input_file}")
    print()

    result = subprocess.run(cmd, shell=True)
    sys.exit(result.returncode)


def run_test():
    """Run complete test workflow."""
    import subprocess

    print("=" * 70)
    print("rRNA-PHYLO TEST WORKFLOW")
    print("=" * 70)
    print()

    # Step 1: Create test data
    print("Step 1/3: Creating test data...")
    result = subprocess.run('python rrna_phylo_test.py', shell=True)
    if result.returncode != 0:
        print("Failed to create test data")
        sys.exit(1)

    print()

    # Step 2: Build trees
    print("Step 2/3: Building phylogenetic trees...")
    result = subprocess.run('python rrna_phylo_cli.py test_sequences.fasta', shell=True)
    if result.returncode != 0:
        print("Failed to build trees")
        sys.exit(1)

    print()

    # Step 3: Verify output
    print("Step 3/3: Verifying output structure...")
    result = subprocess.run('python rrna_phylo_test.py --verify', shell=True)

    print()
    print("=" * 70)
    if result.returncode == 0:
        print("TEST PASSED!")
        print("=" * 70)
        print()
        print("Results saved to: results/test_sequences/")
        print()
        print("To clean up: python app.py clean")
    else:
        print("TEST FAILED!")
        print("=" * 70)

    sys.exit(result.returncode)


def run_clean():
    """Clean up test files."""
    import subprocess

    print("=" * 70)
    print("CLEANING UP TEST FILES")
    print("=" * 70)
    print()

    # Confirm with user
    response = input("This will remove test_sequences.fasta and results/. Continue? (yes/no): ")

    if response.lower() != 'yes':
        print("Cleanup cancelled.")
        return

    print()
    result = subprocess.run('bash cleanup_test.sh', shell=True)

    if result.returncode == 0:
        print()
        print("Cleanup complete!")
    else:
        print()
        print("Cleanup script not found. Manual cleanup:")
        print("  rm -rf test_sequences.fasta results/")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Goodbye!")
        sys.exit(0)

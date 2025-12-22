#!/usr/bin/env python3
"""
rRNA-Phylo Interactive Menu

User-friendly interactive menu for building phylogenetic trees without
needing to remember command-line flags.

Usage:
    python menu.py
"""

import os
import sys
from pathlib import Path
import subprocess

# Prevent Intel OpenMP library conflicts
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'


def clear_screen():
    """Clear the terminal screen."""
    os.system('cls' if os.name == 'nt' else 'clear')


def print_header():
    """Print the menu header."""
    print("=" * 70)
    print("rRNA-PHYLO: PHYLOGENETIC TREE BUILDER")
    print("Interactive Menu")
    print("=" * 70)
    print()


def check_ete3():
    """Check if ETE3 is installed."""
    try:
        import ete3
        return True
    except ImportError:
        return False


def main_menu():
    """Display main menu and handle user selection."""
    while True:
        clear_screen()
        print_header()

        print("Main Menu:")
        print()
        print("  1. Prepare FASTA File (Deduplicate + Clean Headers)")
        print("  2. Build Trees (Quick - All 3 Methods)")
        print("  3. Build Trees (Custom - Choose Options)")
        print("  4. Build ML Tree (Advanced)")
        print("  5. Create Test Dataset")
        print("  6. View Results")
        print("  7. Clean Up Test Files")
        print("  8. Install ETE3 (Visualization)")
        print("  9. Help & Documentation")
        print("  10. Exit")
        print()

        choice = input("Select option (1-10): ").strip()

        if choice == '1':
            prepare_fasta_file()
        elif choice == '2':
            quick_build()
        elif choice == '3':
            custom_build()
        elif choice == '4':
            advanced_ml()
        elif choice == '5':
            create_test_data()
        elif choice == '6':
            view_results()
        elif choice == '7':
            cleanup_results()
        elif choice == '8':
            install_ete3()
        elif choice == '9':
            show_help()
        elif choice == '10':
            print("\nGoodbye!")
            sys.exit(0)
        else:
            print("\nInvalid choice. Press Enter to continue...")
            input()


def select_input_file():
    """Let user select file from data folder or provide custom path."""
    data_dir = Path('data')

    # Check for files in data folder
    if data_dir.exists():
        fasta_files = list(data_dir.glob('*.fasta')) + list(data_dir.glob('*.fa')) + list(data_dir.glob('*.fna'))
    else:
        fasta_files = []

    if fasta_files:
        print("Files in data/ folder:")
        print()
        for i, f in enumerate(fasta_files, 1):
            file_size = f.stat().st_size
            size_kb = file_size / 1024
            print(f"  {i}. {f.name} ({size_kb:.1f} KB)")
        print()
        print(f"  {len(fasta_files) + 1}. Enter custom file path")
        print()

        choice = input(f"Select file (1-{len(fasta_files) + 1}): ").strip()

        if choice.isdigit() and 1 <= int(choice) <= len(fasta_files):
            return str(fasta_files[int(choice) - 1])
        elif choice.isdigit() and int(choice) == len(fasta_files) + 1:
            # Fall through to custom path
            pass
        else:
            return None

    # Custom path
    print("Enter FASTA file path:")
    print("  (Relative to current directory or absolute path)")
    print()
    input_file = input("Path: ").strip()

    if not input_file:
        return None

    return input_file


def prepare_fasta_file():
    """Prepare FASTA file: deduplicate + clean headers."""
    clear_screen()
    print_header()
    print("Prepare FASTA File")
    print("=" * 70)
    print()
    print("This tool prepares raw FASTA files for phylogenetic analysis by:")
    print("  1. Deduplicating sequences (keeps longest per species)")
    print("  2. Cleaning headers to standard ID|Species_name format")
    print()

    # Get input file
    input_file = select_input_file()

    if not input_file:
        print("\nNo file specified.")
        input("\nPress Enter to continue...")
        return

    if not Path(input_file).exists():
        print(f"\nError: File not found: {input_file}")
        input("\nPress Enter to continue...")
        return

    # Get output file
    print()
    print("Output file name:")
    print("  (Suggest adding '_clean' or '_prepared' suffix)")
    print()

    # Suggest output name
    input_path = Path(input_file)
    suggested_output = input_path.parent / f"{input_path.stem}_clean{input_path.suffix}"
    print(f"  Default: {suggested_output}")
    print()

    output_file = input("Output file (Enter for default): ").strip()
    if not output_file:
        output_file = str(suggested_output)

    # Confirm overwrite if exists
    if Path(output_file).exists():
        print()
        print(f"Warning: {output_file} already exists!")
        confirm = input("Overwrite? (y/n): ").strip().lower()
        if confirm != 'y':
            print("\nCancelled.")
            input("\nPress Enter to continue...")
            return

    print()
    print("=" * 70)
    print()

    # Run preparation
    result = subprocess.run(
        f'python prepare_fasta.py "{input_file}" "{output_file}"',
        shell=True
    )

    if result.returncode == 0:
        print()
        print("=" * 70)
        print("SUCCESS!")
        print("=" * 70)
        print()
        print(f"Prepared data saved to: {output_file}")
        print()
        print("You can now use this file for tree building (Options 2, 3, or 4)")
    else:
        print()
        print("Preparation failed. Check the error messages above.")

    input("\nPress Enter to continue...")


def quick_build():
    """Quick build with default settings."""
    clear_screen()
    print_header()
    print("Quick Build - All Three Tree Types")
    print("=" * 70)
    print()

    # Get input file
    input_file = select_input_file()

    if not input_file:
        print("\nNo file specified.")
        input("\nPress Enter to continue...")
        return

    if not Path(input_file).exists():
        print(f"\nError: File not found: {input_file}")
        input("\nPress Enter to continue...")
        return

    # Ask about visualization
    print()
    if check_ete3():
        print("ETE3 is installed. Visualization available.")
        viz = input("Create PDF visualization? (y/n): ").strip().lower()
        visualize = "--visualize pdf" if viz == 'y' else ""
    else:
        print("ETE3 not installed. Visualization disabled.")
        print("(Install with option 7 from main menu)")
        visualize = ""

    # Build command
    cmd = f'python rrna_phylo_cli.py "{input_file}" {visualize}'

    print()
    print("=" * 70)
    print("Building trees...")
    print("=" * 70)
    print()

    # Run command
    result = subprocess.run(cmd, shell=True)

    print()
    if result.returncode == 0:
        print("=" * 70)
        print("SUCCESS!")
        print("=" * 70)
    else:
        print("=" * 70)
        print("Build failed. Check errors above.")
        print("=" * 70)

    input("\nPress Enter to continue...")


def custom_build():
    """Custom build with user-selected options."""
    clear_screen()
    print_header()
    print("Custom Build")
    print("=" * 70)
    print()

    # Get input file
    input_file = select_input_file()

    if not input_file:
        print("\nNo file specified.")
        input("\nPress Enter to continue...")
        return

    if not Path(input_file).exists():
        print(f"\nError: File not found: {input_file}")
        input("\nPress Enter to continue...")
        return

    # Choose method
    print()
    print("Select method:")
    print("  1. All three (UPGMA + BioNJ + ML)")
    print("  2. UPGMA only (fastest)")
    print("  3. BioNJ only (fast)")
    print("  4. ML only (most accurate)")

    method_choice = input("\nChoice (1-4): ").strip()
    method_map = {'1': 'all', '2': 'upgma', '3': 'bionj', '4': 'ml'}
    method = method_map.get(method_choice, 'all')

    # Bootstrap
    print()
    bootstrap = input("Bootstrap replicates (0 = disabled, 100 = publication): ").strip()
    bootstrap = bootstrap if bootstrap and bootstrap.isdigit() else "0"

    # Visualization
    print()
    if check_ete3():
        print("Select visualization format:")
        print("  1. None (ASCII only)")
        print("  2. PDF (vector, best for publication)")
        print("  3. PNG (raster, 300 DPI)")
        print("  4. PNG (high-res, 600 DPI)")
        print("  5. SVG (vector, editable)")

        viz_choice = input("\nChoice (1-5): ").strip()
        viz_map = {
            '1': '',
            '2': '--visualize pdf',
            '3': '--visualize png --dpi 300',
            '4': '--visualize png --dpi 600',
            '5': '--visualize svg'
        }
        visualize = viz_map.get(viz_choice, '')
    else:
        print("ETE3 not installed. Only ASCII visualization available.")
        print("(Install ETE3 with option 7 from main menu)")
        visualize = ''

    # Custom output directory
    print()
    custom_output = input("Custom output directory (Enter = 'results'): ").strip()
    output_dir = f'--output "{custom_output}"' if custom_output else ''

    # Build command
    bootstrap_flag = f'--bootstrap {bootstrap}' if int(bootstrap) > 0 else ''
    cmd = f'python rrna_phylo_cli.py "{input_file}" --method {method} {bootstrap_flag} {visualize} {output_dir}'

    print()
    print("=" * 70)
    print("Command:")
    print(cmd)
    print("=" * 70)
    print()
    print("Building trees...")
    print()

    # Run command
    result = subprocess.run(cmd, shell=True)

    print()
    if result.returncode == 0:
        print("=" * 70)
        print("SUCCESS!")
        print("=" * 70)
    else:
        print("=" * 70)
        print("Build failed. Check errors above.")
        print("=" * 70)

    input("\nPress Enter to continue...")


def advanced_ml():
    """Advanced ML tree building with model selection."""
    clear_screen()
    print_header()
    print("Advanced Maximum Likelihood Tree")
    print("=" * 70)
    print()

    # Get input file
    input_file = select_input_file()

    if not input_file:
        print("\nNo file specified.")
        input("\nPress Enter to continue...")
        return

    if not Path(input_file).exists():
        print(f"\nError: File not found: {input_file}")
        input("\nPress Enter to continue...")
        return

    # Tree search method
    print()
    print("Select tree search method:")
    print("  1. NNI (fast, local search)")
    print("  2. SPR (thorough, global search)")

    method_choice = input("\nChoice (1-2): ").strip()
    method = 'spr' if method_choice == '2' else 'nni'

    # Bootstrap
    print()
    bootstrap = input("Bootstrap replicates (0 = disabled, 100 = publication): ").strip()
    bootstrap = bootstrap if bootstrap and bootstrap.isdigit() else "0"

    # Visualization
    print()
    if check_ete3():
        print("Select visualization format:")
        print("  1. None")
        print("  2. PDF")
        print("  3. PNG (600 DPI)")
        print("  4. SVG")

        viz_choice = input("\nChoice (1-4): ").strip()
        viz_map = {
            '1': '',
            '2': '--visualize pdf',
            '3': '--visualize png --dpi 600',
            '4': '--visualize svg'
        }
        visualize = viz_map.get(viz_choice, '')
    else:
        print("ETE3 not installed. Only ASCII visualization available.")
        visualize = ''

    # Build command
    bootstrap_flag = f'--bootstrap {bootstrap}' if int(bootstrap) > 0 else ''
    cmd = f'python rrna_phylo_ml.py "{input_file}" --method {method} {bootstrap_flag} {visualize}'

    print()
    print("=" * 70)
    print("Building ML tree with automatic model selection...")
    print("=" * 70)
    print()

    # Run command
    result = subprocess.run(cmd, shell=True)

    print()
    if result.returncode == 0:
        print("=" * 70)
        print("SUCCESS!")
        print("=" * 70)
    else:
        print("=" * 70)
        print("Build failed. Check errors above.")
        print("=" * 70)

    input("\nPress Enter to continue...")


def create_test_data():
    """Create test dataset."""
    clear_screen()
    print_header()
    print("Create Test Dataset")
    print("=" * 70)
    print()

    print("This will create test_sequences.fasta with:")
    print("  - 50 sequences")
    print("  - 1000 bp each")
    print("  - Realistic mutations")
    print()

    confirm = input("Create test dataset? (y/n): ").strip().lower()

    if confirm == 'y':
        print()
        result = subprocess.run('python rrna_phylo_test.py', shell=True)

        if result.returncode == 0:
            print()
            print("Test dataset created: test_sequences.fasta")
            print()
            print("You can now:")
            print("  1. Use option 1 or 2 to build trees")
            print("  2. Specify 'test_sequences.fasta' as input file")

    input("\nPress Enter to continue...")


def view_results():
    """View available results."""
    clear_screen()
    print_header()
    print("View Results")
    print("=" * 70)
    print()

    results_dir = Path('results')

    if not results_dir.exists():
        print("No results found.")
        print()
        print("Results are stored in: results/[filename]/")
        input("\nPress Enter to continue...")
        return

    # List all result directories
    subdirs = [d for d in results_dir.iterdir() if d.is_dir()]

    if not subdirs:
        print("No results found.")
        input("\nPress Enter to continue...")
        return

    print("Available results:")
    print()

    for i, subdir in enumerate(subdirs, 1):
        files = list(subdir.glob('*'))
        file_count = len(files)

        # Get file types
        has_newick = any(f.suffix == '.nwk' for f in files)
        has_pdf = any(f.suffix == '.pdf' for f in files)
        has_png = any(f.suffix == '.png' for f in files)
        has_svg = any(f.suffix == '.svg' for f in files)

        types = []
        if has_newick:
            types.append('Newick')
        if has_pdf:
            types.append('PDF')
        if has_png:
            types.append('PNG')
        if has_svg:
            types.append('SVG')

        type_str = ', '.join(types) if types else 'ASCII only'

        print(f"  {i}. {subdir.name}")
        print(f"     Files: {file_count} ({type_str})")
        print(f"     Path: {subdir.absolute()}")
        print()

    print("To view files, use your file manager or:")
    print(f"  cd {results_dir.absolute()}")

    input("\nPress Enter to continue...")


def cleanup_results():
    """Clean up test files and results."""
    clear_screen()
    print_header()
    print("Clean Up Test Files")
    print("=" * 70)
    print()

    print("This will remove:")
    print("  - test_sequences.fasta (generated test data)")
    print("  - results/test_sequences/ (test results only)")
    print()
    print("This will NOT remove:")
    print("  - data/ folder and your prepared FASTA files")
    print("  - results/ for other datasets")
    print()
    print("WARNING: This cannot be undone!")
    print()

    confirm = input("Continue with cleanup? (yes/no): ").strip().lower()

    if confirm == 'yes':
        print()
        result = subprocess.run('bash cleanup_test.sh', shell=True)

        if result.returncode == 0:
            print()
            print("Cleanup complete!")
        else:
            print()
            print("Cleanup script not found. Manual cleanup:")
            print("  rm -rf test_sequences.fasta results/")

    input("\nPress Enter to continue...")


def install_ete3():
    """Install ETE3 for visualization."""
    clear_screen()
    print_header()
    print("Install ETE3")
    print("=" * 70)
    print()

    if check_ete3():
        print("ETE3 is already installed!")
        print()
        print("You can now use visualization options:")
        print("  --visualize pdf")
        print("  --visualize png --dpi 600")
        print("  --visualize svg")
        input("\nPress Enter to continue...")
        return

    print("ETE3 enables publication-quality tree visualization:")
    print("  - PDF (vector graphics)")
    print("  - PNG (high-resolution raster)")
    print("  - SVG (editable vector)")
    print()
    print("Installation command:")
    print("  pip install ete3")
    print()

    confirm = input("Install now? (y/n): ").strip().lower()

    if confirm == 'y':
        print()
        print("Installing ETE3...")
        print()
        result = subprocess.run('pip install ete3', shell=True)

        print()
        if result.returncode == 0:
            print("=" * 70)
            print("ETE3 installed successfully!")
            print("=" * 70)
            print()
            print("You can now create PDF/PNG/SVG visualizations.")
        else:
            print("=" * 70)
            print("Installation failed.")
            print("=" * 70)
            print()
            print("Try manually:")
            print("  pip install ete3")

    input("\nPress Enter to continue...")


def show_help():
    """Show help and documentation."""
    clear_screen()
    print_header()
    print("Help & Documentation")
    print("=" * 70)
    print()

    print("QUICK START:")
    print()
    print("1. Create test data (Option 4)")
    print("2. Build trees (Option 1 or 2)")
    print("3. View results (Option 5)")
    print()
    print("=" * 70)
    print()
    print("TREE BUILDING METHODS:")
    print()
    print("• UPGMA")
    print("  - Fastest method")
    print("  - Assumes molecular clock (constant evolution rate)")
    print("  - Good for closely related sequences")
    print()
    print("• BioNJ")
    print("  - Fast, variance-weighted")
    print("  - No clock assumption")
    print("  - Better for real data")
    print()
    print("• Maximum Likelihood (ML)")
    print("  - Most accurate, but slower")
    print("  - Automatic model selection (JC69, K80, F81, HKY85, GTR)")
    print("  - Statistical inference")
    print()
    print("=" * 70)
    print()
    print("VISUALIZATION FORMATS:")
    print()
    print("• PDF - Best for publications (vector graphics)")
    print("• PNG - High-resolution raster (300-1200 DPI)")
    print("• SVG - Editable vector (Illustrator/Inkscape)")
    print()
    print("Requires: pip install ete3 (Option 7)")
    print()
    print("=" * 70)
    print()
    print("BOOTSTRAP ANALYSIS:")
    print()
    print("• Tests tree reliability")
    print("• 100 replicates recommended for publication")
    print("• 10 replicates for quick testing")
    print("• Values >70% indicate strong support")
    print()
    print("=" * 70)
    print()
    print("OUTPUT STRUCTURE:")
    print()
    print("results/")
    print("└── [filename]/")
    print("    ├── upgma_tree.nwk    (Newick format)")
    print("    ├── bionj_tree.nwk    (Newick format)")
    print("    ├── ml_tree.nwk       (Newick format)")
    print("    ├── upgma_tree.txt    (ASCII visualization)")
    print("    ├── bionj_tree.txt    (ASCII visualization)")
    print("    ├── ml_tree.txt       (ASCII visualization)")
    print("    ├── summary.txt       (Comparison)")
    print("    └── *.pdf/png/svg     (If visualization enabled)")
    print()
    print("=" * 70)
    print()
    print("For more information, see README.md")

    input("\nPress Enter to continue...")


if __name__ == '__main__':
    try:
        main_menu()
    except KeyboardInterrupt:
        print("\n\nInterrupted. Goodbye!")
        sys.exit(0)

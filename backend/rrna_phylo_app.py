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
        print("  7. Clean Up Files (Data + Results)")
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
    test_dir = data_dir / 'test'

    # Check for test datasets
    test_files = []
    if test_dir.exists():
        test_files = list(test_dir.glob('*.fasta')) + list(test_dir.glob('*.fa')) + list(test_dir.glob('*.fna'))

    # Check for user data files (root data folder only)
    user_files = []
    if data_dir.exists():
        user_files = [f for f in data_dir.glob('*.fasta') if f.parent == data_dir]
        user_files += [f for f in data_dir.glob('*.fa') if f.parent == data_dir]
        user_files += [f for f in data_dir.glob('*.fna') if f.parent == data_dir]

    all_files = test_files + user_files

    if all_files:
        # Show test datasets first
        if test_files:
            print("Test datasets (data/test/):")
            print()
            for i, f in enumerate(test_files, 1):
                file_size = f.stat().st_size
                size_kb = file_size / 1024
                print(f"  {i}. {f.name} ({size_kb:.1f} KB)")
            print()

        # Then show user data
        if user_files:
            print("Your data (data/):")
            print()
            offset = len(test_files)
            for i, f in enumerate(user_files, offset + 1):
                file_size = f.stat().st_size
                size_kb = file_size / 1024
                print(f"  {i}. {f.name} ({size_kb:.1f} KB)")
            print()

        print(f"  {len(all_files) + 1}. Enter custom file path")
        print()

        choice = input(f"Select file (1-{len(all_files) + 1}): ").strip()

        if choice.isdigit() and 1 <= int(choice) <= len(all_files):
            return str(all_files[int(choice) - 1])
        elif choice.isdigit() and int(choice) == len(all_files) + 1:
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
    """Clean up files and results with user selection."""
    clear_screen()
    print_header()
    print("Clean Up Files")
    print("=" * 70)
    print()

    print("What would you like to clean?")
    print()
    print("  1. Test files only (test_sequences.fasta + results/test_sequences/)")
    print("  2. Select specific result folders to delete")
    print("  3. Delete ALL results (keeps data/ folder)")
    print("  4. Select specific data/ FASTA files to delete")
    print("  5. Delete ALL data/ FASTA files")
    print("  6. Clean everything (data/ + results/)")
    print("  7. Cancel")
    print()

    choice = input("Select option (1-7): ").strip()

    if choice == '1':
        # Clean test files only
        cleanup_test_files()
    elif choice == '2':
        # Select specific results
        cleanup_select_results()
    elif choice == '3':
        # Clean all results
        cleanup_all_results()
    elif choice == '4':
        # Select specific data files
        cleanup_select_data()
    elif choice == '5':
        # Clean all data
        cleanup_all_data()
    elif choice == '6':
        # Clean everything
        cleanup_everything()
    elif choice == '7':
        pass  # Cancel
    else:
        print("\nInvalid choice.")

    input("\nPress Enter to continue...")


def cleanup_test_files():
    """Clean up test files only."""
    print()
    print("=" * 70)
    print("CLEAN TEST FILES")
    print("=" * 70)
    print()
    print("This will remove:")
    print("  - test_sequences.fasta")
    print("  - results/test_sequences/")
    print()

    confirm = input("Continue? (yes/no): ").strip().lower()
    if confirm == 'yes':
        result = subprocess.run('bash cleanup_test.sh', shell=True)
        if result.returncode == 0:
            print("\nTest files cleaned!")


def cleanup_select_results():
    """Let user select which result folders to delete."""
    print()
    print("=" * 70)
    print("SELECT RESULTS TO DELETE")
    print("=" * 70)
    print()

    results_dir = Path('results')
    if not results_dir.exists() or not any(results_dir.iterdir()):
        print("No results found.")
        return

    folders = [f for f in results_dir.iterdir() if f.is_dir()]
    if not folders:
        print("No result folders found.")
        return

    print("Result folders:")
    print()
    for i, folder in enumerate(folders, 1):
        size = sum(f.stat().st_size for f in folder.rglob('*') if f.is_file())
        size_mb = size / (1024 * 1024)
        print(f"  {i}. {folder.name} ({size_mb:.1f} MB)")
    print()

    selection = input("Enter numbers to delete (comma-separated, or 'all'): ").strip()

    if selection.lower() == 'all':
        to_delete = folders
    else:
        try:
            indices = [int(x.strip()) - 1 for x in selection.split(',')]
            to_delete = [folders[i] for i in indices if 0 <= i < len(folders)]
        except (ValueError, IndexError):
            print("\nInvalid selection.")
            return

    if not to_delete:
        return

    print()
    print("Will delete:")
    for folder in to_delete:
        print(f"  - {folder}")
    print()

    confirm = input("Continue? (yes/no): ").strip().lower()
    if confirm == 'yes':
        import shutil
        for folder in to_delete:
            shutil.rmtree(folder)
            print(f"Deleted: {folder}")
        print("\nSelected results deleted!")


def cleanup_all_results():
    """Delete all results."""
    print()
    print("=" * 70)
    print("DELETE ALL RESULTS")
    print("=" * 70)
    print()

    results_dir = Path('results')
    if not results_dir.exists():
        print("No results directory found.")
        return

    print("WARNING: This will delete the entire results/ directory!")
    print()

    confirm = input("Continue? (yes/no): ").strip().lower()
    if confirm == 'yes':
        import shutil
        shutil.rmtree(results_dir)
        print("\nAll results deleted!")


def cleanup_select_data():
    """Let user select which data files to delete."""
    print()
    print("=" * 70)
    print("SELECT DATA FILES TO DELETE")
    print("=" * 70)
    print()

    data_dir = Path('data')
    if not data_dir.exists():
        print("No data directory found.")
        return

    # Only list files from root data folder (not test/)
    fasta_files = [f for f in data_dir.glob('*.fasta') if f.parent == data_dir]
    fasta_files += [f for f in data_dir.glob('*.fa') if f.parent == data_dir]
    fasta_files += [f for f in data_dir.glob('*.fna') if f.parent == data_dir]

    if not fasta_files:
        print("No FASTA files found in data/.")
        print("(Test datasets in data/test/ are protected)")
        return

    print("FASTA files in data/:")
    print("(Test datasets in data/test/ are protected)")
    print()
    for i, f in enumerate(fasta_files, 1):
        size_kb = f.stat().st_size / 1024
        print(f"  {i}. {f.name} ({size_kb:.1f} KB)")
    print()

    selection = input("Enter numbers to delete (comma-separated, or 'all'): ").strip()

    if selection.lower() == 'all':
        to_delete = fasta_files
    else:
        try:
            indices = [int(x.strip()) - 1 for x in selection.split(',')]
            to_delete = [fasta_files[i] for i in indices if 0 <= i < len(fasta_files)]
        except (ValueError, IndexError):
            print("\nInvalid selection.")
            return

    if not to_delete:
        return

    print()
    print("Will delete:")
    for f in to_delete:
        print(f"  - {f}")
    print()

    confirm = input("Continue? (yes/no): ").strip().lower()
    if confirm == 'yes':
        for f in to_delete:
            f.unlink()
            print(f"Deleted: {f}")
        print("\nSelected files deleted!")


def cleanup_all_data():
    """Delete all FASTA files from data/."""
    print()
    print("=" * 70)
    print("DELETE ALL DATA FILES")
    print("=" * 70)
    print()

    data_dir = Path('data')
    if not data_dir.exists():
        print("No data directory found.")
        return

    # Only delete files from root data folder (not test/)
    fasta_files = [f for f in data_dir.glob('*.fasta') if f.parent == data_dir]
    fasta_files += [f for f in data_dir.glob('*.fa') if f.parent == data_dir]
    fasta_files += [f for f in data_dir.glob('*.fna') if f.parent == data_dir]

    if not fasta_files:
        print("No FASTA files found in data/.")
        print("(Test datasets in data/test/ are protected)")
        return

    print(f"WARNING: This will delete ALL {len(fasta_files)} FASTA files from data/!")
    print("(Test datasets in data/test/ will be preserved)")
    print()

    confirm = input("Continue? (yes/no): ").strip().lower()
    if confirm == 'yes':
        for f in fasta_files:
            f.unlink()
        print(f"\nDeleted {len(fasta_files)} FASTA files!")


def cleanup_everything():
    """Delete everything (data + results)."""
    print()
    print("=" * 70)
    print("DELETE EVERYTHING")
    print("=" * 70)
    print()

    print("WARNING: This will delete:")
    print("  - All FASTA files in data/ (except test/)")
    print("  - All result folders in results/")
    print()
    print("This CANNOT be undone!")
    print()

    confirm = input("Are you SURE? Type 'DELETE EVERYTHING' to confirm: ").strip()
    if confirm == 'DELETE EVERYTHING':
        import shutil

        # Delete data FASTA files (only from root data folder, not test/)
        data_dir = Path('data')
        if data_dir.exists():
            fasta_files = [f for f in data_dir.glob('*.fasta') if f.parent == data_dir]
            fasta_files += [f for f in data_dir.glob('*.fa') if f.parent == data_dir]
            fasta_files += [f for f in data_dir.glob('*.fna') if f.parent == data_dir]
            for f in fasta_files:
                f.unlink()
            print(f"Deleted {len(fasta_files)} data files (test datasets preserved)")

        # Delete results
        results_dir = Path('results')
        if results_dir.exists():
            shutil.rmtree(results_dir)
            print("Deleted results/")

        print("\nEverything deleted!")


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

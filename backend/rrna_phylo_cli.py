#!/usr/bin/env python3
"""
rRNA-Phylo: Build All Tree Types with Visualization

Simple script to build phylogenetic trees using all three methods:
- UPGMA (fast, assumes molecular clock)
- BioNJ (fast, no clock assumption)
- Maximum Likelihood (rigorous, slower)

Usage:
    python build_trees.py sequences.fasta
    python build_trees.py sequences.fasta --bootstrap 100
    python build_trees.py sequences.fasta --method upgma
    python build_trees.py sequences.fasta --output my_results

Output:
    results/[filename]/
        upgma_tree.nwk          - UPGMA tree (Newick)
        bionj_tree.nwk          - BioNJ tree (Newick)
        ml_tree.nwk             - ML tree (Newick)
        upgma_tree.txt          - UPGMA tree (ASCII visualization)
        bionj_tree.txt          - BioNJ tree (ASCII visualization)
        ml_tree.txt             - ML tree (ASCII visualization)
        summary.txt             - Comparison of all methods

    [filename] is the input file name without extension
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

# Prevent Intel OpenMP library conflicts
# This is safe and necessary when multiple packages (NumPy, SciPy, PyTorch) bundle MKL
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

sys.path.insert(0, str(Path(__file__).parent))

from rrna_phylo.io.fasta_parser import FastaParser
from rrna_phylo.core.builder import PhylogeneticTreeBuilder

# Check if ETE3 is available
try:
    from rrna_phylo.visualization.ete3_viz import visualize_tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False


def visualize_tree_ascii(tree, title="PHYLOGENETIC TREE", max_width=120):
    """Create ASCII visualization of tree with proper branch diagram."""
    lines = []
    lines.append("=" * 70)
    lines.append(title)
    lines.append("=" * 70)
    lines.append("")

    # Draw the actual tree structure
    def _draw_tree(node, prefix="", is_last=True, is_root=True):
        """Recursively draw tree structure."""
        result = []

        if node.is_leaf():
            # Leaf node - show species name with branch length
            branch_char = "`-" if is_last else "|-"
            if is_root:
                result.append(f"{node.name} ({node.distance:.6f})")
            else:
                result.append(f"{prefix}{branch_char} {node.name} ({node.distance:.6f})")
        else:
            # Internal node
            if not is_root:
                branch_char = "`-" if is_last else "|-"
                support_str = f" [{node.support:.1f}]" if node.support is not None else ""
                result.append(f"{prefix}{branch_char} Internal{support_str} ({node.distance:.6f})")
                prefix += "   " if is_last else "|  "

            # Draw children (right subtree first for better visual ordering)
            if node.right:
                result.extend(_draw_tree(node.right, prefix, False, False))
            if node.left:
                result.extend(_draw_tree(node.left, prefix, True, False))

        return result

    tree_lines = _draw_tree(tree)
    lines.extend(tree_lines)
    lines.append("")

    # Add summary statistics
    taxa = tree.get_leaf_names()
    lines.append(f"Taxa: {len(taxa)}")
    lines.append("")

    # Add Newick format for reference
    lines.append("Newick Format:")
    lines.append(tree.to_newick() + ";")
    lines.append("")

    return "\n".join(lines)


def build_all_trees(input_file, output_dir="results", method="all", bootstrap=0,
                   visualize=None, dpi=300):
    """
    Build phylogenetic trees using selected method(s).

    Args:
        input_file: Path to FASTA file (aligned or unaligned)
        output_dir: Output directory
        method: Which method(s) to use ('all', 'upgma', 'bionj', 'ml')
        bootstrap: Number of bootstrap replicates (0 = disabled)
        visualize: Visualization format ('pdf', 'png', 'svg', or None)
        dpi: Resolution for PNG output (default: 300)
    """
    # Create output directory with subfolder named after input file
    input_path = Path(input_file)
    input_basename = input_path.stem  # Get filename without extension

    output_path = Path(output_dir) / input_basename
    output_path.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("rRNA-PHYLO: PHYLOGENETIC TREE BUILDER")
    print("=" * 70)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input: {input_file}")
    print(f"Method: {method.upper()}")
    print(f"Bootstrap: {bootstrap if bootstrap > 0 else 'Disabled'}")
    print()

    # Load sequences
    print("[1/3] Loading sequences...")
    parser = FastaParser()
    sequences = parser.parse(input_file)
    print(f"  Loaded: {len(sequences)} sequences")
    print()

    # Build trees
    print("[2/3] Building phylogenetic trees...")
    builder = PhylogeneticTreeBuilder(verbose=True)

    # Detect sequence type and align if needed
    builder.detect_and_validate(sequences)
    aligned_seqs = builder._align_sequences(sequences)

    trees = {}

    if method in ['all', 'upgma']:
        print("\n--- UPGMA Tree ---")
        upgma_tree = builder.build_upgma_tree(aligned_seqs)
        trees['upgma'] = ('UPGMA', upgma_tree, None)
        print()

    if method in ['all', 'bionj']:
        print("\n--- BioNJ Tree ---")
        bionj_tree = builder.build_bionj_tree(aligned_seqs)
        trees['bionj'] = ('BioNJ', bionj_tree, None)
        print()

    if method in ['all', 'ml']:
        print("\n--- Maximum Likelihood Tree ---")
        ml_tree, logL = builder.build_ml_tree(aligned_seqs, alpha=1.0)
        trees['ml'] = ('Maximum Likelihood', ml_tree, logL)
        print()

    # Save results
    print("[3/3] Saving results...")

    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("TREE BUILDING SUMMARY")
    summary_lines.append("=" * 70)
    summary_lines.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append(f"Input: {input_file}")
    summary_lines.append(f"Sequences: {len(sequences)}")
    summary_lines.append(f"Alignment length: {aligned_seqs[0].aligned_length} bp")
    summary_lines.append("")

    for key, (name, tree, logL) in trees.items():
        # Save Newick format
        newick_file = output_path / f"{key}_tree.nwk"
        with open(newick_file, 'w') as f:
            f.write(tree.to_newick() + ";")
        print(f"  Saved: {newick_file}")

        # Save ASCII visualization
        ascii_file = output_path / f"{key}_tree.txt"
        ascii_content = visualize_tree_ascii(tree, title=f"{name.upper()} TREE")
        with open(ascii_file, 'w') as f:
            f.write(ascii_content)
        print(f"  Saved: {ascii_file}")

        # Save ETE3 visualization if requested
        if visualize:
            if not ETE3_AVAILABLE:
                print(f"  Warning: ETE3 not installed, skipping {visualize} visualization")
                print(f"  Install with: pip install ete3")
            else:
                viz_file = output_path / f"{key}_tree.{visualize}"
                try:
                    success = visualize_tree(
                        str(newick_file),
                        str(viz_file),
                        show_bootstrap=bootstrap > 0,
                        show_branch_length=True,
                        dpi=dpi,
                        title=f"{name} Tree"
                    )
                    if success:
                        print(f"  Saved: {viz_file}")
                    else:
                        print(f"  Warning: Failed to create {viz_file}")
                except Exception as e:
                    print(f"  Warning: Visualization failed: {e}")

        # Add to summary
        summary_lines.append(f"{name}:")
        summary_lines.append(f"  File: {key}_tree.nwk")
        if logL is not None:
            summary_lines.append(f"  Log-likelihood: {logL:.4f}")
        summary_lines.append(f"  Taxa: {len(tree.get_leaf_names())}")
        summary_lines.append("")

    # Save summary
    summary_file = output_path / "summary.txt"
    summary_lines.append("=" * 70)
    summary_lines.append(f"All results saved to: {output_path.absolute()}")
    summary_lines.append("=" * 70)

    with open(summary_file, 'w') as f:
        f.write("\n".join(summary_lines))
    print(f"  Saved: {summary_file}")

    print()
    print("=" * 70)
    print("BUILD COMPLETE!")
    print("=" * 70)
    print(f"Results saved to: {output_path.absolute()}")
    print()

    # Print quick comparison
    print("Quick Comparison:")
    for key, (name, tree, logL) in trees.items():
        if logL is not None:
            print(f"  {name:<25} LogL: {logL:>10.4f}")
        else:
            print(f"  {name:<25} (distance-based)")
    print()


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description='Build phylogenetic trees using UPGMA, BioNJ, and/or ML',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build all three tree types (UPGMA, BioNJ, ML)
  python build_trees.py sequences.fasta

  # Build only UPGMA tree (fastest)
  python build_trees.py sequences.fasta --method upgma

  # Build only ML tree (most rigorous)
  python build_trees.py sequences.fasta --method ml

  # Build with bootstrap support (slow)
  python build_trees.py sequences.fasta --bootstrap 100

  # Custom output directory
  python build_trees.py sequences.fasta --output my_results

  # Create PDF visualization (requires: pip install ete3)
  python build_trees.py sequences.fasta --visualize pdf

  # Create high-resolution PNG for publication (600 DPI)
  python build_trees.py sequences.fasta --visualize png --dpi 600

  # Create SVG vector graphics
  python build_trees.py sequences.fasta --visualize svg

Output Files:
  results/[filename]/upgma_tree.nwk    - UPGMA tree (Newick format)
  results/[filename]/bionj_tree.nwk    - BioNJ tree (Newick format)
  results/[filename]/ml_tree.nwk       - ML tree (Newick format)
  results/[filename]/upgma_tree.txt    - UPGMA tree (ASCII visualization)
  results/[filename]/bionj_tree.txt    - BioNJ tree (ASCII visualization)
  results/[filename]/ml_tree.txt       - ML tree (ASCII visualization)
  results/[filename]/summary.txt       - Comparison of all methods

  Note: [filename] is the input file name without extension

Methods:
  UPGMA  - Fast, assumes molecular clock (all lineages evolve at same rate)
  BioNJ  - Fast, no clock assumption (better for real data)
  ML     - Rigorous maximum likelihood (GTR+Gamma model, automatic selection)
        """
    )

    parser.add_argument('input', help='FASTA file (aligned or unaligned)')
    parser.add_argument('--method', choices=['all', 'upgma', 'bionj', 'ml'],
                       default='all',
                       help='Which method(s) to use (default: all)')
    parser.add_argument('--bootstrap', type=int, default=0,
                       help='Bootstrap replicates (0 = disabled, default: 0)')
    parser.add_argument('--output', default='results',
                       help='Output directory (default: results)')
    parser.add_argument('--visualize', choices=['pdf', 'png', 'svg'],
                       help='Create publication-quality visualization (requires ETE3)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for PNG output (default: 300, use 600 for publication)')

    args = parser.parse_args()

    # Validate input
    if not Path(args.input).exists():
        print(f'Error: Input file not found: {args.input}')
        sys.exit(1)

    # Build trees
    try:
        build_all_trees(
            args.input,
            output_dir=args.output,
            method=args.method,
            bootstrap=args.bootstrap,
            visualize=args.visualize,
            dpi=args.dpi
        )
    except Exception as e:
        print(f'\nError: {e}')
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

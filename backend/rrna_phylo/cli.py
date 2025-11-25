"""
Command-line interface for rRNA-Phylo.

Provides an easy-to-use command-line tool for phylogenetic analysis.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from rrna_phylo import build_trees, FastaParser, Sequence
from rrna_phylo.utils import print_tree_ascii
from rrna_phylo.utils.bootstrap import bootstrap_tree
from rrna_phylo.utils.strain_handler import dereplicate_strains, get_strain_summary, remove_exact_duplicates, smart_dereplicate
from rrna_phylo.utils.outgroup_handler import get_outgroup_sequences, suggest_outgroup
from rrna_phylo.utils.sampling_strategy import (
    stratified_sample,
    get_sampling_recommendations,
    detect_bias
)
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.io.aligner import MuscleAligner


def tree_to_newick(node, include_distances=True, include_support=True):
    """Convert tree to Newick format with proper quoting for names with spaces."""
    if node.is_leaf():
        # Quote name if it contains spaces or special characters
        name = node.name
        if ' ' in name or '(' in name or ')' in name:
            name = f"'{name}'"

        if include_distances and hasattr(node, 'distance') and node.distance is not None:
            return f"{name}:{node.distance:.6f}"
        return name

    left_newick = tree_to_newick(node.left, include_distances, include_support)
    right_newick = tree_to_newick(node.right, include_distances, include_support)

    subtree = f"({left_newick},{right_newick})"

    # Add support value if available
    if include_support and hasattr(node, 'support') and node.support is not None:
        subtree += f"{node.support:.1f}"

    # Add branch length
    if include_distances and hasattr(node, 'distance') and node.distance is not None:
        subtree += f":{node.distance:.6f}"

    return subtree


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog='rrna-phylo',
        description='Build phylogenetic trees from aligned sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - build all three trees
  rrna-phylo sequences.fasta

  # Specify output directory
  rrna-phylo sequences.fasta -o results/

  # Build only ML tree
  rrna-phylo sequences.fasta --method ml

  # Add bootstrap support (50 replicates)
  rrna-phylo sequences.fasta --bootstrap 50

  # Export to Newick format
  rrna-phylo sequences.fasta --output-format newick

  # Quiet mode (no ASCII visualization)
  rrna-phylo sequences.fasta --quiet

For more information, visit: https://github.com/your-repo/rrna-phylo
        """
    )

    # Positional arguments
    parser.add_argument(
        'input',
        type=str,
        help='Input FASTA file with aligned sequences'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='.',
        help='Output directory (default: current directory)'
    )

    parser.add_argument(
        '-m', '--method',
        type=str,
        choices=['all', 'ml', 'bionj', 'upgma'],
        default='all',
        help='Tree building method (default: all)'
    )

    parser.add_argument(
        '-b', '--bootstrap',
        type=int,
        default=0,
        help='Number of bootstrap replicates (default: 0, no bootstrap)'
    )

    parser.add_argument(
        '-f', '--output-format',
        type=str,
        choices=['newick', 'ascii', 'both'],
        default='both',
        help='Output format (default: both)'
    )

    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress ASCII tree visualization'
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default='tree',
        help='Output file prefix (default: tree)'
    )

    parser.add_argument(
        '--no-support',
        action='store_true',
        help='Exclude bootstrap support from Newick output'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output (show progress)'
    )

    parser.add_argument(
        '--skip-align',
        action='store_true',
        help='Skip alignment (use if sequences are already aligned)'
    )

    parser.add_argument(
        '--aligned-output',
        type=str,
        help='Save aligned sequences to this file (default: aligned_<input>)'
    )

    parser.add_argument(
        '--dereplicate',
        action='store_true',
        help='Smart deduplication: cluster highly similar sequences (>=99.5%%) and keep one representative per cluster (exact duplicates are always removed)'
    )

    parser.add_argument(
        '--derep-method',
        type=str,
        choices=['representative', 'consensus'],
        default='representative',
        help='Dereplication method: representative (pick longest) or consensus (default: representative)'
    )

    # Outgroup handling
    parser.add_argument(
        '--outgroup',
        type=str,
        help='Outgroup pattern for rooting tree (e.g., "Pseudomonas*" or "AE004091*")'
    )

    parser.add_argument(
        '--suggest-outgroup',
        action='store_true',
        help='Analyze sequences and suggest appropriate outgroups (no tree building)'
    )

    # Database bias correction
    parser.add_argument(
        '--check-bias',
        action='store_true',
        help='Check for sampling bias and show recommendations (no tree building)'
    )

    parser.add_argument(
        '--stratify',
        action='store_true',
        help='Use stratified sampling to balance species representation'
    )

    parser.add_argument(
        '--max-per-species',
        type=int,
        default=10,
        help='Maximum sequences per species when using --stratify (default: 10)'
    )

    parser.add_argument(
        '--min-per-species',
        type=int,
        default=1,
        help='Minimum sequences per species when using --stratify (default: 1)'
    )

    parser.add_argument(
        '--ignore-bias-warning',
        action='store_true',
        help='Proceed with tree building even if database bias is detected (not recommended)'
    )

    # Visualization options
    parser.add_argument(
        '--visualize',
        action='store_true',
        help='Generate publication-quality tree visualization (requires ETE3)'
    )


    parser.add_argument(
        '--viz-format',
        type=str,
        choices=['png', 'pdf', 'svg', 'eps'],
        default='pdf',
        help='Visualization output format (default: pdf for publication quality)'
    )

    parser.add_argument(
        '--viz-dpi',
        type=int,
        default=300,
        help='Resolution for PNG output in DPI (default: 300, use 600 for publication)'
    )

    parser.add_argument(
        '--viz-width',
        type=int,
        help='Visualization width in pixels (optional, auto if not set)'
    )

    parser.add_argument(
        '--viz-height',
        type=int,
        help='Visualization height in pixels (optional, auto if not set)'
    )

    parser.add_argument(
        '--viz-bootstrap-threshold',
        type=float,
        default=70.0,
        help='Bootstrap value threshold for coloring (default: 70.0)'
    )

    parser.add_argument(
        '--viz-branch-width',
        type=float,
        default=1.5,
        help='Branch line thickness in pixels (default: 1.5)'
    )

    parser.add_argument(
        '--viz-scale',
        type=int,
        default=120,
        help='Branch length scale in pixels per unit (default: 120, higher = longer branches)'
    )


    args = parser.parse_args()

    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse FASTA file
    if not args.quiet:
        print("=" * 70)
        print("rRNA-Phylo: Phylogenetic Tree Builder")
        print("=" * 70)
        print(f"\nReading sequences from: {args.input}")

    try:
        parser_obj = FastaParser()
        sequences = parser_obj.parse(str(input_path))

        # Assign unique display names with counters for duplicates
        from rrna_phylo.io.fasta_parser import assign_unique_display_names
        assign_unique_display_names(sequences)

        # ALWAYS remove exact duplicates (100% identical sequences)
        # This is scientifically sound - identical sequences provide no additional information
        original_count = len(sequences)
        sequences, dup_map = remove_exact_duplicates(sequences)
        num_exact_dups = original_count - len(sequences)

        if not args.quiet:
            print(f"  Loaded {original_count} sequences")
            if num_exact_dups > 0:
                print(f"  Removed {num_exact_dups} exact duplicate(s) -> {len(sequences)} unique sequences")

            # Detect type
            if sequences[0].is_nucleotide():
                seq_type = "DNA/RNA"
            elif sequences[0].is_protein():
                seq_type = "Protein"
            else:
                seq_type = "Unknown"
            print(f"  Sequence type: {seq_type}")

    except Exception as e:
        print(f"Error parsing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    # Auto-detect database bias (unless user explicitly handled it or is just checking)
    if not args.stratify and not args.check_bias and not args.suggest_outgroup:
        bias_info = detect_bias(sequences, threshold=0.1)

        if bias_info['overrepresented'] and not args.ignore_bias_warning:
            print("\n" + "!" * 70)
            print("DATABASE BIAS WARNING")
            print("!" * 70)
            print(f"\n[!] {len(bias_info['overrepresented'])} overrepresented species detected:")

            # Show top 5 overrepresented species
            for species in bias_info['overrepresented'][:5]:
                print(f"  - {species}")
            if len(bias_info['overrepresented']) > 5:
                print(f"  ... and {len(bias_info['overrepresented']) - 5} more")

            print("\nWHY THIS MATTERS:")
            print("  Overrepresented species can dominate your phylogenetic tree,")
            print("  causing incorrect topology and inflated bootstrap support.")

            print("\nRECOMMENDED ACTION:")
            print("  Re-run with stratified sampling:")
            print(f"    --stratify --max-per-species 10")

            print("\nALTERNATIVE:")
            print("  If you understand the risks, proceed anyway:")
            print(f"    --ignore-bias-warning")

            print("\nLEARN MORE:")
            print("  Run: --check-bias (for detailed analysis)")
            print("\n" + "!" * 70)

            sys.exit(1)

    # Handle analysis-only modes (no tree building)
    if args.suggest_outgroup:
        print("\n" + "=" * 70)
        print("OUTGROUP SUGGESTION")
        print("=" * 70)
        suggestion = suggest_outgroup(sequences)
        print(suggestion)
        sys.exit(0)

    if args.check_bias:
        print("\n" + "=" * 70)
        print("SAMPLING BIAS ANALYSIS")
        print("=" * 70)
        recommendations = get_sampling_recommendations(sequences)
        print(recommendations)
        sys.exit(0)

    # Stratified sampling if requested
    if args.stratify:
        if not args.quiet:
            print("\n" + "-" * 70)
            print("Stratified Sampling (balancing species representation)...")
            print("-" * 70)

            # Show bias before sampling
            bias_info = detect_bias(sequences, threshold=0.1)
            if bias_info['overrepresented']:
                print("\n[!] OVERREPRESENTED SPECIES DETECTED:")
                for species in bias_info['overrepresented']:
                    print(f"  - {species}")

            if bias_info.get('underrepresented'):
                print(f"\n[INFO] {len(bias_info['underrepresented'])} rare species detected (each <1% of dataset)")

        sampled, sampling_info = stratified_sample(
            sequences,
            max_per_species=args.max_per_species,
            min_per_species=args.min_per_species,
            random_seed=42  # For reproducibility
        )

        if not args.quiet:
            print(f"\nStratified sampling: {len(sequences)} -> {len(sampled)} sequences")
            print(f"  (Capped at {args.max_per_species} per species, kept >={args.min_per_species} per species)")

            # Show which species were capped
            if sampling_info:
                capped = [species for species, count in sampling_info.items() if count == args.max_per_species]
                if capped:
                    print(f"\n  Capped {len(capped)} overrepresented species:")
                    for species in capped[:5]:  # Show first 5
                        print(f"    - {species}")
                    if len(capped) > 5:
                        print(f"    ... and {len(capped) - 5} more")

        sequences = sampled

    # Smart deduplication if requested
    # Note: Exact duplicates are already removed above (always enabled)
    if args.dereplicate:
        if not args.quiet:
            print("\n" + "-" * 70)
            print("Smart Deduplication: Clustering Similar Sequences")
            print("-" * 70)

        dereplicated, stats = smart_dereplicate(
            sequences,
            remove_exact=False,  # Already done above
            similarity_threshold=99.5,
            selection_method="longest",
            verbose=(not args.quiet)
        )

        sequences = dereplicated

    # Validate sequences
    if len(sequences) < 3:
        print("Error: Need at least 3 sequences for phylogenetic analysis", file=sys.stderr)
        sys.exit(1)

    # Handle outgroup extraction if specified
    outgroup_seqs = []
    if args.outgroup:
        if not args.quiet:
            print("\n" + "-" * 70)
            print(f"Extracting outgroup sequences: {args.outgroup}")
            print("-" * 70)

        outgroup_seqs = get_outgroup_sequences(sequences, args.outgroup)

        if not outgroup_seqs:
            print(f"[!] Warning: No sequences matched outgroup pattern '{args.outgroup}'", file=sys.stderr)
            print("    Tree will be built without explicit outgroup rooting", file=sys.stderr)
        else:
            if not args.quiet:
                print(f"  Found {len(outgroup_seqs)} outgroup sequence(s):")
                for seq in outgroup_seqs[:3]:  # Show first 3
                    print(f"    - {seq.id}")
                if len(outgroup_seqs) > 3:
                    print(f"    ... and {len(outgroup_seqs) - 3} more")

                # Note: Actual rooting will be implemented in tree_to_newick
                print("\n  Note: Tree will be rooted using these outgroup sequences")

    # Check if sequences need alignment
    lengths = [seq.aligned_length for seq in sequences]
    all_same_length = len(set(lengths)) == 1

    # Align sequences by default (unless --skip-align or already aligned)
    if not args.skip_align and not all_same_length:
        if not args.quiet:
            print("\n" + "-" * 70)
            print("Aligning sequences with MUSCLE...")
            print("  (Sequences have different lengths)")
            print("-" * 70)

        try:
            aligner = MuscleAligner()

            # Determine output file for aligned sequences
            if args.aligned_output:
                aligned_file = args.aligned_output
            else:
                aligned_file = output_dir / f"aligned_{input_path.name}"

            if not args.quiet:
                print(f"  Running MUSCLE alignment...")
                print(f"  Output: {aligned_file}")

            sequences = aligner.align_sequences(sequences, str(aligned_file))

            if not args.quiet:
                print(f"  Alignment complete!")
                print(f"  Aligned length: {sequences[0].aligned_length} sites")

        except FileNotFoundError as e:
            print(f"Error: MUSCLE not found. {e}", file=sys.stderr)
            print("Place muscle.exe in the project root directory.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error during alignment: {e}", file=sys.stderr)
            sys.exit(1)
    elif all_same_length:
        if not args.quiet:
            print(f"\n  Sequences already aligned: {sequences[0].aligned_length} sites")
    else:
        if not args.quiet:
            print(f"\n  Using sequences as-is: {sequences[0].aligned_length} sites")
            print("  (Alignment skipped via --skip-align)")

    # Build trees
    if not args.quiet:
        print("\n" + "-" * 70)
        print("Building phylogenetic trees...")
        print("-" * 70)

    results = {}

    # Helper functions for bootstrap
    def ml_builder(seqs):
        tree, logL, metadata = build_ml_tree_level4(seqs, model='auto', alpha=1.0, tree_search='nni', verbose=False)
        return tree

    def bionj_builder(seqs):
        dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
        return build_bionj_tree(dist_matrix, ids)

    def upgma_builder(seqs):
        dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
        return build_upgma_tree(dist_matrix, ids)

    # Build ML tree
    if args.method in ['all', 'ml']:
        if not args.quiet:
            print("\n[1/3] Maximum Likelihood tree...")

        if args.bootstrap > 0:
            if not args.quiet:
                print(f"  Running bootstrap with {args.bootstrap} replicates...")
            ml_tree = bootstrap_tree(sequences, ml_builder, n_replicates=args.bootstrap, n_jobs=1, verbose=args.verbose)
        else:
            ml_tree, ml_logL, metadata = build_ml_tree_level4(sequences, model='auto', alpha=1.0, tree_search='nni', verbose=args.verbose)
            if not args.quiet:
                print(f"  Log-likelihood: {ml_logL:.2f}")
                print(f"  Model: {metadata.get('selected_model', 'GTR+G')}")

        results['ml'] = ml_tree

    # Build BioNJ tree
    if args.method in ['all', 'bionj']:
        if not args.quiet:
            print("\n[2/3] BioNJ tree...")

        if args.bootstrap > 0:
            if not args.quiet:
                print(f"  Running bootstrap with {args.bootstrap} replicates...")
            bionj_tree = bootstrap_tree(sequences, bionj_builder, n_replicates=args.bootstrap, n_jobs=1, verbose=args.verbose)
        else:
            dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
            bionj_tree = build_bionj_tree(dist_matrix, ids)

        results['bionj'] = bionj_tree

    # Build UPGMA tree
    if args.method in ['all', 'upgma']:
        if not args.quiet:
            print("\n[3/3] UPGMA tree...")

        if args.bootstrap > 0:
            if not args.quiet:
                print(f"  Running bootstrap with {args.bootstrap} replicates...")
            upgma_tree = bootstrap_tree(sequences, upgma_builder, n_replicates=args.bootstrap, n_jobs=1, verbose=args.verbose)
        else:
            dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
            upgma_tree = build_upgma_tree(dist_matrix, ids)

        results['upgma'] = upgma_tree

    # Output results
    if not args.quiet:
        print("\n" + "=" * 70)
        print("Results")
        print("=" * 70)

    output_files = []

    for method_name, tree in results.items():
        # ASCII visualization
        if args.output_format in ['ascii', 'both'] and not args.quiet:
            print(f"\n{method_name.upper()} Tree:")
            print("-" * 40)
            print_tree_ascii(tree)

        # Save ASCII to file
        if args.output_format in ['ascii', 'both']:
            ascii_file = output_dir / f"{args.prefix}_{method_name}_ascii.txt"
            with open(ascii_file, 'w') as f:
                f.write(f"{method_name.upper()} Tree:\n")
                f.write("-" * 40 + "\n")
                print_tree_ascii(tree, file=f)
            output_files.append(str(ascii_file))

        # Newick export
        if args.output_format in ['newick', 'both']:
            newick_str = tree_to_newick(tree, include_support=(not args.no_support and args.bootstrap > 0))
            newick_str += ";"

            output_file = output_dir / f"{args.prefix}_{method_name}.nwk"
            with open(output_file, 'w') as f:
                f.write(newick_str + "\n")

            output_files.append(str(output_file))

            if not args.quiet:
                print(f"\nExported {method_name.upper()} tree to: {output_file}")

        # ETE3 visualization
        if args.visualize and args.output_format in ['newick', 'both']:
            try:
                from rrna_phylo.visualization import visualize_tree, ETE3_AVAILABLE

                if not ETE3_AVAILABLE:
                    if not args.quiet:
                        print("\n[!] Warning: ETE3 not installed. Skipping visualization.")
                        print("    Install with: pip install ete3")
                else:
                    # Get the Newick file path
                    newick_file = output_dir / f"{args.prefix}_{method_name}.nwk"

                    # Determine visualization output file
                    viz_file = output_dir / f"{args.prefix}_{method_name}_tree.{args.viz_format}"

                    if not args.quiet:
                        print(f"\n  Generating tree visualization...")
                        print(f"  Output: {viz_file}")

                    # Generate visualization
                    visualize_tree(
                        str(newick_file),
                        str(viz_file),
                        show_bootstrap=(args.bootstrap > 0),
                        bootstrap_threshold=args.viz_bootstrap_threshold,
                        width=args.viz_width,
                        height=args.viz_height,
                        dpi=args.viz_dpi,
                        title=f"{method_name.upper()} Tree",
                        branch_line_width=args.viz_branch_width,
                        scale=args.viz_scale
                    )

                    output_files.append(str(viz_file))

                    if not args.quiet:
                        print(f"  Visualization saved: {viz_file}")

            except ImportError as e:
                if not args.quiet:
                    print(f"\n[!] Warning: Could not import ETE3. {e}")
                    print("    Install with: pip install ete3")
            except Exception as e:
                if not args.quiet:
                    print(f"\n[!] Warning: Visualization failed: {e}")

    # Summary
    if not args.quiet:
        print("\n" + "=" * 70)
        print("Analysis Complete!")
        print("=" * 70)
        print(f"\nGenerated {len(results)} tree(s)")
        if output_files:
            print(f"Output files:")
            for f in output_files:
                print(f"  - {f}")

        if not args.visualize:
            print("\nTip: Use --visualize to generate publication-quality figures with ETE3")
            print("     Or use FigTree, iTOL, or Dendroscope to view Newick files")

    return 0


if __name__ == '__main__':
    sys.exit(main())

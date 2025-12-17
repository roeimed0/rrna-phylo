"""
Production-Ready Phylogenetic Tree Builder

This script builds a publication-quality phylogenetic tree from aligned sequences.

Features:
- Automatic model selection (JC69, K80, F81, HKY85, GTR with/without gamma)
- Tree search methods: NNI (fast) or SPR (thorough)
- Bootstrap support calculation
- Branch length optimization
- Multiple output formats (Newick, ASCII)
- Comprehensive metadata and statistics

Usage:
    python build_phylogenetic_tree.py input.fasta [--method nni|spr] [--bootstrap N]

Output:
    results/
        tree.nwk          - Newick format tree
        tree_ascii.txt    - ASCII visualization
        metadata.json     - Tree statistics and model info
        log.txt           - Build log with timing
"""

import sys
import json
import time
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))

from rrna_phylo.io.fasta_parser import FastaParser
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4
from rrna_phylo.models.bootstrap import bootstrap_analysis


def build_tree(input_file: str, method: str = 'nni', bootstrap: int = 0,
               output_dir: str = 'results'):
    """
    Build a phylogenetic tree from aligned sequences.

    Args:
        input_file: Path to aligned FASTA file
        method: Tree search method ('nni' or 'spr')
        bootstrap: Number of bootstrap replicates (0 = no bootstrap)
        output_dir: Output directory for results

    Returns:
        dict: Build results with tree, metadata, and timing
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Initialize log
    log_lines = []
    log_lines.append('='*70)
    log_lines.append('PHYLOGENETIC TREE BUILDER')
    log_lines.append('='*70)
    log_lines.append(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    log_lines.append(f'Input: {input_file}')
    log_lines.append(f'Method: {method.upper()}')
    log_lines.append(f'Bootstrap: {bootstrap if bootstrap > 0 else "None"}')
    log_lines.append('')

    # Load sequences
    log_lines.append('[1/4] Loading sequences...')
    print('\n'.join(log_lines[-2:]))

    parser = FastaParser()
    sequences = parser.parse(input_file)

    # Filter sequences (remove ambiguous bases)
    filtered = []
    for seq in sequences:
        chars = set(seq.sequence.upper()) - {'-', 'A', 'C', 'G', 'T', 'U'}
        if not chars:
            filtered.append(seq)

    n_removed = len(sequences) - len(filtered)
    log_lines.append(f'  Loaded: {len(sequences)} sequences')
    if n_removed > 0:
        log_lines.append(f'  Filtered: {n_removed} sequences (ambiguous bases)')
    log_lines.append(f'  Using: {len(filtered)} sequences')
    log_lines.append(f'  Alignment length: {len(filtered[0].sequence)} bp')
    log_lines.append('')

    print('\n'.join(log_lines[-5:]))

    # Build tree
    if bootstrap > 0:
        log_lines.append('[2/4] Building phylogenetic tree with bootstrap analysis...')
        log_lines.append(f'  Method: {method.upper()}')
        log_lines.append(f'  Bootstrap replicates: {bootstrap}')
        print('\n'.join(log_lines[-3:]))

        start_time = time.time()

        tree, boot_metadata = bootstrap_analysis(
            filtered,
            n_replicates=bootstrap,
            method=method,
            verbose=True
        )

        build_time = time.time() - start_time

        # Extract metadata for compatibility
        logL = boot_metadata['reference_logL']
        metadata = {
            'selected_model': boot_metadata['reference_model'],
            'time_model_selection': boot_metadata['reference_time'],
            'bootstrap_enabled': True,
            'n_replicates': bootstrap,
            'bootstrap_time': boot_metadata['bootstrap_time'],
            'avg_support': boot_metadata['avg_support'],
            'support_distribution': boot_metadata['support_distribution']
        }
    else:
        log_lines.append('[2/4] Building phylogenetic tree...')
        log_lines.append(f'  Method: {method.upper()}')
        print('\n'.join(log_lines[-2:]))

        start_time = time.time()

        tree, logL, metadata = build_ml_tree_level4(
            filtered,
            skip_model_selection=False,
            tree_search=method,
            verbose=True
        )

        build_time = time.time() - start_time
        metadata['bootstrap_enabled'] = False

    log_lines.append('')
    log_lines.append('[3/4] Tree build complete!')
    log_lines.append(f'  Total time: {build_time:.2f}s')
    log_lines.append(f'  Log-likelihood: {logL:.4f}')
    log_lines.append(f'  Model: {metadata["selected_model"]}')

    if metadata.get('bootstrap_enabled'):
        # Bootstrap-specific output
        log_lines.append(f'  Bootstrap replicates: {metadata["n_replicates"]}')
        log_lines.append(f'  Bootstrap time: {metadata["bootstrap_time"]:.2f}s')
        log_lines.append(f'  Average support: {metadata["avg_support"]:.1f}%')

        dist = metadata['support_distribution']
        log_lines.append(f'  Support distribution:')
        log_lines.append(f'    Strong (>=95%):    {dist["strong_95+"]:3d} clades')
        log_lines.append(f'    Moderate (70-94%): {dist["moderate_70-94"]:3d} clades')
        log_lines.append(f'    Weak (50-69%):     {dist["weak_50-69"]:3d} clades')
        log_lines.append(f'    Poor (<50%):       {dist["poor_<50"]:3d} clades')
    else:
        # Regular tree output
        if metadata.get('alpha'):
            log_lines.append(f'  Gamma alpha: {metadata["alpha"]:.3f}')
        log_lines.append(f'  Model selection: {metadata["time_model_selection"]:.2f}s')
        log_lines.append(f'  Tree search: {metadata.get("time_tree_search", 0):.2f}s')

        if method == 'nni':
            log_lines.append(f'  NNI improvements: {metadata.get("n_nni_improvements", 0)}')
        elif method == 'spr':
            log_lines.append(f'  NNI improvements: {metadata.get("n_nni_improvements", 0)}')
            log_lines.append(f'  SPR improvements: {metadata.get("n_spr_improvements", 0)}')

        log_lines.append(f'  Branches collapsed: {metadata.get("n_branches_collapsed", 0)}')

    log_lines.append('')

    print('\n'.join(log_lines[-15:]))

    # Save outputs
    log_lines.append('[4/4] Saving results...')
    print('\n'.join(log_lines[-1:]))

    # Save Newick tree
    newick_file = output_path / 'tree.nwk'
    with open(newick_file, 'w') as f:
        f.write(tree.to_newick())
    log_lines.append(f'  Saved: {newick_file}')
    print(log_lines[-1])

    # Save ASCII visualization
    ascii_file = output_path / 'tree_ascii.txt'
    try:
        # Get taxa names for display
        taxa = tree.get_leaf_names()
        with open(ascii_file, 'w') as f:
            f.write('PHYLOGENETIC TREE VISUALIZATION\n')
            f.write('='*70 + '\n\n')
            f.write(f'Taxa ({len(taxa)}):\n')
            for i, taxon in enumerate(taxa, 1):
                f.write(f'  {i:2d}. {taxon}\n')
            f.write('\n')
            f.write('Newick Format:\n')
            f.write(tree.to_newick() + '\n')
        log_lines.append(f'  Saved: {ascii_file}')
        print(log_lines[-1])
    except Exception as e:
        log_lines.append(f'  Warning: Could not create ASCII tree: {e}')
        print(log_lines[-1])

    # Save metadata
    metadata_file = output_path / 'metadata.json'
    metadata_dict = {
        'build_date': datetime.now().isoformat(),
        'input_file': str(input_file),
        'n_sequences': len(filtered),
        'alignment_length': len(filtered[0].sequence),
        'method': method,
        'log_likelihood': logL,
        'model': metadata['selected_model'],
        'gamma_alpha': metadata.get('alpha'),
        'build_time_seconds': build_time,
        'model_selection_time': metadata['time_model_selection'],
        'tree_search_time': metadata.get('time_tree_search', 0),
        'nni_improvements': metadata.get('n_nni_improvements', 0),
        'spr_improvements': metadata.get('n_spr_improvements', 0) if method == 'spr' else None,
        'branches_collapsed': metadata.get('n_branches_collapsed', 0),
        'taxa': tree.get_leaf_names()
    }

    with open(metadata_file, 'w') as f:
        json.dump(metadata_dict, f, indent=2)
    log_lines.append(f'  Saved: {metadata_file}')
    print(log_lines[-1])

    # Save log
    log_file = output_path / 'build_log.txt'
    log_lines.append('')
    log_lines.append('='*70)
    log_lines.append('BUILD COMPLETE')
    log_lines.append('='*70)
    log_lines.append(f'Results saved to: {output_path.absolute()}')
    log_lines.append('')

    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))

    print('\n'.join(log_lines[-6:]))

    return {
        'tree': tree,
        'log_likelihood': logL,
        'metadata': metadata_dict,
        'output_dir': str(output_path.absolute())
    }


def main():
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Build phylogenetic tree from aligned sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fast tree with NNI (recommended)
  python build_phylogenetic_tree.py sequences.fasta

  # Thorough tree with SPR
  python build_phylogenetic_tree.py sequences.fasta --method spr

  # Custom output directory
  python build_phylogenetic_tree.py sequences.fasta --output my_results

Output files:
  results/tree.nwk       - Newick format tree
  results/tree_ascii.txt - ASCII visualization
  results/metadata.json  - Build statistics
  results/build_log.txt  - Complete build log
        """
    )

    parser.add_argument('input', help='Aligned FASTA file')
    parser.add_argument('--method', choices=['nni', 'spr'], default='nni',
                       help='Tree search method (default: nni)')
    parser.add_argument('--bootstrap', type=int, default=0,
                       help='Number of bootstrap replicates (default: 0)')
    parser.add_argument('--output', default='results',
                       help='Output directory (default: results)')

    args = parser.parse_args()

    # Validate input
    if not Path(args.input).exists():
        print(f'Error: Input file not found: {args.input}')
        sys.exit(1)

    # Build tree
    try:
        result = build_tree(
            args.input,
            method=args.method,
            bootstrap=args.bootstrap,
            output_dir=args.output
        )

        print(f'\nTree successfully built!')
        print(f'Log-likelihood: {result["log_likelihood"]:.4f}')
        print(f'Model: {result["metadata"]["model"]}')
        print(f'\nAll files saved to: {result["output_dir"]}')

    except Exception as e:
        print(f'\nError building tree: {e}')
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

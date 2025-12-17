"""
Bootstrap Support Analysis for Phylogenetic Trees.

Bootstrap resampling (Felsenstein 1985) assesses the reliability of tree topology
by randomly resampling alignment columns with replacement and rebuilding trees.

Bootstrap values represent the percentage of replicate trees that support each clade.
Standard interpretation:
- 95-100%: Strong support
- 70-94%: Moderate support
- 50-69%: Weak support
- <50%: Not supported

Usage:
    from rrna_phylo.models.bootstrap import bootstrap_analysis

    tree_with_support, metadata = bootstrap_analysis(
        sequences,
        n_replicates=100,
        method='nni',
        verbose=True
    )
"""

from typing import List, Tuple, Dict, Optional
import random
import time
from collections import defaultdict

from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4


def resample_alignment(sequences: List[Sequence], seed: Optional[int] = None) -> List[Sequence]:
    """
    Create a bootstrap replicate by resampling alignment columns with replacement.

    Args:
        sequences: Original aligned sequences
        seed: Random seed for reproducibility (optional)

    Returns:
        List of sequences with resampled alignment (same taxa, different columns)
    """
    if seed is not None:
        random.seed(seed)

    alignment_length = len(sequences[0].sequence)
    n_sequences = len(sequences)

    # Sample column indices with replacement
    column_indices = [random.randint(0, alignment_length - 1)
                     for _ in range(alignment_length)]

    # Build resampled sequences
    resampled = []
    for seq in sequences:
        # Resample columns
        new_sequence = ''.join([seq.sequence[i] for i in column_indices])

        # Create new Sequence object
        resampled_seq = Sequence(
            id=seq.id,
            sequence=new_sequence,
            description=seq.description
        )
        resampled.append(resampled_seq)

    return resampled


def get_clades(tree: TreeNode) -> set:
    """
    Extract all clades (groups of taxa) from a tree.

    A clade is represented as a frozenset of leaf names, allowing set operations.

    Args:
        tree: Tree to extract clades from

    Returns:
        Set of frozensets, each representing a clade
    """
    clades = set()

    def traverse(node):
        if node.is_leaf():
            return frozenset([node.name])

        left_taxa = traverse(node.left) if node.left else frozenset()
        right_taxa = traverse(node.right) if node.right else frozenset()

        clade = left_taxa | right_taxa

        # Only store clades with 2+ taxa (exclude trivial clades)
        if len(clade) > 1:
            clades.add(clade)

        return clade

    traverse(tree)
    return clades


def calculate_bootstrap_support(reference_tree: TreeNode,
                                bootstrap_trees: List[TreeNode]) -> Dict[frozenset, float]:
    """
    Calculate bootstrap support values for each clade in the reference tree.

    Bootstrap support = percentage of bootstrap trees that contain the clade.

    Args:
        reference_tree: The main tree (e.g., from original data)
        bootstrap_trees: List of trees from bootstrap replicates

    Returns:
        Dict mapping clades (frozensets of taxa) to support values (0-100)
    """
    # Extract clades from reference tree
    reference_clades = get_clades(reference_tree)

    # Count how many bootstrap trees support each clade
    clade_counts = defaultdict(int)

    for boot_tree in bootstrap_trees:
        boot_clades = get_clades(boot_tree)

        for clade in reference_clades:
            if clade in boot_clades:
                clade_counts[clade] += 1

    # Calculate support percentages
    n_replicates = len(bootstrap_trees)
    support_values = {}

    for clade in reference_clades:
        support = (clade_counts[clade] / n_replicates) * 100
        support_values[clade] = support

    return support_values


def add_bootstrap_to_tree(tree: TreeNode, support_values: Dict[frozenset, float]) -> TreeNode:
    """
    Add bootstrap support values to internal nodes as node names.

    Args:
        tree: Tree to annotate
        support_values: Dict mapping clades to support values

    Returns:
        Annotated tree (modifies in place, but also returns)
    """
    def annotate(node):
        if node.is_leaf():
            return frozenset([node.name])

        # Get taxa for this clade
        left_taxa = annotate(node.left) if node.left else frozenset()
        right_taxa = annotate(node.right) if node.right else frozenset()
        clade = left_taxa | right_taxa

        # Add bootstrap support to internal node
        if len(clade) > 1 and clade in support_values:
            support = support_values[clade]
            # Store as node name (common convention)
            node.name = f'{support:.0f}'

        return clade

    annotate(tree)
    return tree


def bootstrap_analysis(sequences: List[Sequence],
                      n_replicates: int = 100,
                      method: str = 'nni',
                      model: str = 'auto',
                      n_jobs: int = 1,
                      verbose: bool = True) -> Tuple[TreeNode, Dict]:
    """
    Perform bootstrap analysis to assess tree reliability.

    Algorithm:
    1. Build reference tree from original data
    2. Generate N bootstrap replicates by resampling alignment columns
    3. Build tree for each replicate
    4. Count how often each clade appears in replicate trees
    5. Annotate reference tree with support values

    Args:
        sequences: Aligned sequences
        n_replicates: Number of bootstrap replicates (default: 100, publication: 1000)
        method: Tree search method ('nni' or 'spr')
        model: Substitution model ('auto' or specific model)
        n_jobs: Number of parallel jobs (currently sequential, future enhancement)
        verbose: Print progress

    Returns:
        (tree_with_support, metadata)
        - tree_with_support: Reference tree with bootstrap values on internal nodes
        - metadata: Dict with statistics and timing
    """
    if verbose:
        print('='*70)
        print('BOOTSTRAP ANALYSIS')
        print('='*70)
        print(f'Sequences: {len(sequences)}')
        print(f'Alignment length: {len(sequences[0].sequence)} bp')
        print(f'Bootstrap replicates: {n_replicates}')
        print(f'Method: {method.upper()}')
        print(f'Model: {model}')
        print()

    metadata = {
        'n_sequences': len(sequences),
        'alignment_length': len(sequences[0].sequence),
        'n_replicates': n_replicates,
        'method': method,
        'model': model
    }

    # Step 1: Build reference tree
    if verbose:
        print('[1/4] Building reference tree from original data...')

    start_time = time.time()

    reference_tree, reference_logL, ref_meta = build_ml_tree_level4(
        sequences,
        tree_search=method,
        skip_model_selection=(model != 'auto'),
        verbose=False
    )

    ref_time = time.time() - start_time

    if verbose:
        print(f'  [OK] Reference tree built in {ref_time:.1f}s')
        print(f'       LogL: {reference_logL:.2f}')
        print(f'       Model: {ref_meta["selected_model"]}')
        print()

    metadata['reference_logL'] = reference_logL
    metadata['reference_model'] = ref_meta['selected_model']
    metadata['reference_time'] = ref_time

    # Step 2: Generate and analyze bootstrap replicates
    if verbose:
        print(f'[2/4] Generating {n_replicates} bootstrap replicates...')
        print(f'  (This may take a while...)')
        print()

    bootstrap_trees = []
    failed_replicates = 0

    start_time = time.time()

    for i in range(n_replicates):
        if verbose and (i + 1) % max(1, n_replicates // 10) == 0:
            elapsed = time.time() - start_time
            if elapsed > 0:
                rate = (i + 1) / elapsed
                remaining = (n_replicates - i - 1) / rate
                print(f'  Replicate {i+1}/{n_replicates} '
                      f'({100*(i+1)/n_replicates:.0f}%) - '
                      f'~{remaining:.0f}s remaining')
            else:
                print(f'  Replicate {i+1}/{n_replicates} '
                      f'({100*(i+1)/n_replicates:.0f}%)')

        try:
            # Resample alignment
            resampled = resample_alignment(sequences, seed=i)

            # Build tree from resampled data
            boot_tree, boot_logL, boot_meta = build_ml_tree_level4(
                resampled,
                tree_search=method,
                skip_model_selection=True,  # Use reference model for speed
                verbose=False
            )

            bootstrap_trees.append(boot_tree)

        except Exception as e:
            if verbose:
                print(f'  [WARNING] Replicate {i+1} failed: {e}')
            failed_replicates += 1

    boot_time = time.time() - start_time

    if verbose:
        print()
        print(f'  [OK] Bootstrap replicates complete in {boot_time:.1f}s')
        print(f'       Success: {len(bootstrap_trees)}/{n_replicates}')
        if failed_replicates > 0:
            print(f'       Failed: {failed_replicates}')
        print(f'       Average time per replicate: {boot_time/n_replicates:.1f}s')
        print()

    metadata['bootstrap_time'] = boot_time
    metadata['n_success'] = len(bootstrap_trees)
    metadata['n_failed'] = failed_replicates
    metadata['avg_time_per_replicate'] = boot_time / n_replicates

    # Step 3: Calculate bootstrap support
    if verbose:
        print('[3/4] Calculating bootstrap support values...')

    support_values = calculate_bootstrap_support(reference_tree, bootstrap_trees)

    # Get statistics
    support_list = list(support_values.values())
    avg_support = sum(support_list) / len(support_list) if support_list else 0
    min_support = min(support_list) if support_list else 0
    max_support = max(support_list) if support_list else 0

    # Count support levels
    strong = sum(1 for s in support_list if s >= 95)
    moderate = sum(1 for s in support_list if 70 <= s < 95)
    weak = sum(1 for s in support_list if 50 <= s < 70)
    poor = sum(1 for s in support_list if s < 50)

    if verbose:
        print(f'  [OK] Support calculated for {len(support_values)} clades')
        print(f'       Average: {avg_support:.1f}%')
        print(f'       Range: {min_support:.1f}% - {max_support:.1f}%')
        print()
        print(f'  Support distribution:')
        print(f'    Strong (>=95%):   {strong:3d} clades')
        print(f'    Moderate (70-94%): {moderate:3d} clades')
        print(f'    Weak (50-69%):   {weak:3d} clades')
        print(f'    Poor (<50%):     {poor:3d} clades')
        print()

    metadata['n_clades'] = len(support_values)
    metadata['avg_support'] = avg_support
    metadata['min_support'] = min_support
    metadata['max_support'] = max_support
    metadata['support_distribution'] = {
        'strong_95+': strong,
        'moderate_70-94': moderate,
        'weak_50-69': weak,
        'poor_<50': poor
    }

    # Step 4: Annotate tree
    if verbose:
        print('[4/4] Annotating tree with bootstrap values...')

    annotated_tree = add_bootstrap_to_tree(reference_tree, support_values)

    total_time = ref_time + boot_time

    if verbose:
        print(f'  [OK] Tree annotated')
        print()
        print('='*70)
        print('BOOTSTRAP ANALYSIS COMPLETE')
        print('='*70)
        print(f'Total time: {total_time:.1f}s ({total_time/60:.1f} min)')
        print(f'Reference tree: {ref_time:.1f}s')
        print(f'Bootstrap replicates: {boot_time:.1f}s')
        print()

    metadata['total_time'] = total_time

    return annotated_tree, metadata


if __name__ == '__main__':
    # Example usage
    from rrna_phylo.io.fasta_parser import FastaParser

    print('Loading test data...')
    parser = FastaParser()
    sequences = parser.parse('../test_data/bacterial_50_aligned.fasta')

    # Filter sequences
    filtered = []
    for seq in sequences:
        chars = set(seq.sequence.upper()) - {'-', 'A', 'C', 'G', 'T', 'U'}
        if not chars:
            filtered.append(seq)

    # Use small subset for quick test
    test_seqs = filtered[:10]

    print(f'\nTesting bootstrap with {len(test_seqs)} sequences...')

    # Run bootstrap with small number of replicates for testing
    tree, metadata = bootstrap_analysis(
        test_seqs,
        n_replicates=10,
        method='nni',
        verbose=True
    )

    print('\nTree with bootstrap support:')
    print(tree.to_newick())

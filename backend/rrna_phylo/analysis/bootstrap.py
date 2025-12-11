"""Bootstrap analysis for phylogenetic branch support assessment."""

import numpy as np
from typing import List, Callable, Tuple, Dict, Optional, Set
from multiprocessing import Pool, cpu_count
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode


def resample_alignment(sequences: List[Sequence], seed: Optional[int] = None) -> List[Sequence]:
    """Resample alignment columns with replacement (bootstrap sample)."""
    if seed is not None:
        np.random.seed(seed)

    n_sites = len(sequences[0].sequence)

    # Randomly sample column indices WITH replacement
    sampled_indices = np.random.choice(n_sites, size=n_sites, replace=True)

    # Create resampled sequences
    resampled_seqs = []
    for seq in sequences:
        # Build new sequence from sampled columns
        resampled_sequence = ''.join(seq.sequence[i] for i in sampled_indices)
        resampled_seqs.append(
            Sequence(seq.id, seq.description, resampled_sequence)
        )

    return resampled_seqs


def extract_bipartitions(tree: TreeNode) -> Set[frozenset]:
    """Extract all bipartitions (splits) from a tree."""
    bipartitions = set()

    # Get all leaf names (taxa)
    all_taxa = frozenset(tree.get_leaf_names())

    def get_bipartitions_recursive(node: TreeNode):
        """Recursively extract bipartitions from subtree."""
        if node.is_leaf():
            # Leaf contributes its own name
            return frozenset([node.name])

        # Get taxa in left and right subtrees
        left_taxa = get_bipartitions_recursive(node.left) if node.left else frozenset()
        right_taxa = get_bipartitions_recursive(node.right) if node.right else frozenset()

        # Current node's taxa = union of children
        current_taxa = left_taxa | right_taxa

        # Skip root (trivial split of all taxa)
        if len(current_taxa) < len(all_taxa):
            # Store bipartition (use smaller set to avoid duplicates)
            complement = all_taxa - current_taxa
            bipartition = current_taxa if len(current_taxa) <= len(complement) else complement
            bipartitions.add(bipartition)

        return current_taxa

    get_bipartitions_recursive(tree)
    return bipartitions


def calculate_bootstrap_support(
    original_tree: TreeNode,
    bootstrap_trees: List[TreeNode]
) -> TreeNode:
    """Calculate bootstrap support values for branches in original tree."""
    # Extract bipartitions from original tree
    original_bipartitions = extract_bipartitions(original_tree)

    # Count how many bootstrap trees contain each bipartition
    bipartition_counts = {bp: 0 for bp in original_bipartitions}

    for boot_tree in bootstrap_trees:
        boot_bipartitions = extract_bipartitions(boot_tree)
        for bp in original_bipartitions:
            if bp in boot_bipartitions:
                bipartition_counts[bp] += 1

    # Calculate support percentages
    n_replicates = len(bootstrap_trees)
    bipartition_support = {
        bp: (count / n_replicates) * 100
        for bp, count in bipartition_counts.items()
    }

    # Annotate original tree with support values
    def annotate_tree(node: TreeNode) -> frozenset:
        """Recursively annotate tree nodes with support values."""
        if node.is_leaf():
            return frozenset([node.name])

        # Get taxa in this subtree
        left_taxa = annotate_tree(node.left) if node.left else frozenset()
        right_taxa = annotate_tree(node.right) if node.right else frozenset()
        current_taxa = left_taxa | right_taxa

        # Get all taxa
        all_taxa = frozenset(original_tree.get_leaf_names())

        # Skip root (no support value for trivial split)
        # AND skip nodes with only one descendant
        if len(current_taxa) > 1 and len(current_taxa) < len(all_taxa):
            # Find this bipartition's support value
            complement = all_taxa - current_taxa
            bipartition = current_taxa if len(current_taxa) <= len(complement) else complement

            if bipartition in bipartition_support:
                node.support = bipartition_support[bipartition]
            else:
                # This shouldn't happen, but set to 0 as fallback
                node.support = 0.0
        else:
            # Root or single-descendant node - no support value
            node.support = None

        return current_taxa

    annotate_tree(original_tree)
    return original_tree


def _bootstrap_worker(args: Tuple) -> TreeNode:
    """Worker function for parallel bootstrap."""
    sequences, tree_builder, seed = args

    # Resample alignment
    resampled_seqs = resample_alignment(sequences, seed=seed)

    # Build tree from resampled data
    tree = tree_builder(resampled_seqs)

    return tree


def bootstrap_tree(
    sequences: List[Sequence],
    tree_builder: Callable[[List[Sequence]], TreeNode],
    n_replicates: int = 100,
    n_jobs: int = -1,
    verbose: bool = True,
    random_seed: Optional[int] = None
) -> TreeNode:
    """Build tree with bootstrap support values."""
    if verbose:
        print("=" * 70)
        print("BOOTSTRAP ANALYSIS")
        print("=" * 70)
        print(f"Sequences: {len(sequences)} taxa, {len(sequences[0].sequence)} sites")
        print(f"Replicates: {n_replicates}")
        print(f"Parallel jobs: {n_jobs if n_jobs > 0 else cpu_count()}")
        print()

    if random_seed is not None:
        np.random.seed(random_seed)

    original_tree = tree_builder(sequences)

    if verbose:
        print(f"Generating {n_replicates} bootstrap replicates...")

    if n_jobs == -1:
        n_jobs = cpu_count()
    elif n_jobs <= 0:
        n_jobs = 1

    seeds = [random_seed + i if random_seed is not None else None
             for i in range(n_replicates)]
    worker_args = [(sequences, tree_builder, seed) for seed in seeds]

    if n_jobs == 1 or n_replicates < 10:
        if verbose:
            print("Running sequential bootstrap...")
        bootstrap_trees = []
        for i, args in enumerate(worker_args):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Replicate {i + 1}/{n_replicates}")
            tree = _bootstrap_worker(args)
            bootstrap_trees.append(tree)
    else:
        if verbose:
            print(f"Running parallel bootstrap on {n_jobs} cores...")
        with Pool(processes=n_jobs) as pool:
            bootstrap_trees = pool.map(_bootstrap_worker, worker_args)

    if verbose:
        print("\nCalculating branch support values...")
    annotated_tree = calculate_bootstrap_support(original_tree, bootstrap_trees)

    if verbose:
        print("\n" + "=" * 70)
        print("BOOTSTRAP COMPLETE")
        print("=" * 70)

        support_values = []
        def collect_support(node):
            if not node.is_leaf() and hasattr(node, 'support') and node.support is not None:
                support_values.append(node.support)
            if node.left:
                collect_support(node.left)
            if node.right:
                collect_support(node.right)
        collect_support(annotated_tree)

        if support_values:
            print(f"\nSupport value statistics:")
            print(f"  Mean:   {np.mean(support_values):.1f}%")
            print(f"  Median: {np.median(support_values):.1f}%")
            print(f"  Min:    {np.min(support_values):.1f}%")
            print(f"  Max:    {np.max(support_values):.1f}%")

            strong = sum(1 for s in support_values if s >= 95)
            moderate = sum(1 for s in support_values if 75 <= s < 95)
            weak = sum(1 for s in support_values if s < 75)

            print(f"\nSupport categories:")
            print(f"  Strong (>=95%):    {strong}/{len(support_values)}")
            print(f"  Moderate (75-94%): {moderate}/{len(support_values)}")
            print(f"  Weak (<75%):       {weak}/{len(support_values)}")

        print("=" * 70)

    return annotated_tree


def bootstrap_consensus(
    sequences: List[Sequence],
    tree_builders: List[Tuple[str, Callable]],
    n_replicates: int = 100,
    n_jobs: int = -1,
    verbose: bool = True
) -> Dict[str, TreeNode]:
    """Bootstrap analysis for multiple tree building methods."""
    results = {}

    for method_name, tree_builder in tree_builders:
        if verbose:
            print(f"\n{'=' * 70}")
            print(f"BOOTSTRAP: {method_name}")
            print(f"{'=' * 70}\n")

        tree = bootstrap_tree(
            sequences,
            tree_builder,
            n_replicates=n_replicates,
            n_jobs=n_jobs,
            verbose=verbose
        )

        results[method_name] = tree

    return results

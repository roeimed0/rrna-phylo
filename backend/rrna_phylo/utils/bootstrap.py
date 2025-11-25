"""
Bootstrap analysis for phylogenetic trees.

Implements non-parametric bootstrap to assess branch support:
1. Resample alignment columns with replacement
2. Build replicate trees from resampled data
3. Calculate support values (% of replicates containing each split)

Bootstrap is the gold standard for phylogenetic confidence assessment.
"""

import numpy as np
from typing import List, Callable, Tuple, Dict, Optional, Set
from multiprocessing import Pool, cpu_count
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode


def resample_alignment(sequences: List[Sequence], seed: Optional[int] = None) -> List[Sequence]:
    """
    Resample alignment columns with replacement (bootstrap sample).

    Creates a new alignment by randomly sampling columns from the original
    alignment WITH replacement. This is the core of bootstrap methodology.

    Args:
        sequences: Original aligned sequences
        seed: Random seed for reproducibility (optional)

    Returns:
        New sequences with resampled columns

    Example:
        Original (6 sites):
        Seq A: ATGCAT
        Seq B: ATGCAT
        Seq C: GCATGC

        Bootstrap sample (randomly pick 6 sites, e.g., [1,1,3,5,2,4]):
        Seq A: AAGTAT
        Seq B: AAGTAT
        Seq C: GGCAGC
    """
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
    """
    Extract all bipartitions (splits) from a tree.

    A bipartition is a division of taxa into two groups defined by an internal branch.

    Args:
        tree: Phylogenetic tree

    Returns:
        Set of bipartitions, where each bipartition is a frozenset of taxon names

    Example:
        Tree: ((A,B),(C,D))
        Bipartitions:
        - {A, B} | {C, D}
        - {A} | {B, C, D}  (or equivalently {B, C, D})
        - {C} | {A, B, D}
        - {D} | {A, B, C}
    """
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
    """
    Calculate bootstrap support values for branches in original tree.

    For each branch in the original tree, count how many bootstrap replicates
    contain the same split, then divide by total number of replicates.

    Args:
        original_tree: Reference tree to annotate
        bootstrap_trees: List of bootstrap replicate trees

    Returns:
        Original tree with support values added to internal nodes

    Note:
        Modifies original_tree in place by setting support values
    """
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
    """
    Worker function for parallel bootstrap.

    Args:
        args: Tuple of (sequences, tree_builder, seed)

    Returns:
        Bootstrap replicate tree
    """
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
    """
    Build tree with bootstrap support values.

    This is the main function for bootstrap analysis:
    1. Build original tree
    2. Generate n_replicates bootstrap samples
    3. Build tree for each bootstrap sample (in parallel)
    4. Calculate support for each branch

    Args:
        sequences: Aligned sequences
        tree_builder: Function that takes sequences and returns a tree
        n_replicates: Number of bootstrap replicates (default 100, recommend 1000 for publication)
        n_jobs: Number of parallel jobs (-1 = all CPUs, 1 = sequential)
        verbose: Print progress
        random_seed: Seed for reproducibility

    Returns:
        Tree with bootstrap support values on internal branches

    Example:
        >>> from rrna_phylo import Sequence
        >>> from rrna_phylo.methods.bionj import build_bionj_tree
        >>> from rrna_phylo.distance.distance import calculate_distance_matrix
        >>>
        >>> # Define tree builder function
        >>> def my_tree_builder(seqs):
        ...     dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
        ...     return build_bionj_tree(dist_matrix, ids)
        >>>
        >>> # Run bootstrap
        >>> tree = bootstrap_tree(sequences, my_tree_builder, n_replicates=100)
        >>>
        >>> # Tree now has support values
        >>> print(tree.support)  # Bootstrap percentage
    """
    if verbose:
        print("=" * 70)
        print("BOOTSTRAP ANALYSIS")
        print("=" * 70)
        print(f"Sequences: {len(sequences)} taxa, {len(sequences[0].sequence)} sites")
        print(f"Replicates: {n_replicates}")
        print(f"Parallel jobs: {n_jobs if n_jobs > 0 else cpu_count()}")
        print()

    # Set random seed
    if random_seed is not None:
        np.random.seed(random_seed)

    # Step 1: Build original tree
    if verbose:
        print("Building original tree...")
    original_tree = tree_builder(sequences)

    # Step 2: Generate bootstrap replicates
    if verbose:
        print(f"Generating {n_replicates} bootstrap replicates...")

    # Determine number of jobs
    if n_jobs == -1:
        n_jobs = cpu_count()
    elif n_jobs <= 0:
        n_jobs = 1

    # Prepare arguments for workers (with different seeds for each replicate)
    seeds = [random_seed + i if random_seed is not None else None
             for i in range(n_replicates)]
    worker_args = [(sequences, tree_builder, seed) for seed in seeds]

    # Step 3: Build bootstrap trees (parallel if n_jobs > 1)
    if n_jobs == 1 or n_replicates < 10:
        # Sequential execution
        if verbose:
            print("Running sequential bootstrap...")
        bootstrap_trees = []
        for i, args in enumerate(worker_args):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Replicate {i + 1}/{n_replicates}")
            tree = _bootstrap_worker(args)
            bootstrap_trees.append(tree)
    else:
        # Parallel execution
        if verbose:
            print(f"Running parallel bootstrap on {n_jobs} cores...")
        with Pool(processes=n_jobs) as pool:
            bootstrap_trees = pool.map(_bootstrap_worker, worker_args)

    # Step 4: Calculate support values
    if verbose:
        print("\nCalculating branch support values...")
    annotated_tree = calculate_bootstrap_support(original_tree, bootstrap_trees)

    if verbose:
        print("\n" + "=" * 70)
        print("BOOTSTRAP COMPLETE")
        print("=" * 70)

        # Show support value distribution
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

            # Count strong support
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
    """
    Bootstrap analysis for multiple tree building methods.

    Runs bootstrap for each method and returns trees with support values.

    Args:
        sequences: Aligned sequences
        tree_builders: List of (name, builder_function) tuples
        n_replicates: Number of bootstrap replicates per method
        n_jobs: Number of parallel jobs
        verbose: Print progress

    Returns:
        Dictionary mapping method names to bootstrap trees

    Example:
        >>> builders = [
        ...     ('UPGMA', lambda s: build_upgma_tree(...)),
        ...     ('BioNJ', lambda s: build_bionj_tree(...)),
        ...     ('ML', lambda s: build_ml_tree(...))
        ... ]
        >>> trees = bootstrap_consensus(sequences, builders, n_replicates=100)
        >>> upgma_tree = trees['UPGMA']
    """
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

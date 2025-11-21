"""
Test bootstrap analysis implementation.

Tests resampling, support calculation, and parallel execution.
"""

from rrna_phylo import Sequence
from rrna_phylo.utils.bootstrap import (
    resample_alignment,
    extract_bipartitions,
    calculate_bootstrap_support,
    bootstrap_tree
)
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.utils import print_tree_ascii


# Define tree builders at module level (needed for multiprocessing pickle)
def bionj_builder_for_test(seqs):
    """Build BioNJ tree (module-level function for pickling)."""
    dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
    return build_bionj_tree(dist_matrix, ids)


def upgma_builder_for_test(seqs):
    """Build UPGMA tree (module-level function for pickling)."""
    dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
    return build_upgma_tree(dist_matrix, ids)


# Test sequences - primate mitochondrial DNA
primate_seqs = [
    Sequence(
        'human',
        'Human (Homo sapiens)',
        'ATGCTACCGGGCCCATACCCCAAACATGTTGGTTATACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTACTCTATTATCTTTATAACCGTAATATTCGGAACCCTTATCACACT'
        'ATCAAGCTCCCACTGATTATTCGCCTGAGCAGGCCTAGAAATAAACACTCTAGCTAT'
    ),
    Sequence(
        'chimp',
        'Chimpanzee (Pan troglodytes)',
        'ATGCTGCCGGGCCCATGCCCCAAACATGTTGGTCATACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTACTCTATCATCTTTACAACCGTTATATTCGGAACCCTTATCACGCT'
        'ATCAAGCTCCCACTGATTATCCGCCTGAGCAGGCCTAGAAATAAACACTCTAGCCAT'
    ),
    Sequence(
        'gorilla',
        'Gorilla (Gorilla gorilla)',
        'ATGCTGCCGGGCCCATACCCCAAACATGTTGGTTACACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTATCCTATCATCTTTACAACCGTTGTATTCGGAACCTTTATCACACT'
        'GTCAAGCTCCCACTGATTATTCGCCTGAGCGGGCCTAGAAATAAATACTCTAGCCAT'
    ),
    Sequence(
        'orangutan',
        'Orangutan (Pongo pygmaeus)',
        'ATGCTGCCAGGCCCACACCCCCAGCATGTTGGTCACACCCCTTCCCGTACTAATAAA'
        'CCTCACCATCTACTCTATCATCTTCACAACCGTAATATTCGGAACCCTCATTACGCT'
        'ATCAAGTTCCCACTGATTGTTCGCCTGAGCAGGACTAGAAATTAATACTCTAGCTAT'
    ),
]


def test_resample_alignment():
    """Test bootstrap resampling of alignment."""
    print("\n" + "=" * 70)
    print("TEST 1: ALIGNMENT RESAMPLING")
    print("=" * 70)

    # Create simple test alignment
    test_seqs = [
        Sequence('A', 'Seq A', 'ATGC'),
        Sequence('B', 'Seq B', 'ATGC'),
        Sequence('C', 'Seq C', 'TGCA'),
    ]

    print("\nOriginal alignment:")
    for seq in test_seqs:
        print(f"  {seq.id}: {seq.sequence}")

    # Resample with fixed seed for reproducibility
    resampled = resample_alignment(test_seqs, seed=42)

    print("\nResampled alignment (seed=42):")
    for seq in resampled:
        print(f"  {seq.id}: {seq.sequence}")

    # Verify
    assert len(resampled) == len(test_seqs), "Should have same number of sequences"
    assert len(resampled[0].sequence) == len(test_seqs[0].sequence), \
        "Should have same alignment length"
    assert all(seq.id == orig.id for seq, orig in zip(resampled, test_seqs)), \
        "Should preserve sequence IDs"

    # Test multiple resamples give different results
    resample2 = resample_alignment(test_seqs, seed=100)
    assert resample2[0].sequence != resampled[0].sequence, \
        "Different seeds should give different resamples"

    print("\n[OK] Resampling test passed")


def test_extract_bipartitions():
    """Test bipartition extraction from tree."""
    print("\n" + "=" * 70)
    print("TEST 2: BIPARTITION EXTRACTION")
    print("=" * 70)

    # Build a simple tree
    dist_matrix, ids = calculate_distance_matrix(primate_seqs, model="jukes-cantor")
    tree = build_bionj_tree(dist_matrix, ids)

    print("\nTree structure:")
    print_tree_ascii(tree, show_support=False)

    # Extract bipartitions
    bipartitions = extract_bipartitions(tree)

    print(f"\nExtracted {len(bipartitions)} bipartitions:")
    for bp in sorted(bipartitions, key=lambda x: len(x)):
        taxa = sorted(list(bp))
        print(f"  {{{', '.join(taxa)}}}")

    # Verify
    assert len(bipartitions) > 0, "Should extract at least one bipartition"
    assert all(isinstance(bp, frozenset) for bp in bipartitions), \
        "All bipartitions should be frozensets"

    # For 4 taxa, should have 2 internal branches = 2 bipartitions
    # (excluding trivial splits)
    assert len(bipartitions) == 2, f"Should have 2 bipartitions for 4 taxa, got {len(bipartitions)}"

    print("\n[OK] Bipartition extraction test passed")


def test_bootstrap_support():
    """Test bootstrap support calculation."""
    print("\n" + "=" * 70)
    print("TEST 3: BOOTSTRAP SUPPORT CALCULATION")
    print("=" * 70)

    # Build original tree
    dist_matrix, ids = calculate_distance_matrix(primate_seqs, model="jukes-cantor")
    original_tree = build_bionj_tree(dist_matrix, ids)

    # Generate a few bootstrap trees (manually for testing)
    print("\nGenerating 10 bootstrap replicates...")
    bootstrap_trees = []
    for i in range(10):
        resampled = resample_alignment(primate_seqs, seed=i)
        dist_mat, ids_boot = calculate_distance_matrix(resampled, model="jukes-cantor")
        boot_tree = build_bionj_tree(dist_mat, ids_boot)
        bootstrap_trees.append(boot_tree)

    # Calculate support
    print("Calculating support values...")
    annotated_tree = calculate_bootstrap_support(original_tree, bootstrap_trees)

    # Verify
    support_values = []
    def collect_support(node):
        if not node.is_leaf() and hasattr(node, 'support') and node.support is not None:
            support_values.append(node.support)
        if node.left:
            collect_support(node.left)
        if node.right:
            collect_support(node.right)

    collect_support(annotated_tree)

    print(f"\nFound {len(support_values)} support values:")
    for i, val in enumerate(support_values):
        print(f"  Branch {i+1}: {val:.1f}%")

    assert len(support_values) > 0, "Should calculate support for internal branches"
    assert all(0 <= s <= 100 for s in support_values), \
        "Support values should be between 0 and 100"

    print("\n[OK] Bootstrap support test passed")


def test_bootstrap_tree_sequential():
    """Test full bootstrap analysis (sequential)."""
    print("\n" + "=" * 70)
    print("TEST 4: FULL BOOTSTRAP (SEQUENTIAL)")
    print("=" * 70)

    # Run bootstrap (small number for speed)
    print("\nRunning bootstrap with 20 replicates (sequential)...")
    tree = bootstrap_tree(
        primate_seqs,
        bionj_builder_for_test,
        n_replicates=20,
        n_jobs=1,  # Sequential
        verbose=True,
        random_seed=42
    )

    print("\n\nResulting tree with bootstrap support:")
    print_tree_ascii(tree)

    # Verify
    assert tree is not None, "Should return a tree"

    # Count support values
    support_count = 0
    def count_support(node):
        nonlocal support_count
        if not node.is_leaf() and hasattr(node, 'support') and node.support is not None:
            support_count += 1
        if node.left:
            count_support(node.left)
        if node.right:
            count_support(node.right)

    count_support(tree)
    print(f"\n{support_count} branches have bootstrap support values")
    assert support_count > 0, "Should have support values on internal branches"

    print("\n[OK] Sequential bootstrap test passed")


def test_bootstrap_tree_parallel():
    """Test full bootstrap analysis (parallel)."""
    print("\n" + "=" * 70)
    print("TEST 5: FULL BOOTSTRAP (PARALLEL)")
    print("=" * 70)

    # Run bootstrap in parallel
    print("\nRunning bootstrap with 20 replicates (parallel)...")
    tree = bootstrap_tree(
        primate_seqs,
        bionj_builder_for_test,
        n_replicates=20,
        n_jobs=-1,  # Use all cores
        verbose=True,
        random_seed=42
    )

    print("\n\nResulting tree with bootstrap support:")
    print_tree_ascii(tree)

    # Verify
    assert tree is not None, "Should return a tree"

    # Collect support values
    support_values = []
    def collect_support(node):
        if not node.is_leaf() and hasattr(node, 'support') and node.support is not None:
            support_values.append(node.support)
        if node.left:
            collect_support(node.left)
        if node.right:
            collect_support(node.right)

    collect_support(tree)

    # Check reproducibility (same seed should give same result)
    tree2 = bootstrap_tree(
        primate_seqs,
        bionj_builder_for_test,
        n_replicates=20,
        n_jobs=-1,
        verbose=False,
        random_seed=42
    )

    support_values2 = []
    collect_support(tree2)

    print(f"\n\nReproducibility test:")
    print(f"  Run 1 support values: {[f'{s:.1f}' for s in support_values]}")
    print(f"  Run 2 support values: {[f'{s:.1f}' for s in support_values2]}")

    # Support values should be identical with same seed
    assert len(support_values) == len(support_values2), "Should have same number of support values"
    for s1, s2 in zip(support_values, support_values2):
        assert abs(s1 - s2) < 0.1, f"Support values should be identical: {s1} vs {s2}"

    print("  [OK] Results are reproducible with same random seed")

    print("\n[OK] Parallel bootstrap test passed")


def test_bootstrap_comparison():
    """Compare bootstrap support across different methods."""
    print("\n" + "=" * 70)
    print("TEST 6: BOOTSTRAP METHOD COMPARISON")
    print("=" * 70)

    print("\nBootstrapping UPGMA...")
    upgma_tree = bootstrap_tree(
        primate_seqs,
        upgma_builder_for_test,
        n_replicates=50,
        n_jobs=-1,
        verbose=True
    )

    print("\n\nBootstrapping BioNJ...")
    bionj_tree = bootstrap_tree(
        primate_seqs,
        bionj_builder_for_test,
        n_replicates=50,
        n_jobs=-1,
        verbose=True
    )

    # Compare
    print("\n\n" + "=" * 70)
    print("COMPARISON")
    print("=" * 70)

    print("\nUPGMA tree with bootstrap:")
    print_tree_ascii(upgma_tree)

    print("\n\nBioNJ tree with bootstrap:")
    print_tree_ascii(bionj_tree)

    print("\n[OK] Method comparison test passed")


def main():
    """Run all bootstrap tests."""
    print("\n" + "=" * 70)
    print("BOOTSTRAP ANALYSIS - COMPREHENSIVE TESTS")
    print("=" * 70)
    print(f"Testing with {len(primate_seqs)} primate sequences")
    print(f"Alignment length: {len(primate_seqs[0].sequence)} bp")
    print("=" * 70)

    try:
        test_resample_alignment()
        test_extract_bipartitions()
        test_bootstrap_support()
        test_bootstrap_tree_sequential()
        test_bootstrap_tree_parallel()
        test_bootstrap_comparison()

        print("\n\n" + "=" * 70)
        print("ALL BOOTSTRAP TESTS PASSED!")
        print("=" * 70)
        print("\nBootstrap features verified:")
        print("  [OK] Alignment resampling with replacement")
        print("  [OK] Bipartition extraction from trees")
        print("  [OK] Support value calculation")
        print("  [OK] Sequential bootstrap execution")
        print("  [OK] Parallel bootstrap execution (multiprocessing)")
        print("  [OK] Reproducibility with random seeds")
        print("  [OK] Multiple method comparison")
        print("\nReady for publication-quality phylogenetic analysis!")
        print("=" * 70)

    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()

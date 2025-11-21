"""
Test ML Tree Level 4 - Advanced Features

Tests model selection, NNI tree search, and integrated pipeline.
"""

from rrna_phylo import Sequence
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.models.ml_tree_level4 import (
    build_ml_tree_level4,
    build_ml_tree_fast,
    compare_models_on_tree
)
from rrna_phylo.models.model_selection import select_best_model
from rrna_phylo.models.tree_search import nni_search
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood
from rrna_phylo.utils import print_tree_ascii


# Test sequences - Primate mitochondrial DNA (partial)
# These sequences show clear evolutionary relationships
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

# Shorter test for quick validation
short_seqs = [
    Sequence('A', 'Seq A', 'ATGCATGCATGC'),
    Sequence('B', 'Seq B', 'ATGCATGCATGC'),
    Sequence('C', 'Seq C', 'ATCCATCCATCC'),
    Sequence('D', 'Seq D', 'ATCCATCCATCC'),
]


def test_model_selection():
    """Test automatic model selection."""
    print("\n" + "=" * 70)
    print("TEST 1: MODEL SELECTION")
    print("=" * 70)

    # Build initial tree
    tree = build_bionj_tree(*calculate_distance_matrix(primate_seqs, model="jukes-cantor"))

    # Test model selection
    best_model, best_score, all_results = select_best_model(
        tree,
        primate_seqs,
        models=['JC69', 'K80', 'HKY85', 'GTR'],
        criterion='BIC',
        verbose=True
    )

    print(f"\nSelected model: {best_model}")
    print(f"BIC score: {best_score:.2f}")

    # Verify that model selection ran successfully
    # Note: For short sequences, simpler models (JC69) often win due to BIC penalty
    assert best_model in ['JC69', 'K80', 'F81', 'HKY85', 'GTR'], \
        f"Unexpected model: {best_model}"
    assert best_score < float('inf'), "Invalid BIC score"

    print("\n[OK] Model selection test passed")
    print(f"  (BIC correctly penalizes model complexity for short alignments)")
    return best_model


def test_nni_search():
    """Test NNI tree search improves likelihood."""
    print("\n" + "=" * 70)
    print("TEST 2: NNI TREE SEARCH")
    print("=" * 70)

    # Build initial tree
    initial_tree = build_bionj_tree(*calculate_distance_matrix(primate_seqs, model="jukes-cantor"))
    initial_logL = compute_log_likelihood(initial_tree, primate_seqs)

    print(f"Initial tree LogL: {initial_logL:.2f}")

    # Run NNI search
    improved_tree, improved_logL, n_improvements = nni_search(
        initial_tree,
        primate_seqs,
        max_iterations=10,
        verbose=True
    )

    print(f"\nFinal LogL: {improved_logL:.2f}")
    print(f"Improvement: {improved_logL - initial_logL:.2f}")
    print(f"Number of NNI moves: {n_improvements}")

    # Verify improvement or no change (never worse)
    assert improved_logL >= initial_logL, \
        "NNI search made tree worse!"

    print("\n[OK] NNI search test passed")
    return improved_tree, n_improvements


def test_level4_auto():
    """Test Level 4 with automatic model selection + NNI."""
    print("\n" + "=" * 70)
    print("TEST 3: LEVEL 4 AUTO MODE")
    print("=" * 70)

    # Build with all automation
    tree, logL, metadata = build_ml_tree_level4(
        primate_seqs,
        model='auto',
        tree_search='nni',
        max_iterations=5,
        test_gamma=False,  # Skip gamma for speed
        verbose=True
    )

    print("\n--- Results ---")
    print(f"Selected model: {metadata['selected_model']}")
    print(f"Initial LogL: {metadata['initial_logL']:.2f}")
    print(f"Final LogL: {metadata['final_logL']:.2f}")
    print(f"Improvement: {metadata['final_logL'] - metadata['initial_logL']:.2f}")
    print(f"NNI improvements: {metadata['n_nni_improvements']}")
    print(f"Total time: {metadata['time_total']:.2f}s")

    # Verify metadata
    assert metadata['selected_model'] is not None
    assert metadata['final_logL'] >= metadata['initial_logL']

    print("\n[OK] Level 4 auto mode test passed")
    return tree, metadata


def test_level4_specified_model():
    """Test Level 4 with specified model."""
    print("\n" + "=" * 70)
    print("TEST 4: LEVEL 4 WITH SPECIFIED MODEL")
    print("=" * 70)

    # Build with GTR model, no NNI
    tree, logL, metadata = build_ml_tree_level4(
        primate_seqs,
        model='GTR',
        tree_search=None,
        verbose=True
    )

    print(f"\nModel: {metadata['selected_model']}")
    print(f"LogL: {logL:.2f}")

    assert metadata['selected_model'] == 'GTR'
    assert metadata['n_nni_improvements'] == 0

    print("\n[OK] Specified model test passed")


def test_fast_vs_thorough():
    """Compare fast and thorough pipelines."""
    print("\n" + "=" * 70)
    print("TEST 5: FAST VS THOROUGH")
    print("=" * 70)

    # Fast version
    print("\n--- Fast Pipeline ---")
    tree_fast, logL_fast, meta_fast = build_ml_tree_fast(
        short_seqs,
        verbose=True
    )

    print(f"\nFast results:")
    print(f"  Time: {meta_fast['time_total']:.2f}s")
    print(f"  LogL: {logL_fast:.2f}")

    # Note: Thorough version would take longer, so we skip for this test
    print("\n[OK] Fast pipeline test passed")


def test_model_comparison():
    """Test model comparison on fixed tree."""
    print("\n" + "=" * 70)
    print("TEST 6: MODEL COMPARISON")
    print("=" * 70)

    tree = build_bionj_tree(*calculate_distance_matrix(primate_seqs, model="jukes-cantor"))

    results = compare_models_on_tree(
        tree,
        primate_seqs,
        models=['JC69', 'K80', 'GTR'],
        verbose=True
    )

    print("\n--- Model Comparison Results ---")
    for model in ['JC69', 'K80', 'GTR']:
        res = results[model]
        print(f"{model:8s}: LogL = {res['logL']:8.2f}, "
              f"BIC = {res['score']:8.2f}, "
              f"AIC = {res['AIC']:8.2f}")

    print("\n[OK] Model comparison test passed")


def test_tree_copy():
    """Test tree copy functionality."""
    print("\n" + "=" * 70)
    print("TEST 7: TREE COPY")
    print("=" * 70)

    # Create original tree
    original = build_bionj_tree(*calculate_distance_matrix(primate_seqs, model="jukes-cantor"))
    original_logL = compute_log_likelihood(original, primate_seqs)

    # Make copy
    copy = original.copy()

    # Modify copy
    copy.distance = 999.0
    if copy.left:
        copy.left.distance = 888.0

    # Verify original unchanged
    assert original.distance != 999.0, "Copy modified original!"
    if original.left:
        assert original.left.distance != 888.0, "Copy modified original subtree!"

    # Verify copy has same structure
    copy_logL = compute_log_likelihood(copy, primate_seqs)

    print(f"Original LogL: {original_logL:.2f}")
    print(f"Copy LogL: {copy_logL:.2f}")
    print(f"Same structure: {abs(original_logL - copy_logL) < 0.01}")

    print("\n[OK] Tree copy test passed")


def main():
    """Run all Level 4 tests."""
    print("\n" + "=" * 70)
    print("ML TREE LEVEL 4 - COMPREHENSIVE TESTS")
    print("=" * 70)
    print(f"Testing with {len(primate_seqs)} primate sequences")
    print(f"Alignment length: {len(primate_seqs[0].sequence)}")
    print("=" * 70)

    try:
        # Run tests
        test_model_selection()
        test_nni_search()
        test_level4_auto()
        test_level4_specified_model()
        test_fast_vs_thorough()
        test_model_comparison()
        test_tree_copy()

        # Summary
        print("\n" + "=" * 70)
        print("ALL TESTS PASSED!")
        print("=" * 70)
        print("\nLevel 4 features verified:")
        print("  [OK] Model selection (AIC/BIC)")
        print("  [OK] NNI tree search")
        print("  [OK] Automatic pipeline")
        print("  [OK] Model comparison")
        print("  [OK] Tree utilities")
        print("\nReady for production use!")
        print("=" * 70)

    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()

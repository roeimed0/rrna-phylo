"""
Visualize all tree building methods including ML Level 4.

Shows comparison of:
- UPGMA (distance-based)
- BioNJ (distance-based)
- ML Level 3 (GTR+Gamma)
- ML Level 4 (Auto model selection + NNI search)
"""

from rrna_phylo import Sequence, build_trees
from rrna_phylo.models.ml_tree_level4 import (
    build_ml_tree_level4,
    build_ml_tree_fast
)
from rrna_phylo.utils import print_tree_ascii
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood


# Test with protein sequences (cytochrome c)
protein_seqs = [
    Sequence(
        'human',
        'Human Cytochrome C',
        'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE'
    ),
    Sequence(
        'mouse',
        'Mouse Cytochrome C',
        'MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAEGYSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKSERVDLIAYLKKATNE'
    ),
    Sequence(
        'yeast',
        'Yeast Cytochrome C',
        'MGFAAGVAAAPASAATSKKKGPNWHKTGPNLHGLFGRKTGQAEGFKYTDANKSKGIVWNNETLMEYLENPKKYIPGTKMIFAGIKKKSERADLIAYLKKATSS'
    ),
    Sequence(
        'bacteria',
        'Bacterial Cytochrome C',
        'MGDIEKAKKVFKKCQCHTVVKGGKHKTGPNLHGLFGRKTGQAAGFAYTDANKSKGVTWTETTLMEYLENPKKYIPGTKMVFAGLKKAADRDDLIAYLKDATAS'
    )
]

# Primate DNA sequences for better model selection demonstration
primate_seqs = [
    Sequence(
        'human',
        'Human mitochondrial DNA',
        'ATGCTACCGGGCCCATACCCCAAACATGTTGGTTATACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTACTCTATTATCTTTATAACCGTAATATTCGGAACCCTTATCACACT'
        'ATCAAGCTCCCACTGATTATTCGCCTGAGCAGGCCTAGAAATAAACACTCTAGCTAT'
    ),
    Sequence(
        'chimp',
        'Chimpanzee mitochondrial DNA',
        'ATGCTGCCGGGCCCATGCCCCAAACATGTTGGTCATACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTACTCTATCATCTTTACAACCGTTATATTCGGAACCCTTATCACGCT'
        'ATCAAGCTCCCACTGATTATCCGCCTGAGCAGGCCTAGAAATAAACACTCTAGCCAT'
    ),
    Sequence(
        'gorilla',
        'Gorilla mitochondrial DNA',
        'ATGCTGCCGGGCCCATACCCCAAACATGTTGGTTACACCCCTTCCCGTACTAATAAA'
        'CCCCATCATCTATCCTATCATCTTTACAACCGTTGTATTCGGAACCTTTATCACACT'
        'GTCAAGCTCCCACTGATTATTCGCCTGAGCGGGCCTAGAAATAAATACTCTAGCCAT'
    ),
    Sequence(
        'orangutan',
        'Orangutan mitochondrial DNA',
        'ATGCTGCCAGGCCCACACCCCCAGCATGTTGGTCACACCCCTTCCCGTACTAATAAA'
        'CCTCACCATCTACTCTATCATCTTCACAACCGTAATATTCGGAACCCTCATTACGCT'
        'ATCAAGTTCCCACTGATTGTTCGCCTGAGCAGGACTAGAAATTAATACTCTAGCTAT'
    ),
]


def visualize_protein_trees():
    """Build and visualize trees for protein sequences."""
    print("\n" + "=" * 80)
    print("PROTEIN PHYLOGENETIC TREES - ALL METHODS")
    print("=" * 80)
    print(f"Sequences: {len(protein_seqs)} proteins (Cytochrome C)")
    print(f"Sequence length: ~100 amino acids")
    print("=" * 80)

    # Build all trees
    trees = build_trees(protein_seqs, method='all', verbose=True)

    print("\n" + "=" * 80)
    print("TREE VISUALIZATIONS")
    print("=" * 80)

    # UPGMA tree
    print("\n--- UPGMA Tree (Distance-based, assumes molecular clock) ---")
    print_tree_ascii(trees['upgma'])

    # BioNJ tree
    print("\n--- BioNJ Tree (Distance-based, no clock assumption) ---")
    print_tree_ascii(trees['bionj'])

    # ML Level 3 tree
    print("\n--- ML Level 3 Tree (GTR+Gamma model) ---")
    ml_tree, ml_logL = trees['ml']
    print(f"Log-likelihood: {ml_logL:.2f}")
    print_tree_ascii(ml_tree)

    # Consensus tree
    print("\n--- Consensus Tree (Majority-rule from all 3 methods) ---")
    print("Support values show % of methods agreeing on each split")
    print_tree_ascii(trees['consensus'])

    # Tree comparison
    print("\n" + "=" * 80)
    print("TREE COMPARISON (Robinson-Foulds Distance)")
    print("=" * 80)

    distances = trees['tree_distances']
    for comparison, stats in distances.items():
        method1, method2 = comparison.split('_vs_')
        print(f"\n{method1.upper()} vs {method2.upper()}:")
        print(f"  Similarity:       {stats['similarity']:.1%}")
        print(f"  Shared splits:    {stats['shared_splits']}/{stats['total_splits']}")
        print(f"  Identical:        {'Yes' if stats['identical'] else 'No'}")

    # Support values
    print("\n" + "=" * 80)
    print("CONSENSUS SUPPORT VALUES")
    print("=" * 80)

    support_vals = trees['support_values']
    print(f"\nTotal bipartitions: {len(support_vals)}")
    for bp, support in sorted(support_vals.items(), key=lambda x: x[1], reverse=True):
        taxa_str = " + ".join(sorted(bp))
        print(f"  {support:5.1f}%  [{taxa_str}]")


def visualize_ml_level4():
    """Demonstrate ML Level 4 with model selection and NNI search."""
    print("\n\n" + "=" * 80)
    print("ML TREE LEVEL 4 - ADVANCED FEATURES")
    print("=" * 80)
    print(f"Sequences: {len(primate_seqs)} primates (mtDNA)")
    print(f"Alignment length: {len(primate_seqs[0].sequence)} bp")
    print("=" * 80)

    # Level 4: Auto model selection + NNI search
    print("\n--- Building ML Level 4 Tree ---")
    print("Features: Automatic model selection (BIC) + NNI tree search")
    print()

    tree_level4, logL_level4, metadata = build_ml_tree_level4(
        primate_seqs,
        model='auto',
        tree_search='nni',
        max_iterations=5,
        criterion='BIC',
        test_gamma=False,
        verbose=True
    )

    print("\n" + "=" * 80)
    print("LEVEL 4 RESULTS")
    print("=" * 80)

    print("\n--- Final Tree ---")
    print_tree_ascii(tree_level4)

    print("\n--- Performance Summary ---")
    print(f"Selected model:      {metadata['selected_model']}")
    print(f"Initial LogL:        {metadata['initial_logL']:.2f}")
    print(f"Final LogL:          {metadata['final_logL']:.2f}")
    print(f"Improvement:         {metadata['final_logL'] - metadata['initial_logL']:.2f}")
    print(f"NNI improvements:    {metadata['n_nni_improvements']}")
    print(f"\nTiming breakdown:")
    print(f"  Model selection:   {metadata['time_model_selection']:.2f}s")
    print(f"  Tree search:       {metadata['time_tree_search']:.2f}s")
    print(f"  Total:             {metadata['time_total']:.2f}s")


def compare_all_ml_levels():
    """Compare ML Level 3 vs Level 4."""
    print("\n\n" + "=" * 80)
    print("ML LEVEL 3 vs LEVEL 4 COMPARISON")
    print("=" * 80)

    # Build with Level 3 (GTR model, no NNI)
    print("\n--- ML Level 3: GTR+Gamma (fixed topology) ---")
    trees_l3 = build_trees(primate_seqs, method='ml', alpha=1.0, verbose=False)
    ml3_tree, ml3_logL = trees_l3['ml']

    print(f"Model: GTR+Gamma (fixed)")
    print(f"Tree search: None (BioNJ starting tree)")
    print(f"Log-likelihood: {ml3_logL:.2f}")

    # Build with Level 4 (auto model + NNI)
    print("\n--- ML Level 4: Auto model + NNI search ---")
    ml4_tree, ml4_logL, ml4_meta = build_ml_tree_level4(
        primate_seqs,
        model='auto',
        tree_search='nni',
        max_iterations=5,
        test_gamma=False,
        verbose=False
    )

    print(f"Model: {ml4_meta['selected_model']} (auto-selected)")
    print(f"Tree search: NNI ({ml4_meta['n_nni_improvements']} improvements)")
    print(f"Log-likelihood: {ml4_logL:.2f}")

    # Comparison
    print("\n--- Comparison ---")
    improvement = ml4_logL - ml3_logL
    print(f"LogL improvement: {improvement:+.2f}")

    if improvement > 0:
        print("[OK] Level 4 found better tree topology!")
    elif improvement == 0:
        print("[=] Level 4 matched Level 3 (already optimal)")

    print(f"\nLevel 4 advantages:")
    print(f"  - Model selection: Automatically chose {ml4_meta['selected_model']}")
    print(f"  - Tree search: Made {ml4_meta['n_nni_improvements']} topology improvements")
    print(f"  - Total time: {ml4_meta['time_total']:.2f}s")


def fast_pipeline_demo():
    """Demonstrate the fast pipeline."""
    print("\n\n" + "=" * 80)
    print("FAST PIPELINE DEMO")
    print("=" * 80)
    print("Quick ML tree with reasonable defaults")
    print()

    tree, logL, metadata = build_ml_tree_fast(primate_seqs, verbose=True)

    print("\n--- Fast Pipeline Results ---")
    print_tree_ascii(tree)
    print(f"\nModel: {metadata['selected_model']}")
    print(f"LogL: {logL:.2f}")
    print(f"Time: {metadata['time_total']:.2f}s")


def main():
    """Run all visualizations."""
    print("\n" + "=" * 80)
    print("PHYLOGENETIC TREE BUILDER - COMPLETE DEMONSTRATION")
    print("=" * 80)
    print("\nThis demo shows:")
    print("  1. Traditional methods (UPGMA, BioNJ)")
    print("  2. ML Level 3 (GTR+Gamma)")
    print("  3. ML Level 4 (Model selection + NNI)")
    print("  4. Consensus trees with support values")
    print("  5. Tree comparison metrics")
    print("=" * 80)

    # Run all demonstrations
    visualize_protein_trees()
    visualize_ml_level4()
    compare_all_ml_levels()
    fast_pipeline_demo()

    # Final summary
    print("\n\n" + "=" * 80)
    print("DEMONSTRATION COMPLETE!")
    print("=" * 80)
    print("\nAvailable methods:")
    print("  [OK] UPGMA          - Fast, assumes molecular clock")
    print("  [OK] BioNJ          - Fast, no clock assumption")
    print("  [OK] ML Level 3     - GTR+Gamma model")
    print("  [OK] ML Level 4     - Auto model + NNI search")
    print("  [OK] Consensus      - Combine multiple trees")
    print("\nLevel 4 features:")
    print("  [OK] Model selection (JC69, K80, F81, HKY85, GTR)")
    print("  [OK] Information criteria (AIC, BIC)")
    print("  [OK] NNI tree search for better topologies")
    print("  [OK] Automatic parameter optimization")
    print("  [OK] Performance metrics and timing")
    print("\nReady for production phylogenetic analysis!")
    print("=" * 80)


if __name__ == "__main__":
    main()

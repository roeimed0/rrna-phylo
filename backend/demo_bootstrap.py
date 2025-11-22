"""
Quick bootstrap demonstration - shows the feature in action.

Run: python demo_bootstrap.py
"""

from rrna_phylo import Sequence
from rrna_phylo.utils.bootstrap import bootstrap_tree
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4
from rrna_phylo.utils import print_tree_ascii


# Define a tree builder function using Maximum Likelihood
def my_tree_builder(seqs):
    """Build ML tree using Level 4 (automatic model selection + NNI search)."""
    tree, logL, metadata = build_ml_tree_level4(
        seqs,
        model='auto',          # Automatic model selection
        alpha=1.0,             # Gamma rate heterogeneity
        tree_search='nni',     # NNI tree search
        verbose=False          # Silent during bootstrap
    )
    return tree


if __name__ == '__main__':
    # Example sequences - primates
    sequences = [
        Sequence('human', 'Human', 'ATGCTACCGGGCCCATACCCCAAACATGTTGGTTATACCCCTTCCCGTACTAATAAACCCCATCATCTACTCTATTATCTTTATAACCGTAATATTCGGAACCCTTATCACACT'),
        Sequence('chimp', 'Chimpanzee', 'ATGCTGCCGGGCCCATGCCCCAAACATGTTGGTCATACCCCTTCCCGTACTAATAAACCCCATCATCTACTCTATCATCTTTACAACCGTTATATTCGGAACCCTTATCACGCT'),
        Sequence('gorilla', 'Gorilla', 'ATGCTGCCGGGCCCATACCCCAAACATGTTGGTTACACCCCTTCCCGTACTAATAAACCCCATCATCTATCCTATCATCTTTACAACCGTTGTATTCGGAACCTTTATCACACT'),
        Sequence('orangutan', 'Orangutan', 'ATGCTGCCAGGCCCACACCCCCAGCATGTTGGTCACACCCCTTCCCGTACTAATAAACCTCACCATCTACTCTATCATCTTCACAACCGTAATATTCGGAACCCTCATTACGCT'),
    ]

    print("=" * 70)
    print("BOOTSTRAP DEMONSTRATION - ML LEVEL 4")
    print("=" * 70)
    print(f"\nSequences: {len(sequences)} primates")
    print(f"Alignment length: {len(sequences[0].sequence)} bp")
    print("\nMethod: Maximum Likelihood with GTR+Gamma")
    print("Features: Automatic model selection + NNI tree search")
    print("\nRunning bootstrap with 100 replicates (parallel)...")
    print("NOTE: This may take a few minutes due to ML optimization...\n")

    # Run bootstrap analysis
    tree = bootstrap_tree(
        sequences,
        my_tree_builder,
        n_replicates=100,      # Standard for exploratory analysis
        n_jobs=-1,             # Use all CPU cores
        verbose=True,
        random_seed=42         # For reproducibility
    )

    print("\n\nFinal tree with bootstrap support:")
    print("=" * 70)
    print_tree_ascii(tree)
    print("=" * 70)

    print("\nInterpretation:")
    print("  - Strong support (>=95%): Very confident in this branch")
    print("  - Moderate support (75-94%): Reasonably confident")
    print("  - Weak support (<75%): Low confidence, interpret with caution")
    print("\nThis tree is now publication-ready with confidence values!")

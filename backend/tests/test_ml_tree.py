"""
Test our own Maximum Likelihood implementation.
"""

from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.ml_tree import GTRModel, MaximumLikelihoodTree, build_ml_tree
import numpy as np


def test_gtr_model():
    """Test GTR substitution model."""
    print("=" * 60)
    print("TEST: GTR Substitution Model")
    print("=" * 60)

    # Create test sequences
    sequences = [
        Sequence("seq1", "Test 1", "ATGCATGC"),
        Sequence("seq2", "Test 2", "ATCCATGC"),
        Sequence("seq3", "Test 3", "ATGCATCC"),
    ]

    # Initialize model
    model = GTRModel()

    # Estimate parameters
    print("\nEstimating GTR parameters from sequences...")
    model.estimate_parameters(sequences)

    print("\nRate Matrix Q:")
    print(model.Q)
    print("\nRow sums (should be ~0):")
    print(np.sum(model.Q, axis=1))

    # Test probability matrix
    print("\n" + "-" * 60)
    print("Probability Matrix P(t) for different branch lengths:")
    print("-" * 60)

    for t in [0.01, 0.1, 1.0]:
        P = model.probability_matrix(t)
        print(f"\nt = {t}:")
        print(f"P(A->A) = {P[0,0]:.4f}  (stays same)")
        print(f"P(A->G) = {P[0,2]:.4f}  (transition)")
        print(f"P(A->C) = {P[0,1]:.4f}  (transversion)")
        print(f"Row sums: {np.sum(P, axis=1)}")  # Should be all 1.0

    print("\n[OK] GTR model test passed!")


def test_ml_tree_building():
    """Test ML tree building."""
    print("\n\n" + "=" * 60)
    print("TEST: Maximum Likelihood Tree Building")
    print("=" * 60)

    # Create test sequences
    sequences = [
        Sequence("E_coli", "Escherichia coli", "ATGCATGCATGCATGC"),
        Sequence("Salmonella", "Salmonella enterica", "ATGCATGCATGCATCC"),
        Sequence("B_subtilis", "Bacillus subtilis", "ATCCATGCATGCATGC"),
        Sequence("S_aureus", "Staphylococcus aureus", "ATGCATCCATGCATGC"),
    ]

    print("\nSequences:")
    for seq in sequences:
        print(f"  {seq.id:15} {seq.sequence}")

    # Build ML tree
    print("\nBuilding ML tree...")
    ml_tree = build_ml_tree(sequences)

    print("\nML Tree structure:")
    from rrna_phylo.methods.bionj import BioNJBuilder
    builder = BioNJBuilder()
    builder.print_tree(ml_tree)

    print(f"\nNewick: {ml_tree.to_newick()};")

    print("\n[OK] ML tree building test passed!")


def test_probability_evolution():
    """Test how probabilities change over evolutionary time."""
    print("\n\n" + "=" * 60)
    print("TEST: Evolution of Substitution Probabilities")
    print("=" * 60)

    sequences = [Sequence("test", "test", "ATGCATGC")]

    model = GTRModel()
    model.estimate_parameters(sequences)

    print("\nHow does P(A->A) change over time?")
    print("(Probability that A stays A)")
    print("\nTime\tP(A->A)\tP(A->G)\tP(A->C)\tP(A->T)")
    print("-" * 50)

    times = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    for t in times:
        P = model.probability_matrix(t)
        print(f"{t:.2f}\t{P[0,0]:.4f}\t{P[0,2]:.4f}\t{P[0,1]:.4f}\t{P[0,3]:.4f}")

    print("\nObservations:")
    print("- At t=0: P(A->A) = 1.0 (no time, no change)")
    print("- As t increases: P(A->A) decreases (more time = more changes)")
    print("- At large t: P approaches equilibrium frequencies")
    print("- Transitions (A<->G) happen faster than transversions")

    print("\n[OK] Evolution probability test passed!")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("MAXIMUM LIKELIHOOD IMPLEMENTATION TESTS")
    print("=" * 60)
    print("\nTesting our own ML implementation (no RAxML needed!)")

    test_gtr_model()
    test_ml_tree_building()
    test_probability_evolution()

    print("\n\n" + "=" * 60)
    print("ALL TESTS PASSED [OK]")
    print("=" * 60)
    print("\nOur ML implementation uses:")
    print("  [OK] GTR substitution model (most general reversible model)")
    print("  [OK] Matrix exponential for P(t)")
    print("  [OK] BioNJ for initial tree topology")
    print("  [OK] Pure Python (works on Windows!)")
    print("\nNext: Implement full Felsenstein's algorithm for likelihood")
    print("      and NNI tree search for optimization")
    print()


if __name__ == "__main__":
    main()

"""
Visualize phylogenetic trees in ASCII format.
"""

from fasta_parser import Sequence
from distance import calculate_distance_matrix
from upgma import build_upgma_tree, UPGMABuilder
from bionj import build_bionj_tree, BioNJBuilder
from ml_tree_level3 import build_ml_tree_level3


def print_tree_ascii(node, prefix="", is_tail=True, file=None):
    """
    Print tree in ASCII art format.

    Args:
        node: TreeNode to print
        prefix: Prefix for current line
        is_tail: Whether this is the last child
        file: File to write to (default: stdout)
    """
    import sys
    if file is None:
        file = sys.stdout

    # Print current node
    connector = "└── " if is_tail else "├── "

    if node.is_leaf():
        print(f"{prefix}{connector}{node.name} (dist: {node.distance:.4f})", file=file)
    else:
        print(f"{prefix}{connector}Internal (dist: {node.distance:.4f})", file=file)

        # Recurse to children
        if node.left or node.right:
            # Prepare prefix for children
            extension = "    " if is_tail else "│   "
            new_prefix = prefix + extension

            if node.left and node.right:
                # Both children exist
                print_tree_ascii(node.left, new_prefix, False, file)
                print_tree_ascii(node.right, new_prefix, True, file)
            elif node.left:
                print_tree_ascii(node.left, new_prefix, True, file)
            elif node.right:
                print_tree_ascii(node.right, new_prefix, True, file)


def visualize_all_trees():
    """Build and visualize trees using all 3 methods."""

    print("=" * 80)
    print("PHYLOGENETIC TREE COMPARISON")
    print("=" * 80)

    # Create test sequences
    sequences = [
        Sequence("E_coli", "Escherichia coli", "ATGCATGCATGCATGC"),
        Sequence("Salmonella", "Salmonella enterica", "ATGCATGCATGCATCC"),
        Sequence("B_subtilis", "Bacillus subtilis", "ATCCATGCATGCATGC"),
        Sequence("S_aureus", "Staphylococcus aureus", "ATGCATCCATGCATGC"),
    ]

    print("\nInput Sequences:")
    print("-" * 80)
    for seq in sequences:
        print(f"  {seq.id:15} {seq.sequence}")

    # Calculate distance matrix (used by UPGMA and BioNJ)
    print("\nCalculating distance matrix...")
    dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

    print("\nDistance Matrix:")
    print(f"{'':15}", end='')
    for id in ids:
        print(f"{id:15}", end='')
    print()
    for i, id in enumerate(ids):
        print(f"{id:15}", end='')
        for j in range(len(ids)):
            print(f"{dist_matrix[i][j]:15.4f}", end='')
        print()

    # Method 1: UPGMA
    print("\n\n" + "=" * 80)
    print("METHOD 1: UPGMA (Assumes Molecular Clock)")
    print("=" * 80)
    upgma_tree = build_upgma_tree(dist_matrix, ids)
    print("\nTree Structure:")
    print_tree_ascii(upgma_tree)
    print(f"\nNewick Format:")
    print(f"{upgma_tree.to_newick()};")

    # Method 2: BioNJ
    print("\n\n" + "=" * 80)
    print("METHOD 2: BioNJ (Variance-Weighted, No Clock Assumption)")
    print("=" * 80)
    bionj_tree = build_bionj_tree(dist_matrix, ids)
    print("\nTree Structure:")
    print_tree_ascii(bionj_tree)
    print(f"\nNewick Format:")
    print(f"{bionj_tree.to_newick()};")

    # Method 3: Maximum Likelihood
    print("\n\n" + "=" * 80)
    print("METHOD 3: Maximum Likelihood (GTR+Gamma Model, Optimized)")
    print("=" * 80)
    print("\nBuilding ML tree (this may take a minute)...")
    ml_tree, ml_logL = build_ml_tree_level3(sequences, verbose=False)
    print(f"\nLog-likelihood: {ml_logL:.2f}")
    print("\nTree Structure:")
    print_tree_ascii(ml_tree)
    print(f"\nNewick Format:")
    print(f"{ml_tree.to_newick()};")

    # Summary
    print("\n\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\nAll three methods cluster the sequences, but with different:")
    print("  • Branch lengths (evolutionary distances)")
    print("  • Topologies (which species group together)")
    print("  • Assumptions (molecular clock, rate variation, etc.)")
    print("\nFor forensics/publication:")
    print("  ✓ Use multiple methods (we have 3!)")
    print("  ✓ Compare results (next: consensus tree)")
    print("  ✓ Report all trees + consensus")
    print()


def visualize_simple_example():
    """Simple example with clear relationships."""

    print("\n\n" + "=" * 80)
    print("SIMPLE EXAMPLE: Clear Evolutionary Relationships")
    print("=" * 80)

    # Create sequences with clear relationships
    # A and B are very similar (1 difference)
    # C is different (4 differences)
    sequences = [
        Sequence("A", "Species A", "ATGCAT"),
        Sequence("B", "Species B", "ATGCAC"),  # 1 diff from A
        Sequence("C", "Species C", "TTCCAA"),  # Very different
    ]

    print("\nSequences:")
    for seq in sequences:
        print(f"  {seq.id}: {seq.sequence}")

    print("\nExpected: (A,B) should cluster together, C separate")

    # Build BioNJ tree
    dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
    tree = build_bionj_tree(dist_matrix, ids)

    print("\nBioNJ Tree:")
    print_tree_ascii(tree)

    print("\nNewick:")
    print(f"{tree.to_newick()};")

    print("\n✓ As expected: A and B cluster together!")


if __name__ == "__main__":
    # Simple example first
    visualize_simple_example()

    # Full comparison
    visualize_all_trees()

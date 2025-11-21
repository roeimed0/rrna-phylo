"""
Quick test of consensus tree functionality.
"""

from rrna_phylo import build_trees, Sequence
from rrna_phylo.utils import print_tree_ascii

# Sample protein sequences with different lengths (will be aligned)
protein_seqs = [
    Sequence('human', 'Human Cytochrome C', 'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE'),
    Sequence('mouse', 'Mouse Cytochrome C', 'MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAEGYSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKSERVDLIAYLKKATNE'),
    Sequence('yeast', 'Yeast Cytochrome C', 'MGFAAGVAAAPASAATSKKKGPNWHKTGPNLHGLFGRKTGQAEGFKYTDANKSKGIVWNNETLMEYLENPKKYIPGTKMIFAGIKKKSERADLIAYLKKATSS'),
    Sequence('bacteria', 'Bacterial Cytochrome C', 'MGDIEKAKKVFKKCQCHTVVKGGKHKTGPNLHGLFGRKTGQAAGFAYTDANKSKGVTWTETTLMEYLENPKKYIPGTKMVFAGLKKAADRDDLIAYLKDATAS')
]

print("=" * 70)
print("CONSENSUS TREE TEST")
print("=" * 70)

# Build all trees with consensus (verbose=True to see the process)
trees = build_trees(protein_seqs, verbose=True)

print("\n" + "=" * 70)
print("INDIVIDUAL TREES")
print("=" * 70)

print("\n--- UPGMA Tree ---")
print_tree_ascii(trees['upgma'])

print("\n--- BioNJ Tree ---")
print_tree_ascii(trees['bionj'])

print("\n--- ML Tree ---")
ml_tree, logL = trees['ml']
print(f"Log-likelihood: {logL:.2f}")
print_tree_ascii(ml_tree)

print("\n" + "=" * 70)
print("CONSENSUS TREE")
print("=" * 70)

consensus = trees['consensus']
print("\nConsensus tree with support values:")
print_tree_ascii(consensus)

# Print support values
print("\n" + "=" * 70)
print("SUPPORT VALUES")
print("=" * 70)

support_vals = trees['support_values']
print(f"\nTotal bipartitions with support: {len(support_vals)}")
for bp, support in sorted(support_vals.items(), key=lambda x: x[1], reverse=True):
    taxa_str = ", ".join(sorted(bp))
    print(f"  {support:5.1f}%  [{taxa_str}]")

# Print tree distances
print("\n" + "=" * 70)
print("TREE COMPARISON")
print("=" * 70)

distances = trees['tree_distances']
for comparison, stats in distances.items():
    method1, method2 = comparison.split('_vs_')
    print(f"\n{method1.upper()} vs {method2.upper()}:")
    print(f"  RF distance:      {stats['rf_distance']}")
    print(f"  RF normalized:    {stats['rf_normalized']:.3f}")
    print(f"  Similarity:       {stats['similarity']:.1%}")
    print(f"  Shared splits:    {stats['shared_splits']}/{stats['total_splits']}")
    print(f"  Identical:        {stats['identical']}")

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)

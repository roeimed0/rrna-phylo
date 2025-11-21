"""
Visualize all tree building methods for DNA, RNA, and Protein sequences.

This script demonstrates:
- UPGMA, BioNJ, and ML tree building
- For DNA, RNA, and Protein sequence types
- Complete visualization of all 9 combinations
"""

import sys
from rrna_phylo import build_trees, Sequence
from rrna_phylo.utils import print_tree_ascii

def print_section_header(title):
    """Print a formatted section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def print_subsection(title):
    """Print a formatted subsection header."""
    print(f"\n--- {title} ---")

def visualize_dna_trees():
    """Build and visualize trees for DNA sequences."""
    print_section_header("DNA SEQUENCES")

    # Sample DNA sequences (16S rRNA-like)
    dna_seqs = [
        Sequence('ecoli', 'E. coli', 'ATCGATCGATCGATCGATCG'),
        Sequence('salmonella', 'Salmonella', 'ATCGATCGATCGTTCGATCG'),
        Sequence('bacillus', 'Bacillus', 'ATCGATGGTTCGATCGATCG'),
        Sequence('streptococcus', 'Streptococcus', 'ATCGGGCGATCGATCGATCG')
    ]

    print("\nSequences:")
    for seq in dna_seqs:
        print(f"  {seq.id:15} {seq.sequence[:30]}...")

    print("\n[Building trees with UPGMA, BioNJ, and ML (GTR+Gamma)...]")
    trees = build_trees(dna_seqs, verbose=False)

    # UPGMA Tree
    print_subsection("UPGMA Tree")
    print_tree_ascii(trees['upgma'])

    # BioNJ Tree
    print_subsection("BioNJ Tree")
    print_tree_ascii(trees['bionj'])

    # ML Tree
    print_subsection("ML (GTR+Gamma) Tree")
    ml_tree, log_likelihood = trees['ml']
    print(f"Log-likelihood: {log_likelihood:.2f}")
    print_tree_ascii(ml_tree)

def visualize_rna_trees():
    """Build and visualize trees for RNA sequences."""
    print_section_header("RNA SEQUENCES")

    # Sample RNA sequences (18S rRNA-like)
    rna_seqs = [
        Sequence('human', 'Human', 'AUCGAUCGAUCGAUCGAUCG'),
        Sequence('mouse', 'Mouse', 'AUCGAUCGAUCGUUCGAUCG'),
        Sequence('rat', 'Rat', 'AUCGAUCGAUCGUUCGAUCG'),
        Sequence('yeast', 'Yeast', 'AUCGGUCGAUCGAUCGAUCG')
    ]

    print("\nSequences:")
    for seq in rna_seqs:
        print(f"  {seq.id:15} {seq.sequence[:30]}...")

    print("\n[Building trees with UPGMA, BioNJ, and ML (GTR+Gamma)...]")
    trees = build_trees(rna_seqs, verbose=False)

    # UPGMA Tree
    print_subsection("UPGMA Tree")
    print_tree_ascii(trees['upgma'])

    # BioNJ Tree
    print_subsection("BioNJ Tree")
    print_tree_ascii(trees['bionj'])

    # ML Tree
    print_subsection("ML (GTR+Gamma) Tree")
    ml_tree, log_likelihood = trees['ml']
    print(f"Log-likelihood: {log_likelihood:.2f}")
    print_tree_ascii(ml_tree)

def visualize_protein_trees():
    """Build and visualize trees for Protein sequences."""
    print_section_header("PROTEIN SEQUENCES")

    # Sample protein sequences (cytochrome c-like)
    protein_seqs = [
        Sequence('human', 'Human Cytochrome C', 'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE'),
        Sequence('mouse', 'Mouse Cytochrome C', 'MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAEGYSYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKSERVDLIAYLKKATNE'),
        Sequence('yeast', 'Yeast Cytochrome C', 'MGFAAGVAAAPASAATSKKKGPNWHKTGPNLHGLFGRKTGQAEGFKYTDANKSKGIVWNNETLMEYLENPKKYIPGTKMIFAGIKKKSERADLIAYLKKATSS'),
        Sequence('bacteria', 'Bacterial Cytochrome C', 'MGDIEKAKKVFKKCQCHTVVKGGKHKTGPNLHGLFGRKTGQAAGFAYTDANKSKGVTWTETTLMEYLENPKKYIPGTKMVFAGLKKAADRDDLIAYLKDATAS')
    ]

    print("\nSequences:")
    for seq in protein_seqs:
        print(f"  {seq.id:15} {seq.sequence[:40]}...")

    print("\n[Building trees with UPGMA, BioNJ, and ML (WAG+Gamma)...]")
    trees = build_trees(protein_seqs, verbose=False)

    # UPGMA Tree
    print_subsection("UPGMA Tree")
    print_tree_ascii(trees['upgma'])

    # BioNJ Tree
    print_subsection("BioNJ Tree")
    print_tree_ascii(trees['bionj'])

    # ML Tree
    print_subsection("ML (WAG+Gamma) Tree")
    ml_tree, log_likelihood = trees['ml']
    print(f"Log-likelihood: {log_likelihood:.2f}")
    print_tree_ascii(ml_tree)

def main():
    """Main function to visualize all tree types."""
    print("\n" + "=" * 70)
    print("  PHYLOGENETIC TREE VISUALIZATION")
    print("  Testing all methods (UPGMA, BioNJ, ML) on all sequence types")
    print("=" * 70)

    try:
        # DNA Trees
        visualize_dna_trees()

        # RNA Trees
        visualize_rna_trees()

        # Protein Trees
        visualize_protein_trees()

        # Summary
        print_section_header("SUMMARY")
        print("\n[OK] Successfully built and visualized:")
        print("  - 3 DNA trees (UPGMA, BioNJ, ML with GTR+Gamma)")
        print("  - 3 RNA trees (UPGMA, BioNJ, ML with GTR+Gamma)")
        print("  - 3 Protein trees (UPGMA, BioNJ, ML with WAG+Gamma)")
        print("  Total: 9 phylogenetic trees\n")

    except Exception as e:
        print(f"\n[ERROR] Failed to build trees: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

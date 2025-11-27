"""
Complete phylogenetic pipeline test with comprehensive visualization.

This test runs the full pipeline:
1. Load and align sequences
2. Remove exact duplicates (always)
3. Build trees in TWO modes:
   - Regular (exact duplicates only)
   - Deduplicated (smart deduplication)
4. For each mode, build THREE tree types:
   - UPGMA
   - BioNJ
   - ML (Maximum Likelihood)
5. Generate ONE COMBINED OUTPUT FILE showing all 6 trees

Output: test_full_pipeline_output/COMPLETE_TREE_COMPARISON.txt
"""

import shutil
from pathlib import Path

from rrna_phylo.io.fasta_parser import parse_fasta, assign_unique_display_names
from rrna_phylo.io.aligner import MuscleAligner
from rrna_phylo.core.builder import PhylogeneticTreeBuilder
from rrna_phylo.utils.strain_handler import remove_exact_duplicates, smart_dereplicate
from rrna_phylo.utils import print_tree_ascii
from io import StringIO


def test_full_phylogenetic_pipeline():
    """
    Complete test with single comprehensive output file showing all 6 trees.
    """

    print("=" * 80)
    print("FULL PHYLOGENETIC PIPELINE TEST")
    print("=" * 80)
    print()

    # Setup output directory
    output_dir = Path("test_full_pipeline_output")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ========================================================================
    # STEP 1: Load sequences
    # ========================================================================

    print("STEP 1: Loading sequences")
    print("-" * 80)
    fasta_file = "test_real_rrana.fasta"
    sequences = parse_fasta(fasta_file)
    print(f"  Loaded {len(sequences)} sequences from {fasta_file}")
    print()

    # ========================================================================
    # STEP 2: Assign unique display names
    # ========================================================================

    print("STEP 2: Assigning unique display names")
    print("-" * 80)
    assign_unique_display_names(sequences)
    print(f"  Display names assigned")
    print(f"  Example: {sequences[0].display_name}")
    print()

    # ========================================================================
    # STEP 3: Remove exact duplicates (ALWAYS)
    # ========================================================================

    print("STEP 3: Removing exact duplicates")
    print("-" * 80)
    original_count = len(sequences)
    sequences, dup_map = remove_exact_duplicates(sequences)
    num_exact_dups = original_count - len(sequences)
    print(f"  Original: {original_count} sequences")
    print(f"  Removed: {num_exact_dups} exact duplicates")
    print(f"  Remaining: {len(sequences)} unique sequences")
    print()

    # ========================================================================
    # STEP 4: Align sequences
    # ========================================================================

    print("STEP 4: Aligning sequences with MUSCLE")
    print("-" * 80)
    aligner = MuscleAligner()
    aligned_file = output_dir / "aligned_sequences.fasta"
    aligned_seqs = aligner.align_sequences(sequences, str(aligned_file))
    print(f"  Aligned length: {len(aligned_seqs[0].sequence)} sites")
    print(f"  Saved: {aligned_file}")
    print()

    # ========================================================================
    # MODE 1: REGULAR TREES (Exact duplicates only)
    # ========================================================================

    print("=" * 80)
    print(f"MODE 1: REGULAR TREES ({len(aligned_seqs)} sequences)")
    print("=" * 80)
    print()

    builder = PhylogeneticTreeBuilder(verbose=False)

    print("Building trees...")
    print("-" * 80)

    # UPGMA
    print("  [1/3] UPGMA tree...")
    upgma_tree = builder.build_upgma_tree(aligned_seqs)
    print("      Done!")

    # BioNJ
    print("  [2/3] BioNJ tree...")
    bionj_tree = builder.build_bionj_tree(aligned_seqs)
    print("      Done!")

    # ML
    print("  [3/3] ML tree...")
    ml_tree, ml_logl = builder.build_ml_tree(aligned_seqs)
    print("      Done!")

    print()

    # ========================================================================
    # MODE 2: DEDUPLICATED TREES (Smart deduplication)
    # ========================================================================

    print("=" * 80)
    print("MODE 2: DEDUPLICATED TREES (Smart deduplication)")
    print("=" * 80)
    print()

    print("Performing smart deduplication...")
    print("-" * 80)
    dereplicated_seqs, stats = smart_dereplicate(
        aligned_seqs,
        remove_exact=False,  # Already done above
        similarity_threshold=99.5,
        species_aware=True,
        verbose=True
    )
    print(f"  Final: {len(dereplicated_seqs)} representative sequences")
    print(f"  Reduction: {stats['reduction_percentage']:.1f}%")
    print()

    builder_derep = PhylogeneticTreeBuilder(verbose=False)

    print("Building trees...")
    print("-" * 80)

    # UPGMA
    print("  [1/3] UPGMA tree...")
    upgma_tree_derep = builder_derep.build_upgma_tree(dereplicated_seqs)
    print("      Done!")

    # BioNJ
    print("  [2/3] BioNJ tree...")
    bionj_tree_derep = builder_derep.build_bionj_tree(dereplicated_seqs)
    print("      Done!")

    # ML
    print("  [3/3] ML tree...")
    ml_tree_derep, ml_logl_derep = builder_derep.build_ml_tree(dereplicated_seqs)
    print("      Done!")

    print()

    # ========================================================================
    # GENERATE SINGLE COMBINED OUTPUT FILE
    # ========================================================================

    print("=" * 80)
    print("GENERATING COMBINED OUTPUT")
    print("=" * 80)
    print()

    output_file = output_dir / "COMPLETE_TREE_COMPARISON.txt"

    with open(output_file, 'w') as f:
        f.write("=" * 100 + "\n")
        f.write("COMPLETE PHYLOGENETIC TREE COMPARISON\n")
        f.write("=" * 100 + "\n")
        f.write("\n")
        f.write(f"Dataset: {fasta_file}\n")
        f.write(f"Original sequences: {original_count}\n")
        f.write(f"After exact duplicate removal: {len(aligned_seqs)} sequences\n")
        f.write(f"After smart deduplication: {len(dereplicated_seqs)} sequences ({stats['reduction_percentage']:.1f}% reduction)\n")
        f.write("\n")
        f.write("=" * 100 + "\n")
        f.write("MODE 1: REGULAR TREES (19 sequences - exact duplicates removed)\n")
        f.write("=" * 100 + "\n")
        f.write("\n")

        # REGULAR - UPGMA
        f.write("-" * 100 + "\n")
        f.write(f"TREE 1/6: UPGMA (Regular)\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(upgma_tree, file=f)
        f.write("\n")
        f.write(f"Newick: {upgma_tree.to_newick()};\n")
        f.write("\n")

        # REGULAR - BioNJ
        f.write("-" * 100 + "\n")
        f.write(f"TREE 2/6: BioNJ (Regular)\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(bionj_tree, file=f)
        f.write("\n")
        f.write(f"Newick: {bionj_tree.to_newick()};\n")
        f.write("\n")

        # REGULAR - ML
        f.write("-" * 100 + "\n")
        f.write(f"TREE 3/6: ML (Regular) - Log-likelihood: {ml_logl:.2f}\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(ml_tree, file=f)
        f.write("\n")
        f.write(f"Newick: {ml_tree.to_newick()};\n")
        f.write("\n")

        f.write("\n")
        f.write("=" * 100 + "\n")
        f.write(f"MODE 2: DEDUPLICATED TREES ({len(dereplicated_seqs)} sequences - smart deduplication at 99.5%)\n")
        f.write("=" * 100 + "\n")
        f.write("\n")

        # DEDUPLICATED - UPGMA
        f.write("-" * 100 + "\n")
        f.write(f"TREE 4/6: UPGMA (Deduplicated)\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(upgma_tree_derep, file=f)
        f.write("\n")
        f.write(f"Newick: {upgma_tree_derep.to_newick()};\n")
        f.write("\n")

        # DEDUPLICATED - BioNJ
        f.write("-" * 100 + "\n")
        f.write(f"TREE 5/6: BioNJ (Deduplicated)\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(bionj_tree_derep, file=f)
        f.write("\n")
        f.write(f"Newick: {bionj_tree_derep.to_newick()};\n")
        f.write("\n")

        # DEDUPLICATED - ML
        f.write("-" * 100 + "\n")
        f.write(f"TREE 6/6: ML (Deduplicated) - Log-likelihood: {ml_logl_derep:.2f}\n")
        f.write("-" * 100 + "\n")
        f.write("\n")
        print_tree_ascii(ml_tree_derep, file=f)
        f.write("\n")
        f.write(f"Newick: {ml_tree_derep.to_newick()};\n")
        f.write("\n")

        f.write("\n")
        f.write("=" * 100 + "\n")
        f.write("END OF COMPARISON\n")
        f.write("=" * 100 + "\n")

    print(f"  Saved: {output_file}")
    print()

    # Generate individual files: ASCII, Newick, ETE3 PDF, ETE3 PNG for each tree
    print("Generating individual tree files (ASCII + Newick + ETE3)...")
    print("-" * 80)

    from rrna_phylo.visualization.ete3_viz import visualize_tree

    trees_to_save = [
        (upgma_tree, "regular_upgma", "Regular UPGMA", None),
        (bionj_tree, "regular_bionj", "Regular BioNJ", None),
        (ml_tree, "regular_ml", f"Regular ML", ml_logl),
        (upgma_tree_derep, "dedup_upgma", "Deduplicated UPGMA", None),
        (bionj_tree_derep, "dedup_bionj", "Deduplicated BioNJ", None),
        (ml_tree_derep, "dedup_ml", f"Deduplicated ML", ml_logl_derep),
    ]

    for tree, name, desc, logl in trees_to_save:
        # ASCII
        ascii_file = output_dir / f"{name}_ascii.txt"
        with open(ascii_file, 'w') as f:
            f.write(f"{desc}\n")
            f.write("-" * 70 + "\n")
            if logl is not None:
                f.write(f"Log-likelihood: {logl:.2f}\n")
                f.write("-" * 70 + "\n")
            print_tree_ascii(tree, file=f)

        # Newick
        newick_file = output_dir / f"{name}.nwk"
        with open(newick_file, 'w') as f:
            f.write(tree.to_newick() + ";")

        # ETE3 PDF
        pdf_file = output_dir / f"{name}.pdf"
        try:
            visualize_tree(str(newick_file), str(pdf_file), title=desc, dpi=300)
        except Exception as e:
            print(f"  WARNING: ETE3 PDF failed for {name}: {e}")

        # ETE3 PNG
        png_file = output_dir / f"{name}.png"
        try:
            visualize_tree(str(newick_file), str(png_file), title=desc, dpi=300)
        except Exception as e:
            print(f"  WARNING: ETE3 PNG failed for {name}: {e}")

        print(f"  {desc}: ASCII, Newick, PDF, PNG")

    print()

    # ========================================================================
    # SUMMARY
    # ========================================================================

    print("=" * 80)
    print("PIPELINE COMPLETE!")
    print("=" * 80)
    print()
    print(f"Output directory: {output_dir}/")
    print()
    print("Generated files:")
    print(f"  1. COMPLETE_TREE_COMPARISON.txt - All 6 trees combined")
    print(f"  2. aligned_sequences.fasta - Aligned sequences")
    print(f"  3-8.  6 ASCII files (individual tree visualizations)")
    print(f"  9-14. 6 Newick files")
    print(f"  15-20. 6 ETE3 PDF files")
    print(f"  21-26. 6 ETE3 PNG files")
    print()
    print("Tree Summary:")
    print(f"  Regular mode: {len(aligned_seqs)} sequences (3 trees: UPGMA, BioNJ, ML)")
    print(f"  Deduplicated mode: {len(dereplicated_seqs)} sequences (3 trees: UPGMA, BioNJ, ML)")
    print(f"  Total: 6 trees × 4 formats = 24 individual files + combined file + aligned FASTA")
    print()

    # Count files
    all_files = list(output_dir.rglob("*"))
    file_count = len([f for f in all_files if f.is_file()])
    print(f"Generated {file_count} files total")
    print()

    return {
        "output_file": output_file,
        "regular": {
            "sequences": len(aligned_seqs),
            "upgma": upgma_tree,
            "bionj": bionj_tree,
            "ml": ml_tree,
            "ml_logl": ml_logl
        },
        "deduplicated": {
            "sequences": len(dereplicated_seqs),
            "upgma": upgma_tree_derep,
            "bionj": bionj_tree_derep,
            "ml": ml_tree_derep,
            "ml_logl": ml_logl_derep,
            "stats": stats
        }
    }


if __name__ == "__main__":
    result = test_full_phylogenetic_pipeline()
    print("[OK] Full pipeline test completed successfully!")
    print(f"[OK] View results in: {result['output_file']}")

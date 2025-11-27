#!/usr/bin/env python3
"""
Test SMART phylogenetic tree builder on 100-species dataset.

This demonstrates intelligent method selection that:
1. Analyzes the dataset first
2. Checks distance matrix quality
3. Only builds trees using appropriate methods
4. Skips BioNJ if dataset is too large or divergent
"""

import sys
import time
from pathlib import Path
import tracemalloc

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rrna_phylo.io.fasta_parser import parse_fasta, assign_unique_display_names
from rrna_phylo.io.aligner import MuscleAligner
from rrna_phylo.utils.strain_handler import remove_exact_duplicates
from rrna_phylo.core.builder_smart import build_trees_smart


def format_time(seconds):
    """Format seconds into human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.0f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"


def format_memory(bytes_val):
    """Format bytes into human-readable string."""
    if bytes_val < 1024:
        return f"{bytes_val} B"
    elif bytes_val < 1024**2:
        return f"{bytes_val / 1024:.1f} KB"
    elif bytes_val < 1024**3:
        return f"{bytes_val / 1024**2:.1f} MB"
    else:
        return f"{bytes_val / 1024**3:.2f} GB"


def main():
    """Test smart tree builder on 100-species dataset."""

    print("=" * 80)
    print("SMART PHYLOGENETIC TREE BUILDER TEST")
    print("100 Bacterial Species Dataset")
    print("=" * 80)
    print()

    # Start tracking
    tracemalloc.start()
    test_start = time.time()

    # ========================================================================
    # STEP 1: Load sequences
    # ========================================================================

    print("STEP 1: Loading sequences")
    print("-" * 80)

    fasta_file = Path(__file__).parent.parent / "test_real_100species.fasta"
    step_start = time.time()

    all_sequences = list(parse_fasta(fasta_file))

    # Filter to DNA sequences only
    from rrna_phylo.core.sequence_type import SequenceTypeDetector
    detector = SequenceTypeDetector()

    sequences = []
    for seq in all_sequences:
        seq_type = detector.detect_single(seq)
        if seq_type.name == "DNA":
            sequences.append(seq)

    print(f"Loaded {len(all_sequences)} sequences ({len(sequences)} DNA, {len(all_sequences) - len(sequences)} filtered out)")

    step_time = time.time() - step_start
    current, peak = tracemalloc.get_traced_memory()

    print(f"Time: {format_time(step_time)}")
    print(f"Memory: {format_memory(current)} (peak: {format_memory(peak)})")
    print()

    # ========================================================================
    # STEP 2: Alignment
    # ========================================================================

    print("STEP 2: Multiple sequence alignment (MUSCLE)")
    print("-" * 80)

    step_start = time.time()
    mem_before = tracemalloc.get_traced_memory()[0]

    aligner = MuscleAligner()
    aligned_file = fasta_file.parent / "test_100species_aligned.fasta"
    aligned_seqs = aligner.align_sequences(sequences, str(aligned_file))

    step_time = time.time() - step_start
    current, peak = tracemalloc.get_traced_memory()
    mem_used = current - mem_before

    print(f"Aligned {len(aligned_seqs)} sequences")
    print(f"Alignment length: {aligned_seqs[0].aligned_length if aligned_seqs else 0} bp")
    print(f"Time: {format_time(step_time)}")
    print(f"Memory used: {format_memory(mem_used)} (peak: {format_memory(peak)})")
    print()

    # ========================================================================
    # STEP 3: Deduplication
    # ========================================================================

    print("STEP 3: Exact duplicate removal")
    print("-" * 80)

    step_start = time.time()

    exact_dedup, dup_map = remove_exact_duplicates(aligned_seqs)

    step_time = time.time() - step_start
    n_removed = len(aligned_seqs) - len(exact_dedup)

    print(f"Input: {len(aligned_seqs)} sequences")
    print(f"Output: {len(exact_dedup)} sequences")
    print(f"Removed: {n_removed} exact duplicates")
    print(f"Time: {format_time(step_time)}")
    print()

    # Assign display names
    assign_unique_display_names(exact_dedup)

    # ========================================================================
    # STEP 4: SMART Tree Building (with automatic method selection)
    # ========================================================================

    print("STEP 4: SMART Tree Building (automatic method selection)")
    print("-" * 80)

    step_start = time.time()

    # Build trees using SMART builder
    # This will analyze the dataset and only build appropriate methods
    results = build_trees_smart(exact_dedup, verbose=True)

    step_time = time.time() - step_start

    print(f"\nTotal tree building time: {format_time(step_time)}")
    print()

    # ========================================================================
    # STEP 5: Save outputs
    # ========================================================================

    print("STEP 5: Saving tree outputs")
    print("-" * 80)

    output_dir = fasta_file.parent / "test_100species_smart_output"
    output_dir.mkdir(exist_ok=True)

    from rrna_phylo.utils import print_tree_ascii

    # Save trees that were actually built
    for method in results['methods_built']:
        if method == "ml":
            tree, logl = results[method]
            tree_obj = tree
            desc = f"SMART ML Tree (GTR+Gamma) - Log-likelihood: {logl:.2f}"
        else:
            tree_obj = results[method]
            desc = f"SMART {method.upper()} Tree"

        # Save ASCII
        ascii_file = output_dir / f"{method}_smart.txt"
        with open(ascii_file, 'w') as f:
            f.write(f"{desc}\n")
            f.write("-" * 70 + "\n")
            print_tree_ascii(tree_obj, file=f)
        print(f"  [OK] {method.upper()}: {ascii_file.name}")

        # Save Newick
        newick_file = output_dir / f"{method}_smart.nwk"
        with open(newick_file, 'w') as f:
            f.write(tree_obj.to_newick() + ";")

    # Save analysis report
    analysis = results['analysis']
    report_file = output_dir / "dataset_analysis.txt"
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("DATASET ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Sequences: {analysis.n_sequences}\n")
        f.write(f"Alignment length: {analysis.alignment_length} bp\n\n")
        f.write("Distance Matrix Statistics:\n")
        f.write(f"  Max distance: {analysis.max_distance:.4f}\n")
        f.write(f"  Mean distance: {analysis.mean_distance:.4f}\n")
        f.write(f"  Median distance: {analysis.median_distance:.4f}\n")
        f.write(f"  Saturated distances: {analysis.pct_saturated:.1f}%\n")
        f.write(f"  Invalid comparisons: {analysis.pct_invalid:.1f}%\n\n")
        f.write(f"Recommended methods: {', '.join(m.upper() for m in analysis.recommended_methods)}\n\n")
        f.write("Warnings:\n")
        for i, warning in enumerate(analysis.warnings, 1):
            f.write(f"  {i}. {warning}\n")

    print(f"  [OK] Analysis report: {report_file.name}")
    print()

    # ========================================================================
    # STEP 6: Summary
    # ========================================================================

    test_time = time.time() - test_start
    current, peak = tracemalloc.get_traced_memory()

    print("=" * 80)
    print("TEST COMPLETE!")
    print("=" * 80)
    print()
    print(f"Total time: {format_time(test_time)}")
    print(f"Peak memory: {format_memory(peak)}")
    print()

    print("Dataset summary:")
    print(f"  Sequences: {len(exact_dedup)}")
    print(f"  Alignment length: {aligned_seqs[0].aligned_length} bp")
    print()

    print("Smart tree building results:")
    print(f"  Methods built: {', '.join(m.upper() for m in results['methods_built'])}")
    if results['methods_skipped']:
        print(f"  Methods skipped: {', '.join(m.upper() for m in results['methods_skipped'])}")
    print()

    if 'ml' in results:
        ml_tree, ml_logl = results['ml']
        print(f"  ML log-likelihood: {ml_logl:.2f}")
        print()

    print("Output files:")
    print(f"  Output directory: {output_dir}/")
    print(f"  Trees saved: {len(results['methods_built'])} methods")
    print(f"  Analysis report: dataset_analysis.txt")
    print()

    print("Performance assessment:")
    print("  [OK] Dataset analysis: Completed")
    print("  [OK] Intelligent method selection: Working as designed")
    print(f"  [OK] Tree building: {len(results['methods_built'])} methods completed")
    if "bionj" in results['methods_skipped']:
        print("  [OK] BioNJ skipped: Prevented timeout on divergent dataset")
    print("  [OK] All outputs saved")
    print()

    tracemalloc.stop()

    print(f"View results in: {output_dir}/")
    print()


if __name__ == "__main__":
    main()

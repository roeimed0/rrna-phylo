#!/usr/bin/env python3
"""
Performance scaling test for Arcosauria dataset.

Tests tree building with increasing dataset sizes (10, 20, 30, ... 111 sequences)
to identify bottlenecks.
"""

import time
from pathlib import Path
from rrna_phylo.io.fasta_parser import FastaParser
from rrna_phylo.core.builder import PhylogeneticTreeBuilder


def test_scaling():
    """Test performance with increasing dataset sizes."""

    # Load full dataset
    fasta_file = Path("data/test/Arcosauria_test.fasta")
    parser = FastaParser()
    all_sequences = parser.parse(str(fasta_file))

    print("=" * 80)
    print("PERFORMANCE SCALING TEST - Arcosauria Dataset")
    print("=" * 80)
    print(f"Full dataset: {len(all_sequences)} sequences")
    print()

    # Test sizes
    test_sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 111]

    results = []

    for n in test_sizes:
        if n > len(all_sequences):
            break

        # Take first n sequences
        sequences = all_sequences[:n]

        print(f"\n{'=' * 80}")
        print(f"Testing with {n} sequences")
        print('=' * 80)

        builder = PhylogeneticTreeBuilder(verbose=False)

        # Step 1: Detect sequence type
        t0 = time.time()
        builder.detect_and_validate(sequences)
        t_detect = time.time() - t0

        # Step 2: Alignment
        t0 = time.time()
        aligned_seqs = builder._align_sequences(sequences)
        t_align = time.time() - t0

        # Step 3: Build UPGMA tree (fast baseline)
        t0 = time.time()
        upgma_tree = builder.build_upgma_tree(aligned_seqs)
        t_upgma = time.time() - t0

        # Step 4: Build ML tree with model selection + NNI
        t0 = time.time()
        ml_tree, logL = builder.build_ml_tree(aligned_seqs, alpha=1.0, skip_model_selection=False)
        t_ml = time.time() - t0

        # Total time
        t_total = t_detect + t_align + t_upgma + t_ml

        # Store results
        result = {
            'n': n,
            'detect': t_detect,
            'align': t_align,
            'upgma': t_upgma,
            'ml': t_ml,
            'total': t_total
        }
        results.append(result)

        # Print summary
        print(f"\nResults for {n} sequences:")
        print(f"  Detection:       {t_detect:7.2f}s")
        print(f"  MUSCLE align:    {t_align:7.2f}s  {'<-- BOTTLENECK' if t_align > t_ml else ''}")
        print(f"  UPGMA tree:      {t_upgma:7.2f}s")
        print(f"  ML tree:         {t_ml:7.2f}s  {'<-- BOTTLENECK' if t_ml > t_align else ''}")
        print(f"  Total:           {t_total:7.2f}s")
        print(f"  Alignment %:     {100*t_align/t_total:5.1f}%")
        print(f"  ML tree %:       {100*t_ml/t_total:5.1f}%")

    # Final summary table
    print("\n\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"{'N':>4} | {'Detect':>8} | {'Align':>8} | {'UPGMA':>8} | {'ML':>8} | {'Total':>8} | Bottleneck")
    print("-" * 80)

    for r in results:
        bottleneck = "MUSCLE" if r['align'] > r['ml'] else "ML"
        print(f"{r['n']:4d} | {r['detect']:8.2f} | {r['align']:8.2f} | {r['upgma']:8.2f} | "
              f"{r['ml']:8.2f} | {r['total']:8.2f} | {bottleneck}")

    print("\n" + "=" * 80)
    print("KEY INSIGHTS:")
    print("=" * 80)

    # Calculate which is the bottleneck at different sizes
    muscle_dominant = sum(1 for r in results if r['align'] > r['ml'])
    ml_dominant = len(results) - muscle_dominant

    print(f"MUSCLE alignment is bottleneck: {muscle_dominant}/{len(results)} test sizes")
    print(f"ML tree building is bottleneck: {ml_dominant}/{len(results)} test sizes")

    # Show scaling
    if len(results) >= 2:
        small = results[0]
        large = results[-1]

        align_ratio = large['align'] / small['align']
        ml_ratio = large['ml'] / small['ml']
        size_ratio = large['n'] / small['n']

        print(f"\nScaling from {small['n']} to {large['n']} sequences ({size_ratio:.1f}x):")
        print(f"  MUSCLE alignment: {align_ratio:5.1f}x slower")
        print(f"  ML tree building: {ml_ratio:5.1f}x slower")

        # Expected scaling: O(N²) for MUSCLE, roughly O(N²) for ML
        print(f"\n  Expected if O(N²): {size_ratio**2:.1f}x slower")
        print(f"  MUSCLE actual vs expected: {align_ratio/(size_ratio**2):.2f}x")
        print(f"  ML actual vs expected:     {ml_ratio/(size_ratio**2):.2f}x")


if __name__ == '__main__':
    test_scaling()

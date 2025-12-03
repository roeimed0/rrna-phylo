#!/usr/bin/env python3
"""
Profile the ML Level 4 pipeline on 87 sequences to identify hotspots.

This will show exactly where time is spent:
- Model selection
- NNI search
- Branch length optimization
- Likelihood calculations
"""

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import cProfile
import pstats
import time
from io import StringIO

from rrna_phylo.io.fasta_parser import parse_fasta
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4

# Load 87 sequences
fasta_file = "test_100species_aligned.fasta"
all_sequences = parse_fasta(fasta_file)
sequences = all_sequences[:87]

print("=" * 80)
print("PROFILING ML LEVEL 4 - 87 SEQUENCES")
print("=" * 80)
print(f"Dataset: {len(sequences)} sequences")
print(f"Alignment length: {len(sequences[0].sequence)} bp")
print()

# Create profiler
profiler = cProfile.Profile()

# Run with profiling
print("Running ML Level 4 with profiling enabled...")
print()

start = time.time()
profiler.enable()

tree, logL, metadata = build_ml_tree_level4(
    sequences,
    model='auto',
    alpha=None,
    tree_search='nni',
    max_iterations=5,
    criterion='BIC',
    test_gamma=False,
    use_gpu='auto',
    verbose=True
)

profiler.disable()
elapsed = time.time() - start

print()
print("=" * 80)
print("PROFILING RESULTS")
print("=" * 80)
print(f"Total time: {elapsed:.2f}s")
print(f"Log-likelihood: {logL:.2f}")
print(f"Model: {metadata['selected_model']}")
print(f"NNI improvements: {metadata['n_nni_improvements']}")
print()

# Save detailed profile to file
with open('profile_87seqs_full.txt', 'w') as f:
    ps = pstats.Stats(profiler, stream=f)
    ps.sort_stats('cumulative')
    f.write("=" * 80 + "\n")
    f.write("FULL PROFILE (sorted by cumulative time)\n")
    f.write("=" * 80 + "\n\n")
    ps.print_stats()

    f.write("\n\n" + "=" * 80 + "\n")
    f.write("PROFILE (sorted by total time)\n")
    f.write("=" * 80 + "\n\n")
    ps.sort_stats('tottime')
    ps.print_stats()

print("Full profile saved to: profile_87seqs_full.txt")
print()

# Print top 40 functions by cumulative time
print("-" * 80)
print("TOP 40 FUNCTIONS BY CUMULATIVE TIME")
print("-" * 80)
ps = pstats.Stats(profiler)
ps.sort_stats('cumulative')
ps.print_stats(40)

print()
print("-" * 80)
print("TOP 40 FUNCTIONS BY TOTAL TIME")
print("-" * 80)
ps.sort_stats('tottime')
ps.print_stats(40)

print()
print("=" * 80)
print("PROFILING COMPLETE")
print("=" * 80)
print()
print("Analysis:")
print("  - Check 'ncalls' column for functions called many times")
print("  - Check 'tottime' for functions that are slow themselves")
print("  - Check 'cumtime' for functions with slow children")
print()
print("Look for:")
print("  1. compute_likelihood / _conditional_likelihood - core bottleneck")
print("  2. perform_nni_swap / undo_nni_swap - topology changes")
print("  3. optimize_branch_lengths - branch optimization")
print("  4. Model selection functions - initial model comparison")
print()
print("Expected hotspots:")
print("  - Model selection: ~10s (5 models × 87 sequences)")
print("  - NNI search: ~2s (with GPU)")
print("  - If NNI shows many likelihood calls, consider two-stage evaluation")
print()

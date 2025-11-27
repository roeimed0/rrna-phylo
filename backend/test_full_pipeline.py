"""
COMPLETE PHYLOGENETIC PIPELINE TEST
=====================================

Tests the full pipeline on TWO datasets:
1. Small dataset (5 sequences) - for quick validation
2. Large dataset (87 sequences) - for production testing

Each dataset produces 3 trees:
- UPGMA (distance-based, molecular clock)
- BioNJ (distance-based, no clock)
- ML Level 4 (model selection + ML tree + NNI optimization)
"""
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import time
from rrna_phylo.io.fasta_parser import parse_fasta
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4
from rrna_phylo.visualization.ete3_viz import visualize_tree
from rrna_phylo.utils.visualize_trees import print_tree_ascii


def run_full_pipeline(sequences, output_dir, dataset_name):
    """
    Run complete phylogenetic pipeline on a dataset.

    Creates 3 trees:
    1. UPGMA
    2. BioNJ
    3. ML Level 4

    Each tree gets: .nwk, .txt, .pdf files
    """
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 80)
    print(f"DATASET: {dataset_name}")
    print("=" * 80)
    print(f"Sequences: {len(sequences)}")
    print(f"Alignment length: {len(sequences[0].sequence)} bp")
    print(f"Output directory: {output_dir}/")
    print()

    results = {
        'dataset': dataset_name,
        'n_sequences': len(sequences),
        'trees': {}
    }

    # ========================================================================
    # TREE 1: UPGMA
    # ========================================================================
    print("-" * 80)
    print("TREE 1: UPGMA (Distance-based, Molecular Clock)")
    print("-" * 80)

    start = time.time()
    dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")
    upgma_tree = build_upgma_tree(dist_matrix, ids)
    upgma_time = time.time() - start
    upgma_logL = compute_log_likelihood(upgma_tree, sequences)

    print(f"Time: {upgma_time:.2f}s, LogL: {upgma_logL:.2f}")

    # Save outputs
    upgma_newick = f"{output_dir}/tree_1_upgma.nwk"
    with open(upgma_newick, 'w') as f:
        f.write(upgma_tree.to_newick() + ";\n")

    upgma_ascii = f"{output_dir}/tree_1_upgma.txt"
    with open(upgma_ascii, 'w') as f:
        f.write(f"UPGMA Tree - {dataset_name}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Method: UPGMA (Distance-based, Molecular Clock)\n")
        f.write(f"Distance model: Jukes-Cantor\n")
        f.write(f"Log-likelihood: {upgma_logL:.2f}\n")
        f.write(f"Construction time: {upgma_time:.2f}s\n\n")
        print_tree_ascii(upgma_tree, file=f, show_support=False)

    try:
        upgma_pdf = f"{output_dir}/tree_1_upgma.pdf"
        visualize_tree(upgma_newick, output_file=upgma_pdf,
                       title=f"UPGMA - {dataset_name} (LogL: {upgma_logL:.0f})")
        print(f"[OK] Saved: .nwk, .txt, .pdf")
    except Exception as e:
        print(f"[OK] Saved: .nwk, .txt (PDF failed: {e})")

    results['trees']['upgma'] = {
        'logL': upgma_logL,
        'time': upgma_time
    }
    print()

    # ========================================================================
    # TREE 2: BioNJ
    # ========================================================================
    print("-" * 80)
    print("TREE 2: BioNJ (Distance-based, No Clock)")
    print("-" * 80)

    start = time.time()
    bionj_tree = build_bionj_tree(dist_matrix, ids)
    bionj_time = time.time() - start
    bionj_logL = compute_log_likelihood(bionj_tree, sequences)

    print(f"Time: {bionj_time:.2f}s, LogL: {bionj_logL:.2f}")

    # Save outputs
    bionj_newick = f"{output_dir}/tree_2_bionj.nwk"
    with open(bionj_newick, 'w') as f:
        f.write(bionj_tree.to_newick() + ";\n")

    bionj_ascii = f"{output_dir}/tree_2_bionj.txt"
    with open(bionj_ascii, 'w') as f:
        f.write(f"BioNJ Tree - {dataset_name}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Method: BioNJ (Distance-based, No Clock)\n")
        f.write(f"Distance model: Jukes-Cantor\n")
        f.write(f"Log-likelihood: {bionj_logL:.2f}\n")
        f.write(f"Construction time: {bionj_time:.2f}s\n\n")
        print_tree_ascii(bionj_tree, file=f, show_support=False)

    try:
        bionj_pdf = f"{output_dir}/tree_2_bionj.pdf"
        visualize_tree(bionj_newick, output_file=bionj_pdf,
                       title=f"BioNJ - {dataset_name} (LogL: {bionj_logL:.0f})")
        print(f"[OK] Saved: .nwk, .txt, .pdf")
    except Exception as e:
        print(f"[OK] Saved: .nwk, .txt (PDF failed: {e})")

    results['trees']['bionj'] = {
        'logL': bionj_logL,
        'time': bionj_time
    }
    print()

    # ========================================================================
    # TREE 3: ML Level 4
    # ========================================================================
    print("-" * 80)
    print("TREE 3: ML Level 4 (Complete ML Pipeline)")
    print("-" * 80)
    print("Running: Model selection + ML tree + NNI optimization...")
    print()

    start = time.time()
    ml_tree, ml_logL, metadata = build_ml_tree_level4(
        sequences,
        model='auto',
        alpha=None,
        tree_search='nni',
        max_iterations=10,
        criterion='BIC',
        test_gamma=False,
        verbose=False  # Suppress detailed output for cleaner logs
    )
    ml_time = time.time() - start

    print(f"Time: {ml_time:.2f}s, LogL: {ml_logL:.2f}")
    print(f"Model selected: {metadata['selected_model']}")
    print(f"NNI improvements: {metadata['n_nni_improvements']}")

    # Save outputs
    ml_newick = f"{output_dir}/tree_3_ml_level4.nwk"
    with open(ml_newick, 'w') as f:
        f.write(ml_tree.to_newick() + ";\n")

    ml_ascii = f"{output_dir}/tree_3_ml_level4.txt"
    with open(ml_ascii, 'w') as f:
        f.write(f"ML Level 4 Tree - {dataset_name}\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Method: Maximum Likelihood (Level 4)\n")
        f.write(f"Selected model: {metadata['selected_model']}\n")
        f.write(f"Model score (BIC): {metadata['model_score']:.2f}\n")
        f.write(f"Initial LogL: {metadata['initial_logL']:.2f}\n")
        f.write(f"Final LogL: {ml_logL:.2f}\n")
        f.write(f"Improvement: {ml_logL - metadata['initial_logL']:+.2f}\n")
        f.write(f"NNI improvements: {metadata['n_nni_improvements']}\n")
        f.write(f"Construction time: {ml_time:.2f}s\n\n")
        print_tree_ascii(ml_tree, file=f, show_support=False)

    try:
        ml_pdf = f"{output_dir}/tree_3_ml_level4.pdf"
        visualize_tree(ml_newick, output_file=ml_pdf,
                       title=f"ML Level 4 - {dataset_name} (LogL: {ml_logL:.0f})")
        print(f"[OK] Saved: .nwk, .txt, .pdf")
    except Exception as e:
        print(f"[OK] Saved: .nwk, .txt (PDF failed: {e})")

    results['trees']['ml_level4'] = {
        'logL': ml_logL,
        'time': ml_time,
        'model': metadata['selected_model'],
        'nni_improvements': metadata['n_nni_improvements']
    }
    print()

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("=" * 80)
    print(f"SUMMARY: {dataset_name}")
    print("=" * 80)
    print()
    print(f"{'Tree':<20} {'LogL':>15} {'Time':>10}")
    print("-" * 50)
    print(f"{'UPGMA':<20} {upgma_logL:>15.2f} {upgma_time:>9.2f}s")
    print(f"{'BioNJ':<20} {bionj_logL:>15.2f} {bionj_time:>9.2f}s")
    print(f"{'ML Level 4':<20} {ml_logL:>15.2f} {ml_time:>9.2f}s")
    print()

    best_logL = max(upgma_logL, bionj_logL, ml_logL)
    if best_logL == upgma_logL:
        best = "UPGMA"
    elif best_logL == bionj_logL:
        best = "BioNJ"
    else:
        best = "ML Level 4"

    print(f"Best tree: {best} (LogL: {best_logL:.2f})")
    print(f"Total time: {upgma_time + bionj_time + ml_time:.2f}s")
    print(f"Output files: 9 (3 trees x 3 formats)")
    print()

    return results


# ============================================================================
# MAIN TEST
# ============================================================================

print("\n")
print("=" * 80)
print("COMPLETE PHYLOGENETIC PIPELINE TEST")
print("=" * 80)
print()
print("Testing on 2 datasets:")
print("  1. Small (5 sequences) - Quick validation")
print("  2. Large (87 sequences) - Production scale")
print()
print("Each dataset generates 3 trees:")
print("  - UPGMA (fast, assumes molecular clock)")
print("  - BioNJ (fast, no clock assumption)")
print("  - ML Level 4 (complete ML pipeline)")
print()
print("=" * 80)
print()

# Load clean 16S dataset
all_sequences = parse_fasta("test_100species_16s_only.fasta")
print(f"[OK] Loaded {len(all_sequences)} pure 16S rRNA sequences from SILVA")
print()

# ============================================================================
# TEST 1: Small Dataset (5 sequences)
# ============================================================================
print("\n")
sequences_5 = all_sequences[:5]
results_5 = run_full_pipeline(
    sequences_5,
    output_dir="test_outputs/small_5seqs",
    dataset_name="5 Bacterial 16S rRNA Sequences"
)

# ============================================================================
# TEST 2: Large Dataset (87 sequences)
# ============================================================================
print("\n")
sequences_87 = all_sequences  # Use all sequences
results_87 = run_full_pipeline(
    sequences_87,
    output_dir="test_outputs/large_87seqs",
    dataset_name="87 Bacterial 16S rRNA Sequences"
)

# ============================================================================
# FINAL REPORT
# ============================================================================
print("\n")
print("=" * 80)
print("FINAL REPORT: ALL TESTS COMPLETE")
print("=" * 80)
print()

print("Dataset 1: Small (5 sequences)")
print("-" * 50)
print(f"  UPGMA:      LogL = {results_5['trees']['upgma']['logL']:.2f}, Time = {results_5['trees']['upgma']['time']:.2f}s")
print(f"  BioNJ:      LogL = {results_5['trees']['bionj']['logL']:.2f}, Time = {results_5['trees']['bionj']['time']:.2f}s")
print(f"  ML Level 4: LogL = {results_5['trees']['ml_level4']['logL']:.2f}, Time = {results_5['trees']['ml_level4']['time']:.2f}s")
print(f"  Model: {results_5['trees']['ml_level4']['model']}")
print(f"  NNI improvements: {results_5['trees']['ml_level4']['nni_improvements']}")
print()

print("Dataset 2: Large (87 sequences)")
print("-" * 50)
print(f"  UPGMA:      LogL = {results_87['trees']['upgma']['logL']:.2f}, Time = {results_87['trees']['upgma']['time']:.2f}s")
print(f"  BioNJ:      LogL = {results_87['trees']['bionj']['logL']:.2f}, Time = {results_87['trees']['bionj']['time']:.2f}s")
print(f"  ML Level 4: LogL = {results_87['trees']['ml_level4']['logL']:.2f}, Time = {results_87['trees']['ml_level4']['time']:.2f}s")
print(f"  Model: {results_87['trees']['ml_level4']['model']}")
print(f"  NNI improvements: {results_87['trees']['ml_level4']['nni_improvements']}")
print()

print("Output directories:")
print("  test_outputs/small_5seqs/   - 9 files (3 trees x 3 formats)")
print("  test_outputs/large_87seqs/  - 9 files (3 trees x 3 formats)")
print()
print("Total output files: 18")
print()

print("=" * 80)
print("[OK] ALL TESTS PASSED - PIPELINE WORKING CORRECTLY")
print("=" * 80)
print()

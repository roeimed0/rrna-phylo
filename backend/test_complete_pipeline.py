"""
Complete Pipeline Test: FASTA → Alignment → 3 Trees → Visualization

Tests the entire workflow:
1. Load FASTA file
2. Auto-align with MUSCLE
3. Build 3 tree types (UPGMA, BioNJ, ML)
4. Generate ASCII and Newick outputs
5. Test with and without bootstrap
"""

import subprocess
import sys
from pathlib import Path

def test_pipeline(test_name, args, check_files):
    """Run a pipeline test and verify outputs."""
    print("\n" + "=" * 80)
    print(f"TEST: {test_name}")
    print("=" * 80)
    print(f"Command: rrna-phylo {' '.join(args)}")
    print("-" * 80)

    cmd = [sys.executable, "-m", "rrna_phylo.cli"] + args
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)

    # Check for expected output files
    success = True
    for file_path in check_files:
        if Path(file_path).exists():
            size = Path(file_path).stat().st_size
            print(f"  [OK] {file_path} ({size} bytes)")
        else:
            print(f"  [FAIL] {file_path} NOT FOUND")
            success = False

    if result.returncode != 0:
        print(f"  [FAIL] Exit code: {result.returncode}")
        success = False

    return success

def main():
    print("""
================================================================================
                     COMPLETE PIPELINE TESTS

  Testing: FASTA -> Alignment -> Tree Building -> Visualization
  Methods: UPGMA, BioNJ, ML
  Outputs: ASCII + Newick
================================================================================
""")

    test_file = "test_real_rrana.fasta"
    if not Path(test_file).exists():
        print(f"ERROR: Test file {test_file} not found!")
        return 1

    results = {}

    # Test 1: All methods, no bootstrap, ASCII only
    results['all_methods_ascii'] = test_pipeline(
        "All Methods - ASCII Only",
        [
            test_file,
            "--method", "all",
            "--output-format", "ascii",
            "--ignore-bias-warning",
            "-o", "test_pipeline_1/"
        ],
        [
            "test_pipeline_1/tree_upgma_ascii.txt",
            "test_pipeline_1/tree_bionj_ascii.txt",
            "test_pipeline_1/tree_ml_ascii.txt"
        ]
    )

    # Test 2: All methods, no bootstrap, Newick only
    results['all_methods_newick'] = test_pipeline(
        "All Methods - Newick Only",
        [
            test_file,
            "--method", "all",
            "--output-format", "newick",
            "--ignore-bias-warning",
            "-o", "test_pipeline_2/"
        ],
        [
            "test_pipeline_2/tree_upgma.nwk",
            "test_pipeline_2/tree_bionj.nwk",
            "test_pipeline_2/tree_ml.nwk"
        ]
    )

    # Test 3: All methods, no bootstrap, both formats
    results['all_methods_both'] = test_pipeline(
        "All Methods - Both Formats",
        [
            test_file,
            "--method", "all",
            "--output-format", "both",
            "--ignore-bias-warning",
            "-o", "test_pipeline_3/"
        ],
        [
            "test_pipeline_3/tree_upgma_ascii.txt",
            "test_pipeline_3/tree_upgma.nwk",
            "test_pipeline_3/tree_bionj_ascii.txt",
            "test_pipeline_3/tree_bionj.nwk",
            "test_pipeline_3/tree_ml_ascii.txt",
            "test_pipeline_3/tree_ml.nwk"
        ]
    )

    # Test 4: ML with bootstrap
    results['ml_bootstrap'] = test_pipeline(
        "ML with Bootstrap (10 replicates)",
        [
            test_file,
            "--method", "ml",
            "--bootstrap", "10",
            "--output-format", "both",
            "--ignore-bias-warning",
            "-o", "test_pipeline_4/"
        ],
        [
            "test_pipeline_4/tree_ml_ascii.txt",
            "test_pipeline_4/tree_ml.nwk"
        ]
    )

    # Test 5: BioNJ with bootstrap
    results['bionj_bootstrap'] = test_pipeline(
        "BioNJ with Bootstrap (10 replicates)",
        [
            test_file,
            "--method", "bionj",
            "--bootstrap", "10",
            "--output-format", "both",
            "--ignore-bias-warning",
            "-o", "test_pipeline_5/"
        ],
        [
            "test_pipeline_5/tree_bionj_ascii.txt",
            "test_pipeline_5/tree_bionj.nwk"
        ]
    )

    # Test 6: UPGMA with bootstrap
    results['upgma_bootstrap'] = test_pipeline(
        "UPGMA with Bootstrap (10 replicates)",
        [
            test_file,
            "--method", "upgma",
            "--bootstrap", "10",
            "--output-format", "both",
            "--ignore-bias-warning",
            "-o", "test_pipeline_6/"
        ],
        [
            "test_pipeline_6/tree_upgma_ascii.txt",
            "test_pipeline_6/tree_upgma.nwk"
        ]
    )

    # Test 7: Full workflow with all quality features
    results['full_workflow'] = test_pipeline(
        "Full Workflow - Dereplicate + Outgroup + Bootstrap",
        [
            test_file,
            "--dereplicate",
            "--outgroup", "AE004091*",
            "--method", "ml",
            "--bootstrap", "10",
            "--output-format", "both",
            "-o", "test_pipeline_7/"
        ],
        [
            "test_pipeline_7/tree_ml_ascii.txt",
            "test_pipeline_7/tree_ml.nwk",
            "test_pipeline_7/aligned_test_real_rrana.fasta"
        ]
    )

    # Test 8: Stratified sampling workflow
    results['stratified_workflow'] = test_pipeline(
        "Stratified Sampling Workflow",
        [
            test_file,
            "--stratify",
            "--max-per-species", "3",
            "--method", "all",
            "--output-format", "both",
            "-o", "test_pipeline_8/"
        ],
        [
            "test_pipeline_8/tree_upgma_ascii.txt",
            "test_pipeline_8/tree_bionj_ascii.txt",
            "test_pipeline_8/tree_ml_ascii.txt",
            "test_pipeline_8/tree_upgma.nwk",
            "test_pipeline_8/tree_bionj.nwk",
            "test_pipeline_8/tree_ml.nwk"
        ]
    )

    # Print summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    for test_name, result in results.items():
        status = "[PASS]" if result else "[FAIL]"
        print(f"{status} {test_name}")

    print("\n" + "=" * 80)
    print(f"RESULTS: {passed}/{total} tests passed")
    print("=" * 80)

    if passed == total:
        print("\n[OK] ALL TESTS PASSED! Core pipeline is complete.")
        return 0
    else:
        print(f"\n[FAIL] {total - passed} test(s) failed. Review errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())

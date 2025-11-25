"""
Test script demonstrating outgroup and database bias handling.

This shows all three features:
1. Stratified sampling (--stratify)
2. Dereplication (--dereplicate)
3. Outgroup handling (--outgroup, --suggest-outgroup)
"""

import subprocess
import sys
from pathlib import Path

def run_cli(args, description):
    """Run CLI command and show output."""
    print("\n" + "=" * 80)
    print(f"TEST: {description}")
    print("=" * 80)
    print(f"Command: rrna-phylo {' '.join(args)}")
    print("-" * 80)

    cmd = [sys.executable, "-m", "rrna_phylo.cli"] + args
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr, file=sys.stderr)

    return result.returncode == 0

def main():
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║              DATABASE BIAS AND OUTGROUP HANDLING TESTS                       ║
║                                                                              ║
║  This demonstrates the three key features for handling real-world data:     ║
║  1. Bias detection and stratified sampling                                  ║
║  2. Strain dereplication (multiple rRNA copies)                             ║
║  3. Outgroup selection for tree rooting                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

    test_file = Path("test_real_rrana.fasta")

    if not test_file.exists():
        print(f"Error: Test file {test_file} not found!")
        print("Please create it first with the 24 bacterial sequences")
        return 1

    # Test 1: Check for bias
    print("\n" + "#" * 80)
    print("# TEST 1: Detect sampling bias")
    print("#" * 80)
    run_cli([
        str(test_file),
        "--check-bias"
    ], "Analyze dataset for sampling bias")

    # Test 2: Suggest outgroup
    print("\n" + "#" * 80)
    print("# TEST 2: Suggest appropriate outgroups")
    print("#" * 80)
    run_cli([
        str(test_file),
        "--suggest-outgroup"
    ], "Get outgroup recommendations")

    # Test 3: Build tree with dereplication only
    print("\n" + "#" * 80)
    print("# TEST 3: Dereplication (24 sequences → 5 representatives)")
    print("#" * 80)
    run_cli([
        str(test_file),
        "--dereplicate",
        "-o", "test_dereplicated_only/",
        "--method", "ml",
        "--output-format", "ascii"
    ], "Remove multiple rRNA copies from same genome")

    # Test 4: Build tree with stratified sampling
    # (Won't do much for this dataset since it's already balanced)
    print("\n" + "#" * 80)
    print("# TEST 4: Stratified sampling (cap at 3 per species)")
    print("#" * 80)
    run_cli([
        str(test_file),
        "--stratify",
        "--max-per-species", "3",
        "-o", "test_stratified/",
        "--method", "ml",
        "--output-format", "ascii"
    ], "Balance species representation")

    # Test 5: Combined workflow (RECOMMENDED)
    print("\n" + "#" * 80)
    print("# TEST 5: RECOMMENDED WORKFLOW (dereplicate + outgroup)")
    print("#" * 80)
    run_cli([
        str(test_file),
        "--dereplicate",
        "--outgroup", "AE004091*",  # Pseudomonas aeruginosa
        "-o", "test_full_pipeline/",
        "--method", "ml",
        "--bootstrap", "10",
        "--output-format", "both"
    ], "Full pipeline: dereplication + outgroup rooting + bootstrap")

    print("\n" + "=" * 80)
    print("ALL TESTS COMPLETE!")
    print("=" * 80)
    print("""
Key Takeaways:
1. --check-bias: Always run this first to understand your data
2. --suggest-outgroup: Get recommendations for tree rooting
3. --dereplicate: Essential for genomes with multiple rRNA copies
4. --stratify: Use when you have overrepresented species (e.g., 1000 E. coli)
5. --outgroup: Specify outgroup pattern for rooting

Recommended workflow for publication:
  rrna-phylo sequences.fasta \\
    --dereplicate \\
    --stratify --max-per-species 10 \\
    --outgroup "Pseudomonas*" \\
    --bootstrap 100 \\
    -o results/
""")

    return 0

if __name__ == "__main__":
    sys.exit(main())

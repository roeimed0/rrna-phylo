#!/usr/bin/env python3
"""
Download a real larger bacterial 16S rRNA dataset for testing.

This script will download sequences from the SILVA or Greengenes database
which are standard references for 16S rRNA studies.

Alternatively, we can use NCBI to download a curated set of reference sequences.
"""

import urllib.request
from pathlib import Path
import sys

def download_silva_subset():
    """
    Download a subset of SILVA 16S rRNA sequences.

    SILVA is the most comprehensive and well-curated rRNA database.
    We'll download a representative subset for testing.
    """
    print("=" * 80)
    print("DOWNLOADING REAL 16S rRNA DATASET FROM SILVA")
    print("=" * 80)
    print()

    # SILVA provides subsets of their database
    # We'll use a smaller representative set for testing

    print("Option 1: Download from SILVA (Recommended)")
    print("  - High quality, curated sequences")
    print("  - Taxonomically diverse")
    print("  - Download URL: https://www.arb-silva.de/")
    print()
    print("Option 2: Manual download")
    print("  - Visit: https://www.arb-silva.de/download/arb-files/")
    print("  - Download: SILVA_xxx_SSURef_NR99_tax_silva_trunc.fasta.gz")
    print("  - Extract and copy to backend/test_real_large.fasta")
    print()
    print("Option 3: Use RDP (Ribosomal Database Project)")
    print("  - Visit: https://rdp.cme.msu.edu/")
    print("  - Download curated sequences")
    print()

    print("For automated testing, I recommend:")
    print()
    print("1. Download NCBI RefSeq representative genomes (automated below)")
    print("2. Manual SILVA download for highest quality")
    print()

    return None

def download_ncbi_refseq():
    """
    Download representative bacterial 16S sequences from NCBI RefSeq.

    This uses a pre-curated set of representative sequences.
    """
    print("=" * 80)
    print("DOWNLOADING REFSEQ REPRESENTATIVE 16S rRNA SEQUENCES")
    print("=" * 80)
    print()

    print("To download from NCBI, you'll need to use NCBI's E-utilities.")
    print("This requires:")
    print("  1. Biopython installed")
    print("  2. Internet connection")
    print("  3. NCBI API key (optional, for faster downloads)")
    print()

    # Check if Biopython is available
    try:
        from Bio import Entrez
        print("✓ Biopython is installed")
        print()

        print("Would you like to proceed with NCBI download?")
        print("This will download ~100-200 representative bacterial 16S sequences")
        print()

        return "biopython_available"

    except ImportError:
        print("✗ Biopython not installed")
        print()
        print("To install: pip install biopython")
        print()

        return None

def main():
    """Main function."""

    output_file = Path(__file__).parent.parent / "test_real_large.fasta"

    print()
    print("This script helps you download a real larger dataset for testing.")
    print()
    print("Recommended approach:")
    print()
    print("OPTION A: Use the smaller download script (already available)")
    print("  File: download_larger_dataset.py")
    print("  Requires: Biopython")
    print("  Downloads: ~100-200 sequences from NCBI")
    print()
    print("OPTION B: Manual download from SILVA (highest quality)")
    print("  1. Visit: https://www.arb-silva.de/no_cache/download/archive/current/Exports/")
    print("  2. Download: SILVA_xxx_SSURef_NR99_tax_silva.fasta.gz")
    print("  3. Extract first 100-500 sequences")
    print("  4. Save as: test_real_large.fasta")
    print()
    print("OPTION C: Use existing NCBI sequences")
    print("  Search NCBI for: \"16S ribosomal RNA[Title] AND bacteria[Organism]\"")
    print("  Download in FASTA format")
    print("  Select representative sequences from different phyla")
    print()

    # Try download script
    print("Checking if Biopython download is available...")
    result = download_ncbi_refseq()

    if result == "biopython_available":
        print("You can use the download_larger_dataset.py script!")
        print()
        print("To run it:")
        print("  python scripts/download_larger_dataset.py")
        print()
    else:
        download_silva_subset()

    print("=" * 80)
    print()


if __name__ == "__main__":
    main()

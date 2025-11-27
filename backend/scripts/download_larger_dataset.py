#!/usr/bin/env python3
"""
Download a larger bacterial 16S rRNA dataset for testing phylogenetic pipeline at scale.

This script creates a diverse test dataset with:
- 100-200 sequences from different bacterial phyla
- Mix of well-characterized reference strains
- Includes some closely related strains to test deduplication
- Real-world complexity (varying lengths, quality)

Dataset composition:
- Proteobacteria (40-50 sequences): E. coli, Salmonella, Pseudomonas, etc.
- Firmicutes (20-30 sequences): Bacillus, Staphylococcus, Streptococcus, etc.
- Actinobacteria (20-30 sequences): Mycobacterium, Corynebacterium, etc.
- Bacteroidetes (10-20 sequences): Bacteroides, Prevotella, etc.
- Other phyla (10-20 sequences): Cyanobacteria, Spirochaetes, etc.
"""

from Bio import Entrez, SeqIO
import time
from pathlib import Path

# NCBI requires email for Entrez
Entrez.email = "your.email@example.com"  # Replace with your email

def download_sequences(query, max_sequences, delay=0.5):
    """Download sequences from NCBI using Entrez."""
    print(f"Searching: {query}")
    print(f"Max sequences: {max_sequences}")

    try:
        # Search for sequences
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_sequences, idtype="acc")
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"Found {len(id_list)} sequences")

        if not id_list:
            return []

        # Fetch sequences
        time.sleep(delay)  # Be nice to NCBI servers
        handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        return sequences

    except Exception as e:
        print(f"Error: {e}")
        return []

def main():
    """Create larger test dataset."""

    output_file = Path(__file__).parent.parent / "test_larger_dataset.fasta"

    # Define queries for diverse bacterial groups
    queries = [
        # Proteobacteria (most diverse, well-studied)
        ("Escherichia coli[Organism] AND 16S ribosomal RNA[Title]", 15),
        ("Salmonella enterica[Organism] AND 16S ribosomal RNA[Title]", 10),
        ("Pseudomonas aeruginosa[Organism] AND 16S ribosomal RNA[Title]", 10),
        ("Klebsiella pneumoniae[Organism] AND 16S ribosomal RNA[Title]", 8),
        ("Vibrio cholerae[Organism] AND 16S ribosomal RNA[Title]", 5),
        ("Helicobacter pylori[Organism] AND 16S ribosomal RNA[Title]", 5),

        # Firmicutes (important pathogens and commensals)
        ("Bacillus subtilis[Organism] AND 16S ribosomal RNA[Title]", 10),
        ("Staphylococcus aureus[Organism] AND 16S ribosomal RNA[Title]", 8),
        ("Streptococcus pneumoniae[Organism] AND 16S ribosomal RNA[Title]", 8),
        ("Clostridium difficile[Organism] AND 16S ribosomal RNA[Title]", 5),
        ("Listeria monocytogenes[Organism] AND 16S ribosomal RNA[Title]", 5),

        # Actinobacteria (high GC content)
        ("Mycobacterium tuberculosis[Organism] AND 16S ribosomal RNA[Title]", 10),
        ("Corynebacterium diphtheriae[Organism] AND 16S ribosomal RNA[Title]", 5),
        ("Streptomyces coelicolor[Organism] AND 16S ribosomal RNA[Title]", 5),

        # Bacteroidetes (anaerobes, gut microbiota)
        ("Bacteroides fragilis[Organism] AND 16S ribosomal RNA[Title]", 8),
        ("Prevotella[Organism] AND 16S ribosomal RNA[Title]", 5),

        # Other important phyla
        ("Chlamydia trachomatis[Organism] AND 16S ribosomal RNA[Title]", 5),
        ("Cyanobacteria[Organism] AND 16S ribosomal RNA[Title]", 5),
        ("Treponema pallidum[Organism] AND 16S ribosomal RNA[Title]", 5),
    ]

    all_sequences = []

    print("=" * 80)
    print("DOWNLOADING LARGER TEST DATASET")
    print("=" * 80)
    print()

    for query, max_seqs in queries:
        seqs = download_sequences(query, max_seqs)
        all_sequences.extend(seqs)
        print(f"Total collected: {len(all_sequences)}")
        print()
        time.sleep(1)  # Be respectful to NCBI servers

    # Write to file
    print(f"Writing {len(all_sequences)} sequences to {output_file}")
    SeqIO.write(all_sequences, output_file, "fasta")

    print()
    print("=" * 80)
    print("DOWNLOAD COMPLETE!")
    print("=" * 80)
    print()
    print(f"Total sequences: {len(all_sequences)}")
    print(f"Output file: {output_file}")
    print()

    # Summary statistics
    lengths = [len(seq.seq) for seq in all_sequences]
    print("Sequence length statistics:")
    print(f"  Min: {min(lengths)} bp")
    print(f"  Max: {max(lengths)} bp")
    print(f"  Mean: {sum(lengths) / len(lengths):.0f} bp")
    print()

if __name__ == "__main__":
    main()

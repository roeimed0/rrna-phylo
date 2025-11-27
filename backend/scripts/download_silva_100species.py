#!/usr/bin/env python3
"""
Download 16S rRNA sequences from 100 diverse bacterial species.

This script downloads high-quality 16S rRNA sequences representing
100 different bacterial species across all major phyla.

Strategy:
- Use NCBI to download representative sequences
- Select type strains when possible
- Ensure phylogenetic diversity across:
  * Proteobacteria (Alpha, Beta, Gamma, Delta, Epsilon)
  * Firmicutes (Bacilli, Clostridia)
  * Actinobacteria
  * Bacteroidetes
  * Cyanobacteria
  * Other phyla

Target: ~100-150 sequences for robust phylogenetic testing
"""

from Bio import Entrez, SeqIO
import time
from pathlib import Path
from collections import defaultdict

# IMPORTANT: Set your email for NCBI
Entrez.email = "your.email@example.com"

def download_sequences_for_taxon(taxon_query, max_sequences, delay=0.5):
    """
    Download sequences for a specific taxonomic group.

    Args:
        taxon_query: NCBI search query
        max_sequences: Maximum number to download
        delay: Delay between requests (seconds)

    Returns:
        List of SeqRecord objects
    """
    print(f"  Searching: {taxon_query}")

    try:
        # Search for sequences
        handle = Entrez.esearch(
            db="nucleotide",
            term=taxon_query,
            retmax=max_sequences,
            idtype="acc"
        )
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"    Found {len(id_list)} sequences")

        if not id_list:
            return []

        # Fetch sequences in batches
        sequences = []
        batch_size = 20  # NCBI recommends smaller batches

        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i+batch_size]

            time.sleep(delay)  # Be nice to NCBI servers

            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=batch_ids,
                    rettype="fasta",
                    retmode="text"
                )
                batch_seqs = list(SeqIO.parse(handle, "fasta"))
                handle.close()

                sequences.extend(batch_seqs)
                print(f"    Downloaded {len(sequences)}/{len(id_list)}")

            except Exception as e:
                print(f"    Warning: Batch {i//batch_size + 1} failed: {e}")
                continue

        return sequences

    except Exception as e:
        print(f"    Error: {e}")
        return []


def main():
    """Download 100 diverse bacterial 16S rRNA sequences."""

    output_file = Path(__file__).parent.parent / "test_real_100species.fasta"

    print("=" * 80)
    print("DOWNLOADING 100 DIVERSE BACTERIAL 16S rRNA SEQUENCES")
    print("=" * 80)
    print()

    # Define queries for taxonomic diversity
    # Each tuple: (query, number_of_sequences, description)
    queries = [
        # ====================================================================
        # PROTEOBACTERIA (most diverse phylum) - 30 species
        # ====================================================================
        ("Escherichia[Organism] AND 16S ribosomal RNA[Title] AND type strain", 3,
         "Escherichia (type strains)"),

        ("Salmonella enterica[Organism] AND 16S ribosomal RNA[Title]", 3,
         "Salmonella enterica"),

        ("Pseudomonas[Organism] AND 16S ribosomal RNA[Title] AND type strain", 5,
         "Pseudomonas (type strains)"),

        ("Vibrio[Organism] AND 16S ribosomal RNA[Title] AND type strain", 3,
         "Vibrio (type strains)"),

        ("Klebsiella[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Klebsiella"),

        ("Helicobacter[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Helicobacter"),

        ("Rhizobium[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Rhizobium (Alpha)"),

        ("Caulobacter[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Caulobacter (Alpha)"),

        ("Neisseria[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Neisseria (Beta)"),

        ("Burkholderia[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Burkholderia (Beta)"),

        ("Desulfovibrio[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Desulfovibrio (Delta)"),

        ("Campylobacter[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Campylobacter (Epsilon)"),

        # ====================================================================
        # FIRMICUTES - 25 species
        # ====================================================================
        ("Bacillus[Organism] AND 16S ribosomal RNA[Title] AND type strain", 5,
         "Bacillus (type strains)"),

        ("Staphylococcus[Organism] AND 16S ribosomal RNA[Title]", 4,
         "Staphylococcus"),

        ("Streptococcus[Organism] AND 16S ribosomal RNA[Title]", 4,
         "Streptococcus"),

        ("Listeria[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Listeria"),

        ("Clostridium[Organism] AND 16S ribosomal RNA[Title]", 3,
         "Clostridium"),

        ("Enterococcus[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Enterococcus"),

        ("Lactobacillus[Organism] AND 16S ribosomal RNA[Title]", 3,
         "Lactobacillus"),

        ("Mycoplasma[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Mycoplasma"),

        # ====================================================================
        # ACTINOBACTERIA - 15 species
        # ====================================================================
        ("Mycobacterium[Organism] AND 16S ribosomal RNA[Title] AND type strain", 5,
         "Mycobacterium (type strains)"),

        ("Corynebacterium[Organism] AND 16S ribosomal RNA[Title]", 3,
         "Corynebacterium"),

        ("Streptomyces[Organism] AND 16S ribosomal RNA[Title]", 4,
         "Streptomyces"),

        ("Actinomyces[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Actinomyces"),

        ("Bifidobacterium[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Bifidobacterium"),

        # ====================================================================
        # BACTEROIDETES - 10 species
        # ====================================================================
        ("Bacteroides[Organism] AND 16S ribosomal RNA[Title]", 4,
         "Bacteroides"),

        ("Prevotella[Organism] AND 16S ribosomal RNA[Title]", 3,
         "Prevotella"),

        ("Flavobacterium[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Flavobacterium"),

        ("Porphyromonas[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Porphyromonas"),

        # ====================================================================
        # CYANOBACTERIA - 5 species
        # ====================================================================
        ("Synechococcus[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Synechococcus"),

        ("Anabaena[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Anabaena"),

        ("Nostoc[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Nostoc"),

        ("Prochlorococcus[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Prochlorococcus"),

        # ====================================================================
        # OTHER IMPORTANT PHYLA - 15 species
        # ====================================================================
        ("Chlamydia[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Chlamydia (Chlamydiae)"),

        ("Treponema[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Treponema (Spirochaetes)"),

        ("Borrelia[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Borrelia (Spirochaetes)"),

        ("Deinococcus[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Deinococcus (Deinococcus-Thermus)"),

        ("Thermus[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Thermus (Deinococcus-Thermus)"),

        ("Aquifex[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Aquifex (Aquificae)"),

        ("Thermotoga[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Thermotoga (Thermotogae)"),

        ("Fusobacterium[Organism] AND 16S ribosomal RNA[Title]", 2,
         "Fusobacterium (Fusobacteria)"),

        ("Chlorobium[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Chlorobium (Chlorobi)"),

        ("Chloroflexus[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Chloroflexus (Chloroflexi)"),

        ("Planctomyces[Organism] AND 16S ribosomal RNA[Title]", 1,
         "Planctomyces (Planctomycetes)"),
    ]

    print(f"Target: {sum(q[1] for q in queries)} sequences across {len(queries)} taxonomic groups")
    print()

    all_sequences = []
    stats_by_phylum = defaultdict(int)

    print("Downloading sequences by taxonomic group:")
    print("-" * 80)

    for query, max_seqs, description in queries:
        print(f"\n{description}:")

        seqs = download_sequences_for_taxon(query, max_seqs, delay=0.5)

        if seqs:
            all_sequences.extend(seqs)
            # Approximate phylum from description
            phylum = description.split('(')[-1].rstrip(')') if '(' in description else "Other"
            stats_by_phylum[phylum] += len(seqs)

        print(f"    Total collected so far: {len(all_sequences)}")

        time.sleep(1)  # Be respectful to NCBI servers

    print()
    print("=" * 80)
    print("DOWNLOAD COMPLETE!")
    print("=" * 80)
    print()

    # Write sequences
    print(f"Writing {len(all_sequences)} sequences to: {output_file}")
    SeqIO.write(all_sequences, output_file, "fasta")

    print()
    print("Dataset Summary:")
    print(f"  Total sequences: {len(all_sequences)}")
    print()

    # Sequence statistics
    lengths = [len(seq.seq) for seq in all_sequences]
    print("  Sequence lengths:")
    print(f"    Min: {min(lengths)} bp")
    print(f"    Max: {max(lengths)} bp")
    print(f"    Mean: {sum(lengths) / len(lengths):.0f} bp")
    print()

    print("  Sequences by group:")
    for phylum, count in sorted(stats_by_phylum.items(), key=lambda x: -x[1]):
        print(f"    {phylum}: {count}")
    print()

    print(f"Output file: {output_file}")
    print()
    print("Next step: Run performance test with this dataset")
    print()


if __name__ == "__main__":
    main()

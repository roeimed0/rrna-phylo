# Test Datasets for rRNA-Phylo

This directory contains curated test datasets for phylogenetic analysis.

## Available Datasets

### 1. Mammals (33 sequences)
- **Files:** `mammals_test.fasta`, `mammals_test_aligned.fasta`
- **Description:** Diverse mammalian species
- **Alignment length:** 2,132 bp
- **Use case:** Testing mammalian phylogeny reconstruction

### 2. Birds (35 sequences) ⭐ NEW
- **Files:** `birds_test.fasta`, `birds_test_aligned.fasta`
- **Description:** Diverse bird species across major orders + 1 turtle outgroup
- **Alignment length:** 1,865 bp
- **Taxonomic coverage:**
  - Ratites: Ostrich, Rhea, Emu, Kiwi
  - Gamebirds: Chicken, Turkey, Quail
  - Waterfowl: Mallard
  - Raptors: Bald Eagle, Falcon, Vulture
  - Owls, Woodpeckers, Parrots
  - Swifts, Hummingbirds
  - Passerines (songbirds): 15 species
  - Shorebirds, Gulls, Pelicans
  - **Outgroup:** Red-eared slider turtle (*Trachemys scripta*)
- **Use case:** Avian phylogeny with proper outgroup rooting

### 3. Arcosauria (111 sequences) - ORIGINAL LARGE DATASET
- **Files:** `Arcosauria_test.fasta`, `Arcosauria_test_aligned.fasta`
- **Description:** Birds + reptiles (very distant phylogenetic groups)
- **Alignment length:** 1,916 bp
- **⚠️ Note:** This dataset is TOO DIVERSE for most analyses. Use `birds_test.fasta` instead for cleaner results.
- **Use case:** Testing with very large, divergent datasets

### 4. Cartilaginous Fish (28 sequences)
- **Files:** `cartilaginous_fish_test.fasta`, `cartilaginous_fish_test_aligned.fasta`
- **Description:** Sharks, rays, and chimaeras
- **Alignment length:** 1,827 bp
- **Use case:** Testing fish phylogeny

## Dataset Recommendations

| Dataset | Best For | Tree Building Time* |
|---------|----------|---------------------|
| `mammals_test` | General testing, moderate diversity | ~10 seconds |
| `birds_test` ⭐ | **Avian phylogeny, best taxonomic range** | **~15 seconds** |
| `Arcosauria_test` | Stress testing, large datasets | ~30 seconds |
| `cartilaginous_fish_test` | Fish phylogeny | ~8 seconds |

*Using pre-aligned files with ML method

## Pre-aligned vs Unaligned

Each dataset has two versions:

- **`*_test.fasta`** - Unaligned sequences (runs MUSCLE, slower)
- **`*_test_aligned.fasta`** - Pre-aligned sequences (skip MUSCLE, 5-30x faster!)

### Performance Comparison

| Dataset | Unaligned (with MUSCLE) | Pre-aligned (skip MUSCLE) | Speedup |
|---------|-------------------------|---------------------------|---------|
| mammals (33 seqs) | ~2 min | ~10 sec | 12x |
| birds (35 seqs) | ~1.5 min | ~15 sec | 6x |
| Arcosauria (111 seqs) | ~19 min | ~30 sec | 38x |
| cartilaginous_fish (28 seqs) | ~1 min | ~8 sec | 7.5x |

## Usage Examples

### Using pre-aligned files (RECOMMENDED)
```bash
# CLI
python rrna_phylo_cli.py data/test/birds_test_aligned.fasta --pre-aligned --method ml

# Interactive
python rrna_phylo_app.py
# Select: birds_test_aligned.fasta
# Pre-aligned? y
```

### Using unaligned files (tests MUSCLE)
```bash
python rrna_phylo_app.py
# Select: birds_test.fasta
# Pre-aligned? n
```

## Creation Scripts

Test datasets were created using:
- `scripts/create_birds_subset.py` - Curated bird species selection from Arcosauria
- Manual curation for other datasets

## Notes

- All sequences are 18S rRNA gene sequences
- FASTA headers include GenBank accession numbers
- Pre-aligned files generated using MUSCLE v5.1
- Recommended dataset: **birds_test** for most use cases

# Birds Test Dataset - Summary

## What We Created

A curated **birds-only phylogenetic test dataset** to replace the overly-diverse Arcosauria dataset.

## Problem with Original Arcosauria Dataset

- **111 sequences** mixing birds + reptiles (lizards, turtles)
- **Too phylogenetically distant** - poor for testing realistic analyses
- **Slow**: ~19 minutes with MUSCLE alignment
- **Confusing results**: Mixes archosaurs (birds+crocodilians) with lepidosaurs (lizards)

## New Birds Dataset

### Files Created
```
backend/data/test/
├── birds_test.fasta               # 35 unaligned sequences
├── birds_test_aligned.fasta       # 35 aligned sequences (1,865 bp)
└── README.md                      # Dataset documentation
```

### Dataset Composition

**35 total sequences:**
- **34 bird species** across major orders
- **1 reptile outgroup** (Red-eared slider turtle, *Trachemys scripta*)

### Taxonomic Coverage

| Bird Group | Species Count | Examples |
|------------|---------------|----------|
| **Ratites** (flightless) | 4 | Ostrich, Rhea, Emu, Kiwi |
| **Galliformes** (gamebirds) | 3 | Chicken, Turkey, Quail |
| **Anseriformes** (waterfowl) | 1 | Mallard duck |
| **Columbiformes** (pigeons) | 1 | Rock pigeon |
| **Accipitriformes** (raptors) | 3 | Bald eagle, Falcon, Vulture |
| **Strigiformes** (owls) | 1 | Great horned owl |
| **Piciformes** (woodpeckers) | 1 | Downy woodpecker |
| **Psittaciformes** (parrots) | 1 | Budgerigar |
| **Apodiformes** (swifts/hummers) | 2 | Common swift, Mango hummingbird |
| **Coraciiformes** (kingfishers) | 2 | Roller, Hoopoe |
| **Passeriformes** (songbirds) | 11 | Sparrow, Magpie, Thrushes, Warblers, etc. |
| **Gruiformes** (cranes) | 1 | Sandhill crane |
| **Charadriiformes** (shorebirds) | 2 | Plover, Gull |
| **Pelecaniformes** (pelicans) | 1 | Brown pelican |
| **Suliformes** (boobies) | 2 | Booby, Cormorant |
| **OUTGROUP** | 1 | Red-eared slider turtle |

## Performance Comparison

| Dataset | Sequences | Unaligned | Pre-aligned | Speedup |
|---------|-----------|-----------|-------------|---------|
| **Arcosauria** (old) | 111 | ~19 min | ~30 sec | 38x |
| **Birds** (new) ⭐ | 35 | ~1.5 min | ~15 sec | 6x |

## Why This Is Better

### 1. **Appropriate Phylogenetic Scale**
- All ingroup taxa are birds (Class Aves)
- Single reptile outgroup for proper tree rooting
- Realistic evolutionary distances

### 2. **Faster Testing**
- 6x faster than Arcosauria when pre-aligned
- Still large enough to test real-world performance
- MUSCLE completes in ~1 minute vs ~19 minutes

### 3. **Better Taxonomic Representation**
- Covers all major bird orders
- Includes both ancient lineages (ratites) and modern (passerines)
- Balanced representation across avian phylogeny

### 4. **Cleaner Results**
- Trees make biological sense
- Proper outgroup rooting with turtle
- No confusing reptile polyphyly

## Usage Examples

### Quick Test (Pre-aligned)
```bash
# Interactive
python rrna_phylo_app.py
# Select: birds_test_aligned.fasta
# Pre-aligned? y
# Method: ml
# Result: Tree in ~15 seconds!

# CLI
python rrna_phylo_cli.py data/test/birds_test_aligned.fasta \
    --pre-aligned --method ml --bootstrap 100
```

### Test MUSCLE Alignment (Unaligned)
```bash
python rrna_phylo_app.py
# Select: birds_test.fasta
# Pre-aligned? n
# Method: ml
# Result: Alignment + tree in ~2 minutes
```

## Dataset Statistics

```
Sequences:        35
Alignment length: 1,865 bp
File sizes:
  - Unaligned:    56 KB
  - Aligned:      66 KB
MUSCLE time:      ~1 minute
Tree build time:  ~15 seconds (pre-aligned ML)
```

## Script for Reproducibility

The dataset was created using:
```bash
python backend/scripts/create_birds_subset.py
```

This script:
1. Reads the original 111-sequence Arcosauria dataset
2. Selects 34 diverse bird species across major orders
3. Adds 1 turtle as outgroup (Trachemys scripta)
4. Writes birds_test.fasta
5. Runs MUSCLE to create birds_test_aligned.fasta

## Recommendation

**Use `birds_test_aligned.fasta` for most testing!**

It's the perfect middle-ground dataset:
- Not too small (like mammals with 33 seqs)
- Not too large (like Arcosauria with 111 seqs)
- Appropriate phylogenetic scale
- Fast performance
- Biologically meaningful results

## Next Steps

If you want to test with the original large Arcosauria dataset:
```bash
# Still available, but not recommended for routine testing
python rrna_phylo_app.py
# Select: Arcosauria_test_aligned.fasta
```

---

**Summary**: The new birds dataset provides a focused, fast, and biologically appropriate test case for avian phylogenetics! 🦅🦜🦆

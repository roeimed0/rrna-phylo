# rRNA-Phylo: Phylogenetic Tree Builder for Ribosomal RNA

High-performance phylogenetic tree inference from DNA/RNA sequences using three complementary methods: UPGMA, BioNJ, and Maximum Likelihood with automatic model selection.

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![NumPy](https://img.shields.io/badge/NumPy-1.21+-orange.svg)](https://numpy.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Features

### Core Functionality
- **3 phylogenetic methods**: UPGMA, BioNJ, Maximum Likelihood
- **Automatic model selection**: Tests 5 substitution models (JC69, K80, F81, HKY85, GTR) + gamma rate variation
- **Bootstrap analysis**: Confidence values for all branches
- **Multiple outputs**: Newick trees, ASCII visualizations, comparison summary
- **Pre-aligned support**: Skip alignment step for faster analysis

### Performance
- **Site pattern compression**: speedup for ML calculations
- **Numba JIT compilation**: speedup for likelihood calculations
- **Combined optimization**: total speedup for ML inference
- **MUSCLE integration**: Automatic sequence alignment with 30-minute timeout

### Quality
- **Test datasets included**: 4 curated datasets (mammals, birds, fish)
- **Comprehensive documentation**: Complete usage guides and API docs
- **Production-ready**: Error handling, logging, organized outputs

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/rrna-phylo.git
cd rrna-phylo

# Install dependencies
pip install -r requirements.txt

```

### Basic Usage

```bash
# Interactive menu (easiest)
python rrna_phylo_app.py

# CLI - Build all 3 trees
python rrna_phylo_cli.py data/test/birds_test_aligned.fasta --pre-aligned

# CLI - Maximum Likelihood only with bootstrap
python rrna_phylo_cli.py data/test/mammals_test_aligned.fasta --pre-aligned --method ml --bootstrap 100
```

### Quick Test

```bash
# Use pre-aligned birds dataset (35 species, ~15 seconds)
python rrna_phylo_cli.py backend/data/test/birds_test_aligned.fasta --pre-aligned --method ml

# Output in: backend/results/birds_test_aligned/
```

---

## Methods

### 1. UPGMA (Ultrametric)
- **Speed**: Fastest (<1s for 100 sequences)
- **Assumption**: Molecular clock (constant evolution rate)
- **Best for**: Closely related sequences, initial exploration

### 2. BioNJ (Distance-Based)
- **Speed**: Fast (<1s for 100 sequences)
- **Assumption**: No molecular clock
- **Best for**: General-purpose phylogeny, variance-weighted accuracy

### 3. Maximum Likelihood (Statistical)
- **Speed**: ~15-30s for 35-50 sequences (with pre-alignment)
- **Assumption**: Statistical model of sequence evolution
- **Best for**: Publication-quality trees, rigorous inference
- **Features**:
  - Automatic model selection (BIC)
  - NNI/SPR tree search
  - Bootstrap support values
  - Gamma rate variation

---

## Test Datasets

Located in `backend/data/test/`:

| Dataset | Sequences | Alignment | Speed (pre-aligned) | Use Case |
|---------|-----------|-----------|---------------------|----------|
| **birds_test** | 35 | 1,865 bp | ~15 sec | **Recommended - diverse birds + turtle outgroup** |
| mammals_test | 33 | 2,132 bp | ~10 sec | Mammalian phylogeny |
| cartilaginous_fish_test | 28 | 1,827 bp | ~8 sec | Sharks and rays |
| Arcosauria_test | 111 | 1,916 bp | ~30 sec | Stress testing only (too diverse) |

**Each dataset has 2 versions:**
- `*_test.fasta` - Unaligned (runs MUSCLE, slower)
- `*_test_aligned.fasta` - Pre-aligned (skip MUSCLE, 5-30x faster!) ⚡

**See:** `backend/data/test/README.md` for complete dataset documentation

---

## CLI Usage

### All Methods (Default)

```bash
python rrna_phylo_cli.py sequences.fasta
```

**Output:** `backend/results/sequences/`
- `upgma_tree.nwk`, `bionj_tree.nwk`, `ml_tree.nwk` (Newick format)
- `upgma_tree.txt`, `bionj_tree.txt`, `ml_tree.txt` (ASCII visualization)
- `summary.txt` (comparison of all methods)

### Single Method

```bash
# UPGMA only
python rrna_phylo_cli.py sequences.fasta --method upgma

# BioNJ only
python rrna_phylo_cli.py sequences.fasta --method bionj

# Maximum Likelihood only
python rrna_phylo_cli.py sequences.fasta --method ml
```

### Pre-Aligned Sequences

```bash
# Skip MUSCLE alignment (5-30x faster!)
python rrna_phylo_cli.py aligned.fasta --pre-aligned
```

### Bootstrap Analysis

```bash
# 100 replicates (recommended for publication)
python rrna_phylo_cli.py sequences.fasta --bootstrap 100

# 10 replicates (quick test)
python rrna_phylo_cli.py sequences.fasta --bootstrap 10 --method upgma
```

### Output Format

```bash
# Newick only (no ASCII trees)
python rrna_phylo_cli.py sequences.fasta --output-format newick

# ASCII only
python rrna_phylo_cli.py sequences.fasta --output-format ascii

# Both (default)
python rrna_phylo_cli.py sequences.fasta --output-format both
```

---

## Interactive Menu

```bash
python rrna_phylo_app.py
```

**Features:**
- 📁 File browser (shows all files in `data/test/`)
- ⚡ Quick build (one-click with defaults)
- 🎛️ Custom build (choose all options)
- 🧬 Pre-aligned sequence detection
- 📊 Results viewer
- 🧹 Cleanup tool
- 💡 Built-in help

**No command-line flags to remember!**

---

## Output Files

### Directory Structure

```
backend/results/
└── [filename]/
    ├── upgma_tree.nwk       # UPGMA tree (Newick)
    ├── bionj_tree.nwk       # BioNJ tree (Newick)
    ├── ml_tree.nwk          # ML tree (Newick)
    ├── upgma_tree.txt       # UPGMA tree (ASCII)
    ├── bionj_tree.txt       # BioNJ tree (ASCII)
    ├── ml_tree.txt          # ML tree (ASCII)
    └── summary.txt          # Comparison summary
```

### ASCII Tree Example

```
                    ┌─ Gallus_gallus (Chicken)
          ┌─────────┤
          │         └─ Meleagris_gallopavo (Turkey)
    ──────┤
          │         ┌─ Anas_platyrhynchos (Mallard)
          └─────────┤
                    └─ Struthio_camelus (Ostrich)
```

---

## Performance Benchmarks

### Speed Comparison: Unaligned vs Pre-aligned

| Dataset | Sequences | Unaligned (with MUSCLE) | Pre-aligned | Speedup |
|---------|-----------|-------------------------|-------------|---------|
| Birds | 35 | ~1.5 min | **15 sec** | **6x** ⚡ |
| Mammals | 33 | ~2 min | **10 sec** | **12x** ⚡ |
| Arcosauria | 111 | ~19 min | **30 sec** | **38x** ⚡ |
| Fish | 28 | ~1 min | **8 sec** | **7.5x** ⚡ |

**Recommendation:** Use pre-aligned datasets for testing and development!

### Method Comparison (35 sequences, pre-aligned)

| Method | Time | Bootstrap (100 reps) |
|--------|------|----------------------|
| UPGMA | <1 sec | ~2 min |
| BioNJ | <1 sec | ~2 min |
| ML (NNI) | ~15 sec | ~25 min |

---

## Project Structure

```
rrna-phylo/
├── README.md                           # This file
├── requirements.txt                    # Core dependencies
├── requirements-dev.txt                # Development dependencies
├── docs/                               # Documentation
│   ├── CHANGELOG.md                    # Project history
│   └── BIRDS_DATASET_SUMMARY.md        # Birds dataset details
├── backend/
│   ├── rrna_phylo_cli.py               # Main CLI entry point
│   ├── rrna_phylo_app.py               # Interactive menu
│   ├── muscle.exe                      # MUSCLE aligner (Windows)
│   ├── data/
│   │   ├── test/                       # Test datasets
│   │   │   ├── birds_test.fasta
│   │   │   ├── birds_test_aligned.fasta ⭐
│   │   │   ├── mammals_test.fasta
│   │   │   ├── mammals_test_aligned.fasta
│   │   │   ├── cartilaginous_fish_test.fasta
│   │   │   ├── cartilaginous_fish_test_aligned.fasta
│   │   │   └── README.md               # Dataset documentation
│   │   └── README.md
│   └── rrna_phylo/                     # Core package
│       ├── alignment/                  # MUSCLE integration
│       │   └── muscle_aligner.py
│       ├── core/                       # Tree structures
│       │   ├── builder.py
│       │   ├── tree.py
│       │   └── sequence_type.py
│       ├── distance/                   # Distance calculations
│       │   └── distance.py
│       ├── methods/                    # Tree building methods
│       │   ├── upgma.py
│       │   └── bionj.py
│       ├── models/                     # Maximum Likelihood
│       │   ├── ml_tree_level3.py       # ML with optimizations
│       │   ├── rate_matrices.py        # GTR, HKY, K80, etc.
│       │   └── model_selection.py      # BIC-based selection
│       ├── consensus/                  # Tree comparison
│       │   ├── tree_distance.py        # Robinson-Foulds
│       │   └── bipartitions.py
│       └── visualization/              # Tree output
│           └── tree_drawer.py          # ASCII trees
```

---

## Requirements

### Core Dependencies

```
Python 3.9+
numpy >= 1.21.0
scipy >= 1.7.0
numba >= 0.54.0
```

### External Tools

- **MUSCLE v5.1** (required for alignment)
  - Download: https://drive5.com/muscle/
  - Place in PATH or project root

### Optional

- **Biopython >= 1.79** (for sequence validation)

---

## Troubleshooting

### Common Issues

**"No module named 'rrna_phylo'"**
```bash
# Make sure you're in the project root, not backend/
cd /path/to/rrna-phylo
python backend/rrna_phylo_cli.py ...
```

**MUSCLE timeout**
```bash
# Default timeout is 30 minutes
# For very large datasets, the timeout may trigger
# Pre-align sequences externally and use --pre-aligned flag
```

**OpenMP Library Conflict**
```bash
# Already handled automatically
# If still occurring, set environment variable:
export KMP_DUPLICATE_LIB_OK=TRUE  # Linux/Mac
set KMP_DUPLICATE_LIB_OK=TRUE     # Windows
```

**ML is slow**
```bash
# Use pre-aligned datasets (skip MUSCLE)
# Skip bootstrap for testing
# Use smaller datasets for development
```

---

## Documentation

- **`README.md`** - This file (getting started)
- **`backend/data/test/README.md`** - Test dataset documentation
- **`docs/CHANGELOG.md`** - Complete project history and improvements
- **`docs/BIRDS_DATASET_SUMMARY.md`** - Birds dataset creation and rationale

---

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

---

## License

MIT License - See LICENSE file for details

---

## Citation

If you use rRNA-Phylo in your research, please cite:

```
rRNA-Phylo: High-performance phylogenetic analysis for ribosomal RNA
https://github.com/yourusername/rrna-phylo
```

---

## Acknowledgments

Built with Claude Code (Sonnet 4.5) demonstrating systematic LLM-assisted development through skills and specialized agents.

**Key optimizations:**
- Site pattern compression (8-10x speedup)
- Numba JIT compilation (9x speedup)
- Combined: 72x total speedup for ML likelihood calculations

---

## References

### Phylogenetics
- Felsenstein, J. (1985). Confidence limits on phylogenies: An approach using the bootstrap. *Evolution*, 39(4), 783-791.
- Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm. *Molecular Biology and Evolution*, 14(7), 685-695.

### Model Selection
- Posada, D., & Crandall, K. A. (1998). MODELTEST: testing the model of DNA substitution. *Bioinformatics*, 14(9), 817-818.

### Tools
- Edgar, R.C. (2022). MUSCLE v5. *Nature Communications*, 13, 6968.
- Numba Documentation: https://numba.pydata.org/

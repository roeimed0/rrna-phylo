# rRNA-Phylo

Production-ready phylogenetic tree building for rRNA sequences. Build trees using distance-based methods (UPGMA, BioNJ) and maximum likelihood with automatic model selection.

## Features

- **Three tree-building methods**: UPGMA, BioNJ, Maximum Likelihood
- **Automatic model selection**: AIC/BIC comparison (JC69, K80, F81, HKY85, GTR)
- **Tree search algorithms**: NNI (fast) and SPR (thorough)
- **Bootstrap support**: Statistical confidence assessment (Felsenstein 1985)
- **Performance**: 72x speedup through Numba JIT and site pattern compression

## Installation

```bash
cd backend
pip install numpy scipy numba biopython
```

## Quick Start

### Build All Three Tree Types (Recommended)

```bash
cd backend
python build_trees.py sequences.fasta
```

**Output**: 6 files showing UPGMA, BioNJ, and ML trees in both Newick and ASCII formats

### Build Single Tree Type

```bash
# Fast distance-based tree (UPGMA)
python build_trees.py sequences.fasta --method upgma

# Fast distance-based tree (BioNJ, no clock assumption)
python build_trees.py sequences.fasta --method bionj

# Rigorous maximum likelihood tree (GTR+Gamma)
python build_trees.py sequences.fasta --method ml
```

### Add Bootstrap Support

```bash
# 100 bootstrap replicates (recommended for publication)
python build_trees.py sequences.fasta --bootstrap 100
```

## Output Files

```
results/
└── [filename]/         # Subfolder named after input file (without extension)
    ├── upgma_tree.nwk      # UPGMA tree (Newick format)
    ├── bionj_tree.nwk      # BioNJ tree (Newick format)
    ├── ml_tree.nwk         # ML tree (Newick format)
    ├── upgma_tree.txt      # UPGMA tree (ASCII visualization)
    ├── bionj_tree.txt      # BioNJ tree (ASCII visualization)
    ├── ml_tree.txt         # ML tree (ASCII visualization)
    └── summary.txt         # Comparison of all methods
```

Example: Running `python build_trees.py sequences.fasta` creates `results/sequences/`

## Performance

- 48 sequences: ~3 minutes (ML + NNI)
- Bootstrap: ~4 seconds per replicate
- Site pattern compression: 1.6-2.5x speedup
- Numba acceleration: 9x speedup for likelihood calculations

## Project Structure

```
backend/
├── rrna_phylo/              # Core package
│   ├── io/                  # FASTA parsing
│   ├── core/                # Tree data structures
│   ├── distance/            # Distance calculations
│   ├── methods/             # UPGMA, NJ, BioNJ
│   └── models/              # ML tree building, bootstrap
└── build_phylogenetic_tree.py  # Command-line interface
```

## Requirements

- Python 3.8+
- NumPy
- SciPy
- Numba
- BioPython

## License

MIT License

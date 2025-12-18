# rRNA-Phylo

Production-ready phylogenetic tree building from rRNA sequences using maximum likelihood methods.

## Features

- **Distance-based methods**: UPGMA, Neighbor-Joining, BioNJ
- **Maximum Likelihood**: GTR+Gamma model with automatic model selection (AIC/BIC)
- **Tree search**: NNI and SPR algorithms
- **Bootstrap support**: Statistical confidence assessment (Felsenstein 1985)
- **Performance**: 72x speedup through Numba JIT and site pattern compression

## Quick Start

```bash
# Install dependencies
cd backend
pip install numpy scipy numba biopython

# Build a tree with bootstrap support
python build_phylogenetic_tree.py sequences.fasta --method nni --bootstrap 100
```

## Output

The tool generates:
- `tree.nwk` - Newick format tree with bootstrap values
- `tree_ascii.txt` - Human-readable visualization
- `metadata.json` - Build statistics and parameters
- `build_log.txt` - Detailed execution log

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

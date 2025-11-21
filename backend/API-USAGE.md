# rRNA-Phylo API Usage Guide

## Installation

The package is now installed in development mode:
```bash
cd backend
pip install -e .
```

## Quick Start

### Simple API - One Line Tree Building

```python
from rrna_phylo import build_trees, Sequence

# Create sequences
sequences = [
    Sequence('human', 'ATCGATCGATCG...', 'Homo sapiens'),
    Sequence('chimp', 'ATCGATCGATCG...', 'Pan troglodytes'),
    Sequence('mouse', 'ATCGTTTTTTCG...', 'Mus musculus'),
]

# Build all three types of trees automatically
upgma_tree, bionj_tree, ml_tree = build_trees(sequences)

# Get Newick format
print(upgma_tree.to_newick())
```

### Builder Class API

```python
from rrna_phylo import PhylogeneticTreeBuilder, Sequence

# Create builder
builder = PhylogeneticTreeBuilder(verbose=True)

# Detect sequence type automatically
seq_type = builder.detect_and_validate(sequences)

# Build individual trees
upgma = builder.build_upgma_tree(sequences)
bionj = builder.build_bionj_tree(sequences)
ml = builder.build_ml_tree(sequences, alpha=1.0)

# Or build all at once
upgma, bionj, ml = builder.build_all_trees(sequences)
```

## Direct Method Calls

### UPGMA (Molecular Clock)

```python
from rrna_phylo import build_upgma_tree
import numpy as np

# From distance matrix
distance_matrix = np.array([[0.0, 0.2, 0.6],
                           [0.2, 0.0, 0.6],
                           [0.6, 0.6, 0.0]])
labels = ['A', 'B', 'C']
tree = build_upgma_tree(distance_matrix, labels)
```

### BioNJ (No Clock Assumption)

```python
from rrna_phylo import build_bionj_tree
import numpy as np

tree = build_bionj_tree(distance_matrix, labels)
```

### Maximum Likelihood

```python
from rrna_phylo import build_ml_tree_level3, Sequence

# DNA/RNA sequences use GTR+Gamma
sequences = [...]
tree, log_likelihood = build_ml_tree_level3(sequences, alpha=1.0)

# Protein sequences automatically use WAG/LG+Gamma
protein_seqs = [...]
from rrna_phylo import build_protein_ml_tree
tree, log_likelihood = build_protein_ml_tree(protein_seqs, model_name="WAG", alpha=1.0)
```

## Distance Calculations

```python
from rrna_phylo import calculate_distance_matrix, calculate_protein_distance_matrix, Sequence

# DNA/RNA distance (Jukes-Cantor or Kimura)
sequences = [...]
dist_matrix, labels = calculate_distance_matrix(sequences, model="jukes_cantor")

# Protein distance (Poisson or Kimura)
protein_seqs = [...]
dist_matrix, labels = calculate_protein_distance_matrix(protein_seqs, model="poisson")
```

## Sequence Type Detection

```python
from rrna_phylo import SequenceType, SequenceTypeDetector, Sequence

detector = SequenceTypeDetector()
sequences = [...]

# Detect type
seq_type = detector.detect_sequences(sequences)

if seq_type == SequenceType.DNA:
    print("DNA sequences detected")
elif seq_type == SequenceType.RNA:
    print("RNA sequences detected")
elif seq_type == SequenceType.PROTEIN:
    print("Protein sequences detected")
```

## Tree Operations

```python
from rrna_phylo import TreeNode

# Create tree
tree = build_upgma_tree(distance_matrix, labels)

# Get Newick format
newick_str = tree.to_newick()

# Check if leaf
if tree.is_leaf():
    print(f"Leaf: {tree.name}")

# Access children
if tree.left:
    left_newick = tree.left.to_newick()
if tree.right:
    right_newick = tree.right.to_newick()
```

## Package Structure

```
rrna_phylo/
├── core/              # Core classes
│   ├── tree.py        # TreeNode
│   ├── sequence_type.py  # Sequence detection
│   └── builder.py     # PhylogeneticTreeBuilder, build_trees()
├── io/                # Input/Output
│   ├── fasta_parser.py   # Sequence, FastaParser
│   └── aligner.py     # Multiple sequence alignment
├── distance/          # Distance calculations
│   ├── distance.py    # DNA/RNA distances
│   └── protein_distance.py  # Protein distances
├── models/            # Substitution models
│   ├── ml_tree_level3.py    # GTR+Gamma (DNA/RNA)
│   └── protein_models.py    # WAG/LG/JTT models
├── methods/           # Tree building methods
│   ├── upgma.py       # UPGMA
│   ├── bionj.py       # BioNJ
│   └── protein_ml.py  # Protein ML
└── utils/             # Utilities
    └── visualize_trees.py  # Tree visualization
```

## Examples

See the [examples/](examples/) directory for more usage examples:
- `compare_trees.py` - Compare UPGMA vs BioNJ on same data

## Testing

```bash
cd backend

# Run individual tests
python tests/test_upgma.py
python tests/test_phylo_builder.py
python tests/test_sequence_type.py

# Test all
python -m pytest tests/
```

## Notes

- Package automatically detects sequence type (DNA/RNA/Protein)
- Automatically selects appropriate substitution model
- GTR+Gamma for DNA/RNA sequences
- WAG/LG/JTT+Gamma for protein sequences
- All functionality works - Unicode display errors are cosmetic only

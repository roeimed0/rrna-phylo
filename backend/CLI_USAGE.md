# rRNA-Phylo Command-Line Interface (CLI)

## Overview

The rRNA-Phylo CLI provides an easy way to build phylogenetic trees from FASTA files without writing any code.

## Installation

```bash
cd backend
pip install -r requirements.txt
```

## Basic Usage

### 1. Prepare Your Input

Create a FASTA file with **aligned** sequences (minimum 3 sequences):

```fasta
>seq1 Escherichia coli
ATGCTACCGGGCCCATACCCCAAACATGTTGGTTATACCCCTTCCCGTACTAATAAAC
>seq2 Salmonella enterica
ATGCTGCCGGGCCCATGCCCCAAACATGTTGGTCATACCCCTTCCCGTACTAATAAAC
>seq3 Bacillus subtilis
ATGCTGCCGGGCCCATACCCCAAACATGTTGGTTACACCCCTTCCCGTACTAATAAAC
>seq4 Staphylococcus aureus
ATGCTGCCAGGCCCACACCCCCAGCATGTTGGTCACACCCCTTCCCGTACTAATAAAC
```

**Important**: Sequences MUST be aligned (same length) before phylogenetic analysis!

### 2. Run the CLI

```bash
# Build all three trees (ML, BioNJ, UPGMA)
python rrna-phylo.py sequences.fasta

# Or if installed system-wide
rrna-phylo sequences.fasta
```

### 3. Get Results

The CLI will:
1. Display ASCII trees in the terminal
2. Export Newick files (.nwk) for external visualization

```
Output files:
  - tree_ml.nwk
  - tree_bionj.nwk
  - tree_upgma.nwk
```

## Common Usage Examples

### Build Only One Method

```bash
# Maximum Likelihood only
python rrna-phylo.py sequences.fasta --method ml

# BioNJ only
python rrna-phylo.py sequences.fasta --method bionj

# UPGMA only
python rrna-phylo.py sequences.fasta --method upgma
```

### Add Bootstrap Support

```bash
# 50 bootstrap replicates (fast, good for exploration)
python rrna-phylo.py sequences.fasta --bootstrap 50

# 100 bootstrap replicates (standard for publications)
python rrna-phylo.py sequences.fasta --bootstrap 100

# Bootstrap with only ML method
python rrna-phylo.py sequences.fasta --method ml --bootstrap 50
```

### Specify Output Directory

```bash
# Save to results/ directory
python rrna-phylo.py sequences.fasta -o results/

# Custom prefix for output files
python rrna-phylo.py sequences.fasta --prefix my_analysis
# Output: my_analysis_ml.nwk, my_analysis_bionj.nwk, my_analysis_upgma.nwk
```

### Output Format Options

```bash
# Only Newick files (no ASCII visualization)
python rrna-phylo.py sequences.fasta --output-format newick

# Only ASCII visualization (no Newick files)
python rrna-phylo.py sequences.fasta --output-format ascii

# Both (default)
python rrna-phylo.py sequences.fasta --output-format both

# Quiet mode (suppress terminal output)
python rrna-phylo.py sequences.fasta --quiet
```

## Complete Command Reference

```
usage: rrna-phylo [-h] [-o OUTPUT] [-m {all,ml,bionj,upgma}]
                  [-b BOOTSTRAP] [-f {newick,ascii,both}] [-q]
                  [--prefix PREFIX] [--no-support] [-v]
                  input

positional arguments:
  input                    Input FASTA file with aligned sequences

options:
  -h, --help              Show help message and exit
  -o OUTPUT, --output     Output directory (default: current directory)
  -m METHOD, --method     Tree building method: all, ml, bionj, upgma (default: all)
  -b N, --bootstrap N     Number of bootstrap replicates (default: 0)
  -f FORMAT              Output format: newick, ascii, both (default: both)
  -q, --quiet            Suppress ASCII tree visualization
  --prefix PREFIX        Output file prefix (default: tree)
  --no-support           Exclude bootstrap support from Newick output
  -v, --verbose          Verbose output (show progress)
```

## Visualizing Results

### Option 1: ASCII Trees (Built-in)

The CLI displays trees directly in your terminal:

```
ML Tree:
----------------------------------------
`-- Internal (dist: 0.0000)
    |-- Internal (dist: 0.0384)
    |   |-- Internal (dist: 0.0133)
    |   |   |-- seq1 (dist: 0.0273)
    |   |   `-- seq2 (dist: 0.0263)
    |   `-- seq3 (dist: 0.0166)
    `-- seq4 (dist: 0.0384)
```

### Option 2: External Tools (Recommended for Publication)

Use the exported Newick files (.nwk) with:

1. **[FigTree](http://tree.bio.ed.ac.uk/software/figtree/)** (Desktop)
   - Download and install FigTree
   - Open your `.nwk` file
   - Customize visualization

2. **[iTOL](https://itol.embl.de/)** (Web-based)
   - Upload your `.nwk` file to iTOL
   - Interactive visualization
   - Publication-quality export

3. **[Dendroscope](http://dendroscope.org/)** (Desktop)
   - Advanced tree visualization
   - Multiple tree comparison

4. **[MEGA](https://www.megasoftware.net/)** (Desktop)
   - Comprehensive phylogenetic software
   - Built-in tree visualization

## Real-World Examples

### Example 1: Bacterial 16S rRNA Analysis

```bash
# Download bacterial 16S sequences from NCBI
# Align them with MUSCLE or MAFFT
# Build trees with bootstrap
python rrna-phylo.py bacterial_16s_aligned.fasta \
  --bootstrap 100 \
  --output-format newick \
  --prefix bacterial_16s \
  -o results/
```

### Example 2: Protein Phylogeny

```bash
# RecA protein sequences
python rrna-phylo.py reca_proteins.fasta \
  --method ml \
  --bootstrap 50 \
  --prefix reca_phylogeny
```

### Example 3: Quick Exploration

```bash
# Fast analysis without bootstrap
python rrna-phylo.py my_sequences.fasta \
  --method bionj \
  --output-format ascii
```

### Example 4: Publication-Quality Analysis

```bash
# Full analysis with all methods and high bootstrap
python rrna-phylo.py my_alignment.fasta \
  --bootstrap 1000 \
  --prefix publication \
  -o publication_results/
```

## Input File Requirements

### 1. FASTA Format

Must be valid FASTA format:
- Headers start with `>`
- Followed by sequence ID and optional description
- Sequences on subsequent lines

### 2. Aligned Sequences

**CRITICAL**: Sequences must be aligned before phylogenetic analysis!

- All sequences must be the same length
- Use alignment tools like:
  - [MUSCLE](https://www.drive5.com/muscle/)
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)
  - [Clustal Omega](http://www.clustal.org/omega/)
  - [MEGA](https://www.megasoftware.net/) (has built-in alignment)

### 3. Minimum Requirements

- At least 3 sequences
- Sequences must be homologous (related by evolution)
- DNA, RNA, or Protein sequences accepted

## Understanding Output

### Newick Format

Standard phylogenetic tree format:

```
((seq1:0.027,seq2:0.026):0.013,seq3:0.017):0.038,seq4:0.038);
```

Structure:
- Parentheses: Internal nodes (ancestors)
- Colons: Branch lengths (evolutionary distance)
- Semicolon: End of tree

With bootstrap support (if `--bootstrap` used):

```
((seq1:0.027,seq2:0.026)95:0.013,seq3:0.017)80:0.038,seq4:0.038);
```

Numbers before colons are bootstrap percentages (confidence values).

### Bootstrap Values

Bootstrap support indicates confidence in each branch:
- **95-100%**: Very strong support - branch is highly reliable
- **75-94%**: Moderate support - branch is reasonably confident
- **<75%**: Weak support - interpret with caution

## Troubleshooting

### Error: "Need at least 3 sequences"

You need minimum 3 sequences for phylogenetic analysis. Add more sequences.

### Error: "Sequences have different lengths"

Your sequences are not aligned. Use MUSCLE, MAFFT, or Clustal Omega to align them first.

### Error: "File not found"

Check that:
1. File path is correct
2. File exists
3. You're in the right directory

### Warning: "Mixed sequence types"

All sequences must be the same type (all DNA, all RNA, or all Protein). Check your FASTA file.

## Performance Tips

### For Short Sequences (<100 bp):

```bash
# Bootstrap will be fast
python rrna-phylo.py sequences.fasta --bootstrap 100
```

### For Medium Sequences (100-500 bp):

```bash
# ML may take a few seconds per replicate
python rrna-phylo.py sequences.fasta --bootstrap 50
```

### For Long Sequences (>500 bp):

```bash
# Use fewer replicates or distance methods
python rrna-phylo.py sequences.fasta --method bionj --bootstrap 100
```

### For Many Sequences (>20 taxa):

```bash
# ML tree search may take longer
python rrna-phylo.py sequences.fasta --method ml --verbose
```

## Getting Help

```bash
# Show help
python rrna-phylo.py --help

# Show version
python -c "import rrna_phylo; print(rrna_phylo.__version__)"
```

## Next Steps

After building your trees:

1. **Visualize** with FigTree or iTOL
2. **Interpret** bootstrap support values
3. **Compare** ML vs BioNJ vs UPGMA results
4. **Publish** using Newick files in your manuscript

## Advanced: Python API

If you need more control, use the Python API directly:

```python
from rrna_phylo import FastaParser, build_trees

# Parse sequences
parser = FastaParser()
sequences = parser.parse("sequences.fasta")

# Build trees
results = build_trees(sequences)

# Access trees
ml_tree, logL = results["ml"]
bionj_tree = results["bionj"]
upgma_tree = results["upgma"]
```

See [API-USAGE.md](API-USAGE.md) for complete API documentation.

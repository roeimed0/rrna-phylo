# rRNA-Phylo: Multi-Method Phylogenetic Tree Builder

High-performance phylogenetic tree inference implementing three complementary methods: distance-based (UPGMA, BioNJ) and maximum likelihood with automatic model selection.

**Dual Purpose Project:**
1. Production-ready phylogenetic analysis toolkit for rRNA sequences
2. Case study for optimizing Claude Code using skills and specialized agents

---

## What This Project Does

Builds phylogenetic trees from DNA/RNA sequences using three different algorithms:

- **UPGMA** - Fast, assumes molecular clock (constant evolution rate)
- **BioNJ** - Fast, variance-weighted, no clock assumption
- **Maximum Likelihood** - Rigorous statistical inference with automatic model selection

**Why three methods?** Compare results, validate with bootstrap, understand trade-offs between speed and accuracy.

---

## Quick Start

### Installation

```bash
# Create conda environment (recommended)
conda create -n rrna_phylo python=3.9
conda activate rrna_phylo

# Install dependencies
conda install numpy scipy numba biopython

# Navigate to project
cd backend
```

### Interactive Application (Easiest - Recommended)

**Step 1: Place FASTA files in data folder (optional)**
```bash
# Copy your FASTA files to the data folder for easy access
cp your_sequences.fasta backend/data/
```

**Step 2: Launch the interactive menu**
```bash
python app.py
```

The menu will display all files from `data/` folder for quick selection, or you can provide a custom path.

**Or use built-in commands:**
```bash
python app.py build sequences.fasta              # Build all 3 trees
python app.py build sequences.fasta --visualize pdf  # With visualization
python app.py test                                # Full test workflow
python app.py clean                              # Clean up test files
python app.py --help                             # Show all options
```

**Interactive menu includes:**
- ✅ **Data preparation** (deduplicate + clean headers)
- ✅ Quick build (one-click with defaults)
- ✅ Custom build (choose all options interactively)
- ✅ Advanced ML (model selection + tree search)
- ✅ Test data creation
- ✅ Results viewer
- ✅ Clean up tool
- ✅ ETE3 installation helper
- ✅ Built-in help

**No command-line flags to remember!**

### Data Preparation (Recommended First Step)

If you have raw FASTA files with duplicates or messy headers:

```bash
# Interactive menu (easiest)
python app.py
# Select Option 1: Prepare FASTA File

# Or direct CLI
python prepare_fasta.py data/raw_data.fasta data/clean_data.fasta
```

**What it does:**
1. Deduplicates sequences (keeps longest per species)
2. Cleans headers to standard `ID|Species_name` format
3. Example: `>AB571241.1|Obazoa|...|Homo_sapiens` → `>AB571241|Homo_sapiens`

**Result:** Clean, standardized FASTA ready for phylogenetic analysis with readable tree labels.

### Build All Three Tree Types (Direct CLI)

```bash
python rrna_phylo_cli.py sequences.fasta
```

**Output:** Creates `results/sequences/` with 7 files:
- 3 Newick trees (`.nwk` - standard format)
- 3 ASCII visualizations (`.txt` - readable text)
- 1 summary comparing all methods

### Build Single Method

```bash
# UPGMA (fastest, assumes clock)
python rrna_phylo_cli.py sequences.fasta --method upgma

# BioNJ (fast, no clock)
python rrna_phylo_cli.py sequences.fasta --method bionj

# Maximum Likelihood (slow, most accurate)
python rrna_phylo_cli.py sequences.fasta --method ml
```

### Add Bootstrap Support

```bash
# 100 replicates (recommended for publication)
python rrna_phylo_cli.py sequences.fasta --bootstrap 100

# 10 replicates (quick test)
python rrna_phylo_cli.py sequences.fasta --bootstrap 10
```

### Create Publication-Quality Visualizations

```bash
# PDF visualization (vector graphics, best for publications)
python rrna_phylo_cli.py sequences.fasta --visualize pdf

# High-resolution PNG (600 DPI for publication)
python rrna_phylo_cli.py sequences.fasta --visualize png --dpi 600

# SVG vector graphics (editable in Illustrator/Inkscape)
python rrna_phylo_cli.py sequences.fasta --visualize svg
```

**Requires:** `pip install ete3`

---

## Features

### Core Functionality

**Methods:**
- UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- BioNJ (Variance-weighted Neighbor-Joining)
- Maximum Likelihood with automatic model selection

**Model Selection (ML):**
- Tests 5 substitution models: JC69, K80, F81, HKY85, GTR
- Tests gamma rate variation (α = 0.5, 1.0, 2.0) for best models
- Selects best using BIC (Bayesian Information Criterion)
- Automatically applied - user doesn't need to know which model to use

**Tree Search (ML):**
- NNI (Nearest Neighbor Interchange) - fast, local topology changes
- SPR (Subtree Pruning and Regrafting) - thorough, larger rearrangements

**Bootstrap Analysis:**
- Felsenstein's bootstrap resampling (1985)
- Confidence values for each branch
- Available for all three methods

**Output Formats:**
- Newick format (standard, for downstream analysis)
- ASCII visualization (readable text tree)
- ETE3 visualization (PDF/PNG/SVG) - integrated in CLI with `--visualize` flag

### Performance Optimizations

**Implemented:**
- Site pattern compression (reduces redundant calculations)
- Numba JIT compilation (accelerates likelihood calculations)
- GPU acceleration (optional CUDA support)

**Improvements over baseline:**
- NumPy vectorization provides initial speedup
- Site compression reduces computation by ~6-10x on typical datasets
- Numba JIT adds additional acceleration for likelihood calculations
- Combined optimizations make ML inference practical for 50+ sequences

### Quality Features

- Automatic OpenMP conflict handling
- Organized output (results grouped by input filename)
- Comprehensive logging and timing
- Pre-aligned sequence detection
- Multiple output formats

---

## Testing

### Verify Installation

```bash
cd backend

# Create test dataset (50 sequences, 1000bp)
python rrna_phylo_test.py

# Build trees (takes ~6 seconds)
python rrna_phylo_cli.py test_sequences.fasta

# Verify output structure
python rrna_phylo_test.py --verify

# Clean up
bash cleanup_test.sh
```

**Expected output:**
```
[OK] Output directory exists: results\test_sequences
[OK] upgma_tree.nwk       (  1602 bytes)
[OK] bionj_tree.nwk       (  1589 bytes)
[OK] ml_tree.nwk          (  1589 bytes)
[OK] upgma_tree.txt       (  2742 bytes)
[OK] bionj_tree.txt       (  2729 bytes)
[OK] ml_tree.txt          (  2742 bytes)
[OK] summary.txt          (   683 bytes)
[OK] All expected files created successfully!
```

---

## Advanced Usage

### Maximum Likelihood with More Control

```bash
# Use NNI tree search (fast, default)
python rrna_phylo_ml.py sequences.fasta --method nni

# Use SPR tree search (thorough, slower)
python rrna_phylo_ml.py sequences.fasta --method spr

# With bootstrap
python rrna_phylo_ml.py sequences.fasta --bootstrap 100
```

**Output:**
- `results/sequences/tree.nwk` - ML tree
- `results/sequences/tree_ascii.txt` - ASCII visualization
- `results/sequences/metadata.json` - Model selection details, log-likelihood
- `results/sequences/build_log.txt` - Detailed build log

### Advanced Visualization Options

Using the `--visualize` flag, you can customize tree visualizations:

```bash
# Maximum likelihood with PDF output
python rrna_phylo_ml.py sequences.fasta --visualize pdf

# With bootstrap support shown
python rrna_phylo_ml.py sequences.fasta --bootstrap 100 --visualize pdf

# Ultra high-resolution for publication (1200 DPI)
python rrna_phylo_ml.py sequences.fasta --visualize png --dpi 1200
```

**Supported formats:** PDF, PNG, SVG, EPS

**Requires:** `pip install ete3`

---

## Output Files

### Directory Structure

```
results/
└── [filename]/              # Named after input (without extension)
    ├── upgma_tree.nwk       # UPGMA tree (Newick)
    ├── bionj_tree.nwk       # BioNJ tree (Newick)
    ├── ml_tree.nwk          # ML tree (Newick)
    ├── upgma_tree.txt       # UPGMA tree (ASCII)
    ├── bionj_tree.txt       # BioNJ tree (ASCII)
    ├── ml_tree.txt          # ML tree (ASCII)
    └── summary.txt          # Comparison of methods
```

**Example:**
```bash
python rrna_phylo_cli.py my_data.fasta
# Creates: results/my_data/
```

### File Formats

**Newick (`.nwk`)** - Standard format for phylogenetic trees:
```
(Species_A:0.01,Species_B:0.02,(Species_C:0.03,Species_D:0.04):0.05);
```

**ASCII (`.txt`)** - Human-readable tree visualization:
```
Species_A         |--------
                  |
Species_B         |-------------
                  |
                  |-------- Species_C
                  |
                  `-------- Species_D
```

---

## Performance Benchmarks

| Dataset | Method | Time | Notes |
|---------|--------|------|-------|
| 10 seqs, 500bp | UPGMA | <1s | Instant |
| 10 seqs, 500bp | BioNJ | <1s | Instant |
| 10 seqs, 500bp | ML | ~2s | Model selection + NNI |
| 50 seqs, 1000bp | UPGMA | <1s | Instant |
| 50 seqs, 1000bp | BioNJ | <1s | Instant |
| 50 seqs, 1000bp | ML | ~6s | Model selection + NNI |
| 50 seqs, 1000bp | ML + bootstrap 100 | ~7min | Parallelizable |

**Notes:**
- Distance-based methods (UPGMA, BioNJ) are instant for datasets up to 100 sequences
- ML is practical for 50-100 sequences with NNI search
- SPR search is slower but more thorough
- Bootstrap performance scales linearly with replicates

---

## Project Structure

```
rrna-phylo/
├── README.md                          # This file
├── backend/
│   ├── app.py                         # Main entry point (CLI & menu)
│   ├── data/                          # Place your FASTA files here
│   │   └── README.md                  # Data folder documentation
│   ├── rrna_phylo/                    # Core package
│   │   ├── io/                        # FASTA parsing
│   │   │   └── fasta_parser.py
│   │   ├── core/                      # Tree structures
│   │   │   ├── sequence.py
│   │   │   ├── tree.py
│   │   │   └── builder.py
│   │   ├── distance/                  # Distance calculations
│   │   │   └── distance.py            # Jukes-Cantor, Kimura
│   │   ├── methods/                   # Distance-based methods
│   │   │   ├── upgma.py
│   │   │   ├── nj.py
│   │   │   └── bionj.py
│   │   ├── models/                    # Maximum Likelihood
│   │   │   ├── ml_tree_level3.py      # ML with site compression
│   │   │   ├── ml_tree_level4.py      # ML with model selection
│   │   │   ├── model_selection.py     # BIC/AIC model comparison
│   │   │   ├── rate_matrices.py       # GTR, HKY, K80, JC69, F81
│   │   │   ├── nni_search.py          # NNI tree search
│   │   │   ├── spr_search.py          # SPR tree search
│   │   │   └── bootstrap.py           # Bootstrap analysis
│   │   └── visualization/             # Tree visualization
│   │       ├── ascii_viz.py           # ASCII trees
│   │       └── ete3_viz.py            # PDF/PNG/SVG visualization
│   ├── rrna_phylo_app.py              # Interactive menu
│   ├── rrna_phylo_cli.py              # Main CLI (all 3 methods)
│   ├── rrna_phylo_ml.py               # Advanced ML CLI
│   ├── rrna_phylo_test.py             # Test suite
│   ├── prepare_fasta.py               # Data preparation (deduplicate + clean)
│   ├── cleanup_test.sh                # Cleanup script
│   ├── diagnose_openmp.py             # OpenMP diagnostics
│   └── OPENMP_CONFLICT_EXPLANATION.md # OpenMP troubleshooting
```

---

## Troubleshooting

### OpenMP Library Conflict

**Symptom:**
```
OMP: Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already initialized.
```

**Solution:** Already handled automatically in both CLI scripts.

If you still see this:
```bash
export KMP_DUPLICATE_LIB_OK=TRUE  # Linux/Mac
set KMP_DUPLICATE_LIB_OK=TRUE     # Windows
```

See [OPENMP_CONFLICT_EXPLANATION.md](backend/OPENMP_CONFLICT_EXPLANATION.md) for details.

### Common Issues

**"ModuleNotFoundError: No module named 'rrna_phylo'"**
- Make sure you're in the `backend/` directory

**Trees look identical**
- Normal for very similar sequences
- Try more diverse data

**ML is slow**
- Use NNI instead of SPR: `--method nni`
- Reduce sequences or alignment length
- Skip bootstrap for testing

---

## Claude Code Optimization Case Study

This project demonstrates **systematic optimization of LLM-assisted development** through skills and specialized agents.

### 1. Skills-Based Knowledge

Custom skills created:
- `phylogenetic-methods` - Comprehensive phylogenetics algorithms
- `ml-tree-methods` - Maximum likelihood patterns
- `ml-tree-level4` - Advanced ML with model selection
- `consensus-tree-methods` - Tree comparison
- `project-architecture-patterns` - FastAPI backend design
- `tree-visualization` - Matplotlib/ggtree/R patterns

**Benefit:** Domain expertise available on-demand without re-explaining concepts

### 2. Agent-Based Workflow

Specialized agents:
- **Explore agent** - Fast codebase navigation and search
- **Plan agent** - Implementation strategy before coding
- **Architecture review agent** - Code quality assessment
- **Refactor planner** - Systematic refactoring plans

**Benefit:** Right tool for the job - exploration vs. implementation vs. review

### 3. Systematic Development Process

**Phases:**
1. **Research** - Understanding phylogenetic algorithms and existing tools
2. **Implementation** - Building from simple (UPGMA) to complex (ML with model selection)
3. **Optimization** - Profiling bottlenecks, applying NumPy/Numba/site compression
4. **Testing** - Verification with real datasets
5. **Documentation** - Comprehensive usage and troubleshooting guides

**Key lesson:** Benchmark before deleting code. We almost removed the faster implementation during cleanup!

### 4. Performance Evolution

**Implementation levels:**
1. **Pure Python** - Baseline, easy to understand, slow
2. **NumPy vectorization** - Matrix operations, faster
3. **Site pattern compression** - Reduce redundant calculations (6-10x improvement on typical data)
4. **Numba JIT** - Compile hot loops, additional acceleration
5. **GPU support** - CUDA for parallel likelihood (optional)

**Result:** ML inference practical for 50+ sequences

### 5. Production-Ready Features

- Automatic model selection (no phylogenetics expertise required)
- Organized outputs (`results/[filename]/`)
- OpenMP conflict handling (automatic)
- Multiple output formats
- Comprehensive error messages
- Performance logging

---

## Implementation Details

### Maximum Likelihood Inference

**Model Selection:**
1. Build initial tree with BioNJ (fast, no assumptions)
2. Test 5 models: JC69, K80, F81, HKY85, GTR
3. Test gamma rate variation for promising models
4. Select best using BIC
5. Use best model for final tree

**Likelihood Calculation:**
- Felsenstein's pruning algorithm
- Site pattern compression
- Numba JIT compilation
- Optional GPU acceleration

**Tree Search:**
- **NNI** - Swaps adjacent branches, fast
- **SPR** - Moves entire subtrees, thorough

### Distance-Based Methods

**UPGMA:**
- O(n²) time complexity
- Assumes molecular clock
- Good for closely related sequences

**BioNJ:**
- Variance-weighted distances
- No clock assumption
- More accurate than standard Neighbor-Joining

### Bootstrap Analysis

- Resamples alignment columns with replacement
- Rebuilds tree for each replicate
- Reports % of replicates supporting each branch
- Standard method since Felsenstein (1985)

---

## Requirements

- Python 3.8+
- NumPy 2.0+
- SciPy 1.13+
- Numba 0.60+
- BioPython 1.80+

**Optional:**
- ETE3 (for PDF/PNG/SVG visualization)
- CUDA (for GPU acceleration)

---

## Future Enhancements

### High Priority
- [x] Integrate ETE3 visualization into CLI (PDF/PNG/SVG output) - DONE
- [ ] Parallel bootstrap (joblib/multiprocessing)

### Medium Priority
- [ ] Consensus tree methods (majority-rule, strict)
- [ ] Tree comparison metrics (Robinson-Foulds distance)
- [ ] Progress bars (tqdm)
- [ ] FastAPI web interface

### Low Priority
- [ ] Interactive visualization (web-based)
- [ ] Distributed computing support
- [ ] ML integration (rRNA classification, tree synthesis)

---

## References

### Phylogenetics
- Felsenstein, J. (1985). Confidence limits on phylogenies: An approach using the bootstrap. *Evolution*, 39(4), 783-791.
- Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm. *Molecular Biology and Evolution*, 14(7), 685-695.
- Stamatakis, A. (2014). RAxML version 8. *Bioinformatics*, 30(9), 1312-1313.

### Model Selection
- Posada, D., & Crandall, K. A. (1998). MODELTEST: testing the model of DNA substitution. *Bioinformatics*, 14(9), 817-818.

### Optimization
- Numba Documentation: https://numba.pydata.org/
- Site Pattern Compression: RAxML implementation notes

---

## License

MIT License

## Acknowledgments

Built with Claude Code (Sonnet 4.5) as a case study for optimizing LLM-assisted development through skills and specialized agents.

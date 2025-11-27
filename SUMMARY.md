# Recent Changes - rRNA-Phylo

**Last Updated**: 2025-11-26

---

## Session 2025-11-26: Display Names & Aligned Visualization

### Completed ✅

#### 1. Display Names in Trees
**Format**: `"Species name (MainAccession) #N"`

**Examples**:
- `Escherichia coli str. K-12 substr. MG1655 (U00096) #1`
- `Salmonella virus Fels2 (AE006468) #1`
- `Pseudomonas aeruginosa PAO1 (AE004091) #1`

**Benefits**:
- Readable species names (not cryptic accessions)
- Simplified accessions (U00096 vs U00096.223771.225312)
- Numbered duplicates for clarity (#1, #2, etc.)

**Implementation**:
- Added `species_name` property to extract from FASTA descriptions
- Added `main_accession` property to simplify accessions
- Added `display_name` property combining both
- Added `assign_unique_display_names()` for #N numbering
- Applied to all tree output formats (ASCII, Newick, ETE3)

#### 2. ETE3 Aligned Leaf Names
**Feature**: All species names aligned in vertical column

**Settings**:
- `force_topology=True` - Align all leaf names
- `draw_guiding_lines=True` - Black solid guide lines
- `guiding_lines_type=0` - Solid lines (no dots)
- `guiding_lines_color="black"` - Professional appearance

**Result**: Publication-quality trees with clean, aligned labels

**Fixed Issues**:
- Newick semicolon missing (caused ETE3 parsing errors)
- Double-quoted names (removed quotes from display_name)
- Added `quoted_node_names=True` to ETE3 Tree constructor

#### 3. Comprehensive Test Pipeline
**File**: `test_full_pipeline.py`

**What it does**:
1. Load and align sequences
2. Remove exact duplicates (always)
3. Build trees in TWO modes:
   - Regular (exact duplicates only) - 19 sequences
   - Deduplicated (smart deduplication) - 8 sequences
4. For each mode, build THREE tree types:
   - UPGMA (distance-based)
   - BioNJ (variance-weighted)
   - ML (Maximum Likelihood with GTR+Gamma)
5. Generate ONE COMBINED OUTPUT FILE showing all 6 trees
6. Generate individual files for each tree

**Total Output**: 26 files
- 1 combined comparison file (COMPLETE_TREE_COMPARISON.txt)
- 6 ASCII files (individual tree visualizations)
- 6 Newick files (standard format)
- 6 PDF files (ETE3 publication-quality)
- 6 PNG files (ETE3 high-resolution)
- 1 aligned FASTA file

### Files Modified
- [fasta_parser.py](backend/rrna_phylo/io/fasta_parser.py) - display_name property
- [ete3_viz.py](backend/rrna_phylo/visualization/ete3_viz.py) - aligned leaf names, guide lines
- [test_full_pipeline.py](backend/tests/test_full_pipeline.py) - comprehensive pipeline
- [tree.py](backend/rrna_phylo/core/tree.py) - Newick semicolon fix

### Understanding Deduplication Impact

**Regular mode** (19 sequences):
- Shows all unique sequences after exact duplicate removal
- E. coli: 7 rRNA operons (slightly different)
- Salmonella: 7 rRNA operons
- Good for studying intra-strain variation

**Deduplicated mode** (8 sequences):
- Smart deduplication at 99.5% similarity
- E. coli: 3 representatives (showing meaningful diversity)
- Salmonella: 2 representatives
- Good for inter-species phylogenetic relationships

**Key Insight**: Deduplicated trees keep MULTIPLE representatives per species when they show meaningful variation (not just 1 per species)

---

## Session 2025-11-25: Smart Deduplication & Testing

### Completed ✅

#### 1. Smart Deduplication
**Two-tier strategy**:
- **Tier 1**: Exact duplicate removal (always applied)
- **Tier 2**: Similarity clustering at 99.5% (optional with `--dereplicate`)

**Species-aware clustering**:
- Only clusters sequences from the SAME genome/species
- Prevents accidental merging of phylogenetically distinct sequences
- Critical for accurate phylogenetic inference

**Results with test data**:
- 24 sequences → 19 after exact duplicate removal
- 24 sequences → 11 after smart deduplication (54% reduction)
- Phylogenetic relationships preserved

#### 2. Comprehensive Test Coverage
**File**: `test_strain_handler.py`

**28 tests covering**:
- Exact duplicate removal (5 tests)
- Sequence similarity calculation (4 tests)
- Similarity-based clustering (4 tests)
- Smart deduplication pipeline (4 tests)
- Representative selection (5 tests)
- Strain grouping (2 tests)
- Legacy functions (1 test)
- Strain summary (1 test)
- Full integration (2 tests)

**All 28 tests PASSED** ✅

### Files Created
- [strain_handler.py](backend/rrna_phylo/utils/strain_handler.py) - All deduplication logic
- [test_strain_handler.py](backend/tests/test_strain_handler.py) - 28 comprehensive tests

---

## Session 2025-11-21: Bootstrap Analysis

### Completed ✅

#### Bootstrap Support Values
**Implementation**:
- 100-1000 replicates supported
- Parallel processing using all CPU cores
- Works with UPGMA, BioNJ, and ML methods
- Reproducible with controlled random seeds

**Performance**:
- UPGMA: ~2 seconds (100 replicates, 19 sequences)
- BioNJ: ~3 seconds (100 replicates, 19 sequences)
- ML: ~120 seconds (100 replicates, 19 sequences with Numba)

**Testing**:
- Verified with 2-200 replicates
- All three methods tested
- Bootstrap values in valid range (0-100%)

### Files Modified
- [builder.py](backend/rrna_phylo/core/builder.py) - Bootstrap integration
- [bootstrap.py](backend/rrna_phylo/consensus/bootstrap.py) - Bootstrap implementation

---

## Earlier Work (November 2025)

### ETE3 Visualization (2025-11-25)
- Publication-quality tree output (Nature/Science/Cell standard)
- Multiple formats: PDF, PNG, SVG, EPS
- Bootstrap support display with color coding
- CLI integration with 7 new flags

### ML Level 4 (2025-11-20)
- Automatic model selection (AIC/BIC)
- NNI tree search
- Numba acceleration (9x speedup)
- Model comparison utilities

### Performance Optimization (2025-11-19)
- Site pattern compression (10-100x speedup)
- Numba JIT compilation (9x speedup for likelihood)
- Fast matrix exponential (Pade approximation)
- Combined: ~300x vs naive implementation

---

## Current Status

### Production-Ready Features ✅
1. **Phylogenetic Tree Building**
   - UPGMA (distance-based, molecular clock)
   - BioNJ (variance-weighted neighbor-joining)
   - ML (Maximum Likelihood with GTR+Gamma)

2. **Bootstrap Analysis**
   - 100-1000 replicates
   - Parallel processing (all CPU cores)
   - Integrated with all methods

3. **ETE3 Visualization**
   - Publication-quality output
   - Aligned leaf names
   - Multiple formats (PDF, PNG, SVG, EPS)
   - Bootstrap display

4. **Smart Deduplication**
   - Two-tier strategy (exact + similarity)
   - Species-aware clustering
   - 40-60% dataset reduction
   - 28 comprehensive tests

5. **Display Names**
   - Human-readable species names
   - Simplified accessions
   - Numbered duplicates

### Known Issues ⚠️
1. **Consensus Trees** - BROKEN (disabled since 2025-11-20)
   - Algorithm creates incorrect bipartitions
   - Needs complete redesign
   - Not blocking other work

### Next Priorities

**Immediate** (This Week):
1. Fix Consensus Trees - Proper algorithm (3-5 days)
2. Progress Bars - Add tqdm (2-3 hours)

**Short Term** (Next 2-4 Weeks):
3. **Generative AI Phase 1** - ML-Enhanced Quality Scoring (1 week)
4. **Generative AI Phase 2** - Multi-Tree Ensemble Methods (1-2 weeks)
5. Tree Comparison CLI - RF distance (1 day)
6. Parallel Bootstrap - joblib for Windows (2-3 hours)

**Medium Term** (Next 1-2 Months):
7. **Generative AI Phase 3** - GNN/Transformer Tree Generation (2-3 weeks)
8. SPR Tree Search - Better than NNI (6-8 hours)
9. Advanced Visualization - Branch coloring, clade highlighting (1-2 days)
10. Performance Profiling - Optimize bottlenecks (2-3 days)

---

## Documentation Files

### User Documentation
- [USAGE_GUIDE.md](backend/docs/user-guide/usage-guide.md) - Complete user guide
- [DEDUPLICATION_SUMMARY.md](backend/DEDUPLICATION_SUMMARY.md) - Deduplication details
- [visualization.md](backend/docs/features/visualization.md) - ETE3 visualization

### Developer Documentation
- [STATUS.md](STATUS.md) - Complete project status (UPDATED 2025-11-26)
- [ROADMAP.md](ROADMAP.md) - Development roadmap
- [DOCUMENTATION_REVIEW.md](DOCUMENTATION_REVIEW.md) - Doc improvement plan (NEW)
- [architecture.md](backend/docs/development/architecture.md) - Project architecture
- [performance.md](backend/docs/development/performance.md) - Performance optimization

### Skills
9 Claude Code skills available for development guidance:
- skill-developer
- rrna-prediction-patterns
- phylogenetic-methods
- project-architecture-patterns
- ml-integration-patterns
- ml-tree-methods
- consensus-tree-methods
- ml-tree-level4
- tree-visualization

---

## Quick Start

### Run Comprehensive Test
```bash
cd backend
python tests/test_full_pipeline.py
```

**Output**: 26 files showing all 6 trees (regular + deduplicated × UPGMA + BioNJ + ML)

### Build a Tree
```bash
cd backend
python -m rrna_phylo.cli test_real_rrana.fasta \
    --method ml \
    --bootstrap 100 \
    --dereplicate \
    --visualize \
    --output-format both \
    -o results/
```

**Output**:
- `results/tree_ml.nwk` - Newick format
- `results/tree_ml_ascii.txt` - ASCII visualization
- `results/tree_ml_tree.pdf` - Publication-quality PDF

---

**Generated**: 2025-11-26
**Status**: Documentation updated to reflect all recent work

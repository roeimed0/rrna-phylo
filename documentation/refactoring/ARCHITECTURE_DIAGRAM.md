# rRNA-Phylo Architecture - Before & After Refactoring

---

## Current Architecture (Before)

```
backend/rrna_phylo/
│
├── core/                       [4 files, ~1,000 lines] ✅ GOOD
│   ├── tree.py                 # TreeNode class
│   ├── sequence_type.py        # DNA/RNA/Protein detection
│   ├── builder.py              # Main tree builder
│   └── builder_smart.py        # Smart builder (dataset analysis)
│
├── io/                         [3 files, ~400 lines] ✅ GOOD
│   ├── fasta_parser.py         # FASTA parsing
│   ├── aligner.py              # MUSCLE wrapper
│   └── __init__.py
│
├── distance/                   [3 files, ~400 lines] ✅ GOOD
│   ├── distance.py             # Jukes-Cantor for DNA/RNA
│   ├── protein_distance.py     # Protein distances
│   └── __init__.py
│
├── methods/                    [4 files, ~900 lines] ✅ GOOD
│   ├── upgma.py                # UPGMA algorithm
│   ├── bionj.py                # BioNJ algorithm
│   ├── protein_ml.py           # Protein ML
│   └── __init__.py
│
├── models/                     [13 files, ~5,000 lines] ⚠️ LARGE
│   ├── ml_tree.py              # Level 1 (basic GTR)
│   ├── ml_tree_level2.py       # Level 2 (Felsenstein)
│   ├── ml_tree_level3.py       # Level 3 (production)
│   ├── ml_tree_level4.py       # Level 4 (full features)
│   ├── model_selection.py      # AIC/BIC selection
│   ├── substitution_models.py  # DNA models (JC69, K80, HKY, GTR)
│   ├── rate_matrices.py        # Rate matrix calculations
│   ├── branch_length_optimizer.py  # Branch optimization
│   ├── tree_search.py          # NNI topology search
│   ├── numba_likelihood.py     # Numba acceleration
│   ├── gpu_likelihood_torch.py # GPU acceleration (PyTorch)
│   ├── protein_models.py       # WAG, LG, JTT models
│   └── __init__.py
│
├── utils/                      [7 files, ~2,500 lines] 🔴 PROBLEMATIC
│   ├── strain_handler.py       # 🔴 NOT UTILITY - Deduplication logic
│   ├── bootstrap.py            # 🔴 NOT UTILITY - Statistical analysis
│   ├── dataset_analyzer.py     # 🔴 NOT UTILITY - Dataset analysis
│   ├── sampling_strategy.py    # 🔴 NOT UTILITY - Sampling strategies
│   ├── outgroup_handler.py     # 🔴 NOT UTILITY - Outgroup selection
│   ├── visualize_trees.py      # 🔴 NOT UTILITY - Tree visualization
│   └── console.py              # ✅ ACTUAL UTILITY - Console formatting
│
├── consensus/                  [3 files, ~400 lines] ✅ GOOD
│   ├── bipartitions.py         # Bipartition extraction
│   ├── tree_distance.py        # Robinson-Foulds distance
│   └── __init__.py
│
├── visualization/              [2 files, ~100 lines] ✅ GOOD
│   ├── ete3_viz.py             # ETE3 publication-quality viz
│   └── __init__.py
│
├── cli.py                      [686 lines] ✅ ACCEPTABLE
└── config.py                   [~30 lines] ⚠️ MINIMAL

Total: 44 Python files, ~11,400 lines
```

### Issues Highlighted

```
❌ Problem: 6 files (~2,200 lines) misplaced in utils/
❌ Problem: utils/ is a catch-all for non-utility code
⚠️  Issue: models/ is large (13 files, 5,000 lines) but manageable
⚠️  Issue: ML level progression unclear
⚠️  Issue: Two builder classes with overlapping functionality
```

---

## Proposed Architecture (After Phase 1)

```
backend/rrna_phylo/
│
├── core/                       [4 files, ~1,000 lines] ✅
│   ├── tree.py
│   ├── sequence_type.py
│   ├── builder.py
│   └── builder_smart.py        # TODO: Document or merge
│
├── io/                         [4 files, ~450 lines] ✅
│   ├── fasta_parser.py
│   ├── aligner.py
│   ├── console.py              # 🟢 MOVED from utils/
│   └── __init__.py
│
├── distance/                   [3 files, ~400 lines] ✅
│   ├── distance.py
│   ├── protein_distance.py
│   └── __init__.py
│
├── methods/                    [4 files, ~900 lines] ✅
│   ├── upgma.py
│   ├── bionj.py
│   ├── protein_ml.py
│   └── __init__.py
│
├── models/                     [13 files, ~5,000 lines] ✅ BETTER ORGANIZED
│   ├── ml_tree_level4.py       # 🟢 Main: Use this for production
│   ├── gpu_likelihood.py       # 🟢 RENAMED (was gpu_likelihood_torch.py)
│   ├── model_selection.py
│   ├── substitution_models.py
│   ├── rate_matrices.py
│   ├── branch_length_optimizer.py
│   ├── tree_search.py
│   ├── numba_likelihood.py
│   ├── protein_models.py
│   ├── reference/              # 🟢 NEW: Reference implementations
│   │   ├── README.md           # Explains ML progression
│   │   ├── ml_tree_level1.py   # Educational: Basic GTR
│   │   ├── ml_tree_level2.py   # Educational: Felsenstein
│   │   └── ml_tree_level3.py   # Reference: Production without model selection
│   └── __init__.py
│
├── preprocessing/              # 🟢 NEW PACKAGE
│   ├── __init__.py             # [4 files, ~1,100 lines]
│   ├── strain_handler.py       # 🟢 MOVED from utils/ (512 lines)
│   ├── outgroup_handler.py     # 🟢 MOVED from utils/ (~150 lines)
│   ├── sampling_strategy.py    # 🟢 MOVED from utils/ (350 lines)
│   └── README.md               # Purpose and usage
│
├── analysis/                   # 🟢 NEW PACKAGE
│   ├── __init__.py             # [3 files, ~700 lines]
│   ├── bootstrap.py            # 🟢 MOVED from utils/ (398 lines)
│   ├── dataset_analyzer.py     # 🟢 MOVED from utils/ (294 lines)
│   └── README.md               # Purpose and usage
│
├── consensus/                  [3 files, ~400 lines] ✅
│   ├── bipartitions.py
│   ├── tree_distance.py
│   └── __init__.py
│
├── visualization/              [3 files, ~250 lines] ✅
│   ├── ete3_viz.py
│   ├── ascii_viz.py            # 🟢 MOVED from utils/visualize_trees.py
│   └── __init__.py
│
├── cli.py                      [686 lines] ✅
└── config.py                   [~30 lines]

Total: 48 Python files (~11,400 lines + new READMEs)
```

### Improvements

```
✅ Fixed: No non-utility modules in utils/ (removed entirely or moved console to io/)
✅ Fixed: Clear preprocessing package for data preparation
✅ Fixed: Clear analysis package for statistical tools
✅ Fixed: Visualization consolidated in one package
✅ Fixed: ML reference implementations archived with documentation
✅ Fixed: Implementation detail removed from filename (gpu)
✅ Better: Package structure matches functional boundaries
```

---

## Package Dependency Graph

### Before (Current)

```
                    ┌──────────────┐
                    │     core/    │ (Foundation)
                    │  TreeNode    │
                    │  Sequence    │
                    └──────┬───────┘
                           │
            ┌──────────────┼──────────────┐
            │              │              │
       ┌────▼───┐    ┌────▼────┐   ┌────▼────┐
       │  io/   │    │distance/│   │methods/ │
       │ FASTA  │    │ Jukes-  │   │ UPGMA   │
       │ MUSCLE │    │ Cantor  │   │ BioNJ   │
       └────┬───┘    └────┬────┘   └────┬────┘
            │             │              │
            └─────────┬───┴──────────────┘
                      │
                 ┌────▼─────┐
                 │  models/ │ (ML, GPU, optimization)
                 │   5000+  │
                 │   lines  │
                 └────┬─────┘
                      │
        ┌─────────────┼─────────────┐
        │             │             │
   ┌────▼────┐  ┌────▼─────┐  ┌───▼──────┐
   │ utils/  │  │consensus/│  │visualiza-│
   │  MESS!  │  │  trees   │  │  tion/   │
   │  🔴     │  │          │  │          │
   └─────────┘  └──────────┘  └──────────┘
```

### After (Proposed)

```
                    ┌──────────────┐
                    │     core/    │ (Foundation)
                    │  TreeNode    │
                    │  Sequence    │
                    └──────┬───────┘
                           │
            ┌──────────────┼──────────────┐
            │              │              │
       ┌────▼───┐    ┌────▼────┐   ┌────▼────┐
       │  io/   │    │distance/│   │methods/ │
       │ FASTA  │    │ Jukes-  │   │ UPGMA   │
       │ MUSCLE │    │ Cantor  │   │ BioNJ   │
       └────────┘    └─────────┘   └────┬────┘
                                         │
                 ┌───────────────────────┴────────┐
                 │                                │
            ┌────▼─────┐                    ┌────▼─────┐
            │  models/ │                    │analysis/ │
            │    ML    │                    │bootstrap │
            │   GPU    │                    │ dataset  │
            │optimize  │                    │ analyzer │
            └────┬─────┘                    └──────────┘
                 │
     ┌───────────┼───────────┐
     │           │           │
┌────▼──────┐ ┌─▼────────┐ ┌▼────────┐
│preprocess-│ │consensus/│ │visualiza-│
│   ing/    │ │  trees   │ │  tion/   │
│ 🟢 CLEAR  │ │          │ │ ascii+   │
│  PURPOSE  │ │          │ │ ete3     │
└───────────┘ └──────────┘ └─────────┘

Legend:
  🔴 = Problem area
  🟢 = Fixed/improved
  ✅ = Already good
```

---

## Data Flow - Tree Building Pipeline

### Full Pipeline (All Three Methods)

```
┌─────────────────┐
│  Input FASTA    │
│  sequences.fasta│
└────────┬────────┘
         │
         ▼
┌─────────────────────────┐
│  preprocessing/         │  🟢 NEW PACKAGE
│  - Deduplication        │
│  - Strain handling      │
│  - Sampling (if large)  │
└────────┬────────────────┘
         │
         ▼
┌─────────────────────────┐
│  io/aligner.py          │
│  MUSCLE alignment       │
└────────┬────────────────┘
         │
         ▼
┌─────────────────────────┐
│  core/sequence_type.py  │
│  Detect DNA/RNA/Protein │
└────────┬────────────────┘
         │
         ▼
┌─────────────────────────┐
│  distance/distance.py   │
│  Calculate distance     │
│  matrix (Jukes-Cantor)  │
└────────┬────────────────┘
         │
         ├──────┬──────────────┬──────────┐
         │      │              │          │
         ▼      ▼              ▼          ▼
   ┌─────────┐ ┌──────────┐ ┌──────────────────┐
   │ UPGMA   │ │  BioNJ   │ │  ML (GTR+Gamma)  │
   │ methods/│ │ methods/ │ │  models/ml_tree_ │
   │ upgma.py│ │ bionj.py │ │  level4.py       │
   └────┬────┘ └────┬─────┘ └────┬─────────────┘
        │           │            │
        │           │            ▼
        │           │      ┌─────────────────┐
        │           │      │ models/tree_    │
        │           │      │ search.py (NNI) │
        │           │      └────┬────────────┘
        │           │           │
        └───────────┴───────────┘
                    │
                    ▼
         ┌──────────────────────┐
         │ analysis/bootstrap.py│  🟢 NEW PACKAGE
         │ (Optional)            │
         │ n=100 replicates     │
         └──────────┬───────────┘
                    │
                    ▼
         ┌──────────────────────┐
         │ consensus/           │
         │ Compare trees        │
         │ Robinson-Foulds dist │
         └──────────┬───────────┘
                    │
                    ▼
         ┌──────────────────────┐
         │ visualization/       │  🟢 IMPROVED
         │ - ASCII (terminal)   │
         │ - ETE3 (publication) │
         └──────────────────────┘
```

---

## Import Path Changes

### Phase 1 Changes

```python
# OLD (Before)                          # NEW (After)
from rrna_phylo.utils.strain_handler   →  from rrna_phylo.preprocessing.strain_handler
from rrna_phylo.utils.bootstrap        →  from rrna_phylo.analysis.bootstrap
from rrna_phylo.utils.dataset_analyzer →  from rrna_phylo.analysis.dataset_analyzer
from rrna_phylo.utils.visualize_trees  →  from rrna_phylo.visualization.ascii_viz
from rrna_phylo.utils.console          →  from rrna_phylo.io.console

# GPU module rename
from rrna_phylo.models.gpu_likelihood_torch  →  from rrna_phylo.models.gpu_likelihood
```

### Backward Compatibility (Transition Period)

```python
# rrna_phylo/__init__.py
# Maintain old imports with deprecation warnings for 1-2 versions

from rrna_phylo.preprocessing import strain_handler as _sh
import warnings

# For backward compatibility (deprecated)
def __getattr__(name):
    if name == "strain_handler":
        warnings.warn("Import from rrna_phylo.preprocessing instead", DeprecationWarning)
        return _sh
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
```

---

## CLI Workflow Mapping

### CLI Command Flow (Example: Full Pipeline)

```bash
python -m rrna_phylo.cli sequences.fasta \
    --method all \
    --bootstrap 100 \
    --dereplicate \
    --visualize \
    -o results/
```

**Execution Path**:

```
cli.py (main entry point)
  │
  ├─► preprocessing/strain_handler.py (--dereplicate)
  │
  ├─► io/fasta_parser.py (read sequences)
  │
  ├─► io/aligner.py (align if needed)
  │
  ├─► core/builder.py or core/builder_smart.py
  │     │
  │     ├─► methods/upgma.py (build UPGMA tree)
  │     ├─► methods/bionj.py (build BioNJ tree)
  │     └─► models/ml_tree_level4.py (build ML tree)
  │           │
  │           ├─► models/model_selection.py (select best model)
  │           └─► models/tree_search.py (NNI optimization)
  │
  ├─► analysis/bootstrap.py (--bootstrap 100)
  │
  ├─► consensus/tree_distance.py (compare trees)
  │
  └─► visualization/ete3_viz.py (--visualize)
        or
      visualization/ascii_viz.py (terminal output)
```

---

## File Size Distribution

### Before

```
Largest files (>400 lines):
  686  cli.py                        ✅ Acceptable (many CLI options)
  593  models/ml_tree_level3.py      ⚠️ Large but complex ML
  538  models/gpu_likelihood_torch.py⚠️ GPU implementation
  512  utils/strain_handler.py       🔴 Misplaced, should move
  491  models/model_selection.py     ⚠️ Complex model selection
  451  core/builder.py               ⚠️ Could extract alignment
  443  models/ml_tree_level2.py      ⚠️ Reference implementation
  442  models/substitution_models.py ✅ Reasonable (many models)
  429  methods/protein_ml.py         ✅ Reasonable (protein ML)
  420  models/tree_search.py         ⚠️ NNI implementation
  409  models/ml_tree_level4.py      ✅ Main ML interface
```

### After (Potential Splits)

```
If we extract from large files:

builder.py (451 lines)
  → builder.py (350 lines) + alignment_utils.py (100 lines)

ml_tree_level3.py (593 lines)
  → ml_tree_level3.py (450 lines) + gamma_rates.py (150 lines)

But: Current sizes are acceptable for complex functionality
Recommendation: Keep current sizes, improve documentation instead
```

---

## Testing Impact

### Tests to Update After Phase 1

```python
# Test files that will need import updates:

tests/test_strain_handler.py
  - Update: from rrna_phylo.preprocessing.strain_handler import ...

tests/test_bootstrap.py
  - Update: from rrna_phylo.analysis.bootstrap import ...

tests/test_dataset_analyzer.py
  - Update: from rrna_phylo.analysis.dataset_analyzer import ...

tests/test_visualization.py
  - Update: from rrna_phylo.visualization.ascii_viz import ...

tests/integration/test_full_pipeline.py
  - Update: Multiple imports for moved modules
```

### Test Coverage Goals

```
Target coverage after refactoring:
  core/          >90%  (critical functionality)
  methods/       >85%  (algorithms)
  models/        >80%  (ML is complex, focus on critical paths)
  preprocessing/ >80%  (important for data quality)
  analysis/      >85%  (statistical correctness)
  io/            >75%  (I/O operations)
  visualization/ >60%  (visual output, harder to test)

Overall target: >80%
```

---

## Summary

### Before Refactoring
- 44 files in 8 packages
- utils/ contains non-utility code (2,500 lines misplaced)
- ML progression unclear
- Some confusion about which modules to use

### After Refactoring (Phase 1)
- 48 files in 9 packages (added preprocessing/, analysis/)
- Clear package purposes
- Better discoverability
- Maintained backward compatibility
- All tests pass

### Key Improvements
1. ✅ Logical package organization
2. ✅ Clear functional boundaries
3. ✅ Better discoverability for new users
4. ✅ Maintained all functionality
5. ✅ Low-risk changes (mostly file moves)

---

**Next**: Implement Phase 1 (package reorganization) first, then assess if further changes are needed.

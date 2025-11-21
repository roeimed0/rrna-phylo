# Backend Refactoring Plan - rRNA-Phylo

**Date**: 2025-11-21
**Status**: Draft
**Author**: Claude Code (Refactor Planner Agent)

---

## Executive Summary

The rRNA-Phylo backend has grown organically to 3,500+ lines of working code with comprehensive DNA/RNA/Protein phylogenetic analysis capabilities. All 12 tests pass, and the functionality is production-ready. However, the codebase suffers from:

- Flat directory structure (26 files in c:\Users\User\Desktop\projects\rrna-phylo\backend\)
- Tests mixed with source code
- Code duplication (TreeNode, GTRModel, import inconsistencies)
- Partial migration artifacts (files exist in both old and new locations)
- No proper Python package structure

**Proposed Solution**: Systematic reorganization into a proper Python package (rrna_phylo) with clean module boundaries, eliminated duplication, and separated tests - all while preserving 100% functionality.

**Timeline**: 4-6 hours of careful, incremental work
**Risk Level**: Low-Medium (with proper testing at each step)
**Success Criteria**: All tests pass, cleaner imports, no duplication, proper package structure

---

## Table of Contents

1. [Current State Analysis](#current-state-analysis)
2. [Identified Issues](#identified-issues)
3. [Target Architecture](#target-architecture)
4. [Code Consolidation Opportunities](#code-consolidation-opportunities)
5. [Refactoring Phases](#refactoring-phases)
6. [Risk Assessment](#risk-assessment)
7. [Testing Strategy](#testing-strategy)
8. [Import Update Strategy](#import-update-strategy)
9. [Success Metrics](#success-metrics)

---

## Current State Analysis

### File Inventory

**Location**: c:\Users\User\Desktop\projects\rrna-phylo\backend\

#### Source Files (Root Level - 14 files)
```
aligner.py              (64 lines)   - MUSCLE wrapper
bionj.py                (282 lines)  - BioNJ tree builder
compare_trees.py        (205 lines)  - Tree comparison utility
distance.py             (193 lines)  - Distance calculation (JC)
fasta_parser.py         (212 lines)  - FASTA parser
ml_tree.py              (322 lines)  - ML Level 1 (GTR model)
ml_tree_level2.py       (443 lines)  - ML Level 2 (Felsenstein)
ml_tree_level3.py       (448 lines)  - ML Level 3 (GTR+Gamma)
phylo_builder.py        (339 lines)  - Unified builder
protein_distance.py     (291 lines)  - Protein distances
protein_ml.py           (429 lines)  - Protein ML
protein_models.py       (356 lines)  - Protein substitution models
sequence_type.py        (238 lines)  - Type detection
upgma.py                (188 lines)  - UPGMA builder
visualize_trees.py      (172 lines)  - Visualization
migrate_structure.py    (108 lines)  - Migration script
```

#### Test Files (Root Level - 12 files)
```
test_aligner.py         (132 lines)
test_bionj.py           (199 lines)
test_distance.py        (158 lines)
test_ml_level2.py       (212 lines)
test_ml_level3.py       (159 lines)
test_ml_tree.py         (138 lines)
test_parser.py          (165 lines)
test_phylo_builder.py   (296 lines)
test_protein_phylo.py   (257 lines)
test_sequence_type.py   (234 lines)
test_upgma.py           (163 lines)
```

#### Partial Migration (rrna_phylo/ directory exists)
```
rrna_phylo/
├── __init__.py                    (52 lines)  - Package root with exports
├── core/
│   ├── sequence_type.py           (238 lines) - DUPLICATE of root
│   └── tree.py                    (122 lines) - TreeNode class (extracted)
└── io/
    ├── __init__.py                (16 lines)
    ├── aligner.py                 (64 lines)  - DUPLICATE of root
    └── fasta_parser.py            (212 lines) - DUPLICATE of root
```

**Additional empty directories**: tests/, rrna_phylo/distance/, rrna_phylo/methods/, rrna_phylo/models/, rrna_phylo/utils/

### Code Architecture Patterns

#### Pattern 1: TreeNode Duplication
- **upgma.py** (lines 12-53): Full TreeNode class definition
- **rrna_phylo/core/tree.py** (lines 8-123): Enhanced TreeNode class with better docs
- **Problem**: Two implementations, old code imports from upgma.py
- **Impact**: Changes to TreeNode need to be made twice

#### Pattern 2: Import Inconsistency
```python
# Old style (most files)
from upgma import TreeNode
from fasta_parser import Sequence
from sequence_type import SequenceType

# New style (rrna_phylo/__init__.py)
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.sequence_type import SequenceType
```

#### Pattern 3: ML Evolution (Three Levels)
The ML implementation evolved through three iterations:
- **ml_tree.py** (Level 1): Basic GTR model, placeholder likelihood
- **ml_tree_level2.py** (Level 2): Full Felsenstein pruning, branch optimization
- **ml_tree_level3.py** (Level 3): GTR+Gamma, site pattern compression

**Issue**: Levels 1 & 2 are superseded but still exist. Level 3 imports from both Level 2 and Level 1.

#### Pattern 4: Sequence Type Detection
Three separate paths for DNA/RNA/Protein:
- DNA/RNA: Use GTR model (4x4 matrix), share code path
- Protein: Use empirical models (20x20 matrix), separate code path
- This is **correct architecture** (different biological models)

#### Pattern 5: Builder Pattern
phylo_builder.py provides unified interface:
```python
builder = PhylogeneticTreeBuilder()
builder.detect_and_validate(sequences)  # Auto-detect type
upgma = builder.build_upgma_tree(sequences)
bionj = builder.build_bionj_tree(sequences)
ml = builder.build_ml_tree(sequences)
```

### Dependency Analysis

#### Core Dependencies (No dependencies)
1. fasta_parser.py → (none)
2. tree.py (rrna_phylo/core/) → (none)

#### Second-Level Dependencies
3. sequence_type.py → fasta_parser
4. aligner.py → fasta_parser
5. distance.py → fasta_parser, sequence_type
6. protein_distance.py → fasta_parser

#### Tree Building Dependencies
7. upgma.py → (defines TreeNode internally)
8. bionj.py → distance, fasta_parser, upgma.TreeNode
9. ml_tree.py → fasta_parser, upgma.TreeNode, bionj, distance
10. ml_tree_level2.py → fasta_parser, upgma.TreeNode, ml_tree.GTRModel
11. ml_tree_level3.py → fasta_parser, upgma.TreeNode, ml_tree.GTRModel, ml_tree_level2.LikelihoodCalculator

#### Protein Dependencies
12. protein_models.py → fasta_parser
13. protein_ml.py → fasta_parser, upgma.TreeNode, protein_models, ml_tree_level3.GammaRates

#### High-Level Dependencies
14. phylo_builder.py → ALL of the above

### Test Coverage

All 12 test files use root-level imports:
```python
from fasta_parser import Sequence
from upgma import build_upgma_tree
from ml_tree_level3 import build_ml_tree_level3
```

Tests are comprehensive but not organized:
- No conftest.py for shared fixtures
- No test data directory
- Tests mixed with source code

---

## Identified Issues

### Critical Issues (Block Clean Architecture)

#### C1. TreeNode Duplication
- **Location**: upgma.py (lines 12-53) vs rrna_phylo/core/tree.py (lines 8-123)
- **Impact**: All tree-building code imports from upgma.py
- **Risk**: Changes to TreeNode class require updating two places
- **Severity**: High

#### C2. Partial Migration State
- **Location**: Files exist in both root and rrna_phylo/
- **Files**: aligner.py, fasta_parser.py, sequence_type.py
- **Impact**: Confusion about which version is canonical
- **Risk**: Edits to wrong file, diverging implementations
- **Severity**: High

#### C3. Tests Mixed with Source
- **Location**: All 12 test_*.py files in backend/
- **Impact**: Cannot distinguish source from tests, harder to package
- **Risk**: Tests get imported as code, confusing for users
- **Severity**: Medium

### Major Issues (Reduce Maintainability)

#### M1. ML Level Progression Confusion
- **Files**: ml_tree.py (Level 1), ml_tree_level2.py (Level 2), ml_tree_level3.py (Level 3)
- **Issue**: Level 1 & 2 are superseded but still imported by tests
- **Impact**: Dead code that must be maintained
- **Recommendation**: Archive Level 1 & 2, refactor Level 3 to be standalone
- **Severity**: Medium

#### M2. GTRModel Duplication
- **Location**: ml_tree.py defines GTRModel
- **Usage**: ml_tree_level2.py and ml_tree_level3.py import it
- **Issue**: GTRModel should be in models/ module, not in ml_tree.py
- **Severity**: Medium

#### M3. Flat Import Namespace
- **Current**: All modules at root level, no organization
- **Impact**: Hard to understand what's core vs utility vs method
- **Example**: `from upgma import build_upgma_tree` - unclear this is a tree method
- **Severity**: Medium

#### M4. No Package Metadata
- **Missing**: setup.py or pyproject.toml
- **Impact**: Cannot `pip install -e .`, cannot import from outside
- **Risk**: Cannot use as a library in other projects
- **Severity**: Low (future-proofing)

### Minor Issues (Quality of Life)

#### m1. Inconsistent File Naming
- Some files: `ml_tree_level3.py` (snake_case with version)
- Others: `phylo_builder.py` (snake_case, no version)
- Recommendation: Use consistent naming, drop version numbers

#### m2. Missing __init__.py in Tests
- tests/ directory exists but is empty
- Should have conftest.py for fixtures

#### m3. No Documentation Directory Structure
- README.md and ARCHITECTURE.md exist (good!)
- But no docs/ directory for user guides, API docs

---

## Target Architecture

### Proposed Directory Structure

```
backend/
├── rrna_phylo/                    # Main package
│   ├── __init__.py                # Public API exports
│   │
│   ├── core/                      # Core data structures
│   │   ├── __init__.py
│   │   ├── tree.py                # TreeNode (canonical)
│   │   └── sequence_type.py      # Type detection
│   │
│   ├── io/                        # Input/Output
│   │   ├── __init__.py
│   │   ├── fasta_parser.py        # FASTA parsing
│   │   └── aligner.py             # MUSCLE wrapper
│   │
│   ├── distance/                  # Distance calculations
│   │   ├── __init__.py
│   │   ├── nucleotide.py          # DNA/RNA distances (JC, K2P)
│   │   └── protein.py             # Protein distances (Poisson)
│   │
│   ├── models/                    # Substitution models
│   │   ├── __init__.py
│   │   ├── gtr.py                 # GTR model (DNA/RNA)
│   │   ├── gamma.py               # Gamma rate heterogeneity
│   │   └── protein_models.py      # WAG, LG, JTT
│   │
│   ├── methods/                   # Tree building methods
│   │   ├── __init__.py
│   │   ├── upgma.py               # UPGMA
│   │   ├── bionj.py               # BioNJ
│   │   ├── ml_nucleotide.py       # ML for DNA/RNA (GTR+Gamma)
│   │   └── ml_protein.py          # ML for Protein
│   │
│   ├── utils/                     # Utilities
│   │   ├── __init__.py
│   │   ├── visualization.py       # Tree visualization
│   │   └── comparison.py          # Tree comparison
│   │
│   └── builder.py                 # High-level unified interface
│
├── tests/                         # All tests (separate)
│   ├── __init__.py
│   ├── conftest.py                # Pytest fixtures
│   │
│   ├── core/
│   │   ├── test_tree.py
│   │   └── test_sequence_type.py
│   │
│   ├── io/
│   │   ├── test_fasta_parser.py
│   │   └── test_aligner.py
│   │
│   ├── distance/
│   │   ├── test_nucleotide.py
│   │   └── test_protein.py
│   │
│   ├── methods/
│   │   ├── test_upgma.py
│   │   ├── test_bionj.py
│   │   ├── test_ml_nucleotide.py
│   │   └── test_ml_protein.py
│   │
│   └── integration/
│       ├── test_builder.py         # End-to-end tests
│       └── test_dna_rna_protein.py # Full pipeline tests
│
├── examples/                      # Usage examples
│   ├── basic_usage.py
│   ├── 16s_rrna_example.py
│   └── protein_example.py
│
├── docs/                          # Documentation
│   ├── architecture.md            # Move ARCHITECTURE.md here
│   ├── api_reference.md
│   └── user_guide.md
│
├── pyproject.toml                 # Package metadata (modern)
├── README.md                      # Project overview
└── CHANGELOG.md                   # Version history
```

### Module Responsibilities

#### Core (rrna_phylo/core/)
- **tree.py**: TreeNode class, Newick format, tree utilities
- **sequence_type.py**: SequenceType enum, SequenceTypeDetector, model selection

#### IO (rrna_phylo/io/)
- **fasta_parser.py**: Sequence dataclass, FastaParser, validation
- **aligner.py**: MuscleAligner wrapper

#### Distance (rrna_phylo/distance/)
- **nucleotide.py**: DNA/RNA distance calculations (Jukes-Cantor, Kimura 2-parameter)
- **protein.py**: Protein distance calculations (Poisson, Kimura)

#### Models (rrna_phylo/models/)
- **gtr.py**: GTRModel class (rate matrix, probability matrices)
- **gamma.py**: GammaRates class (rate heterogeneity)
- **protein_models.py**: ProteinModel, WAG, LG, JTT matrices

#### Methods (rrna_phylo/methods/)
- **upgma.py**: UPGMABuilder, build_upgma_tree()
- **bionj.py**: BioNJBuilder, build_bionj_tree()
- **ml_nucleotide.py**: Maximum Likelihood for DNA/RNA (GTR+Gamma)
- **ml_protein.py**: Maximum Likelihood for Protein (WAG/LG/JTT)

#### Utils (rrna_phylo/utils/)
- **visualization.py**: Tree visualization, ASCII printing
- **comparison.py**: Tree comparison, consensus methods

#### Builder (rrna_phylo/builder.py)
- **PhylogeneticTreeBuilder**: Unified interface, auto-detection
- **build_trees()**: Convenience function for all three methods

### Public API (rrna_phylo/__init__.py)

Users should be able to:
```python
# High-level interface
from rrna_phylo import build_trees, Sequence, TreeNode, SequenceType

sequences = [
    Sequence("seq1", "Species 1", "ATGCAT"),
    Sequence("seq2", "Species 2", "ATGCAC"),
]

results = build_trees(sequences)
upgma_tree = results["upgma"]
bionj_tree = results["bionj"]
ml_tree, log_likelihood = results["ml"]
```

Advanced users can import specific modules:
```python
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.models.gtr import GTRModel
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
```

---

## Code Consolidation Opportunities

### Consolidation 1: TreeNode Class

**Current State**:
- upgma.py (lines 12-53): Original TreeNode implementation
- rrna_phylo/core/tree.py (lines 8-123): Enhanced TreeNode with more methods

**Proposed**:
- Use rrna_phylo/core/tree.py as canonical source
- Remove TreeNode from upgma.py
- Update all imports: `from rrna_phylo.core.tree import TreeNode`

**Files to Update** (11 files):
1. upgma.py - remove TreeNode class, import from core.tree
2. bionj.py - change `from upgma import TreeNode` to `from rrna_phylo.core.tree import TreeNode`
3. ml_tree.py - change import
4. ml_tree_level2.py - change import
5. ml_tree_level3.py - change import
6. protein_ml.py - change import
7. phylo_builder.py - already imports correctly (update if needed)
8. All 12 test files - update imports

**Benefits**:
- Single source of truth for TreeNode
- Enhanced TreeNode available everywhere (get_leaves(), count_leaves())
- Clearer architecture

**Risk**: Low (TreeNode API is stable)

### Consolidation 2: GTRModel Extraction

**Current State**:
- ml_tree.py (lines 28-189): GTRModel class definition
- ml_tree_level2.py: imports GTRModel from ml_tree
- ml_tree_level3.py: imports GTRModel from ml_tree

**Proposed**:
- Create rrna_phylo/models/gtr.py with GTRModel
- Remove GTRModel from ml_tree.py (keep only placeholder ML code or archive)
- Update imports in ml_tree_level2.py and ml_tree_level3.py

**Benefits**:
- Models separated from methods (clean architecture)
- GTRModel can be used independently
- Easier to test models separately

**Risk**: Low (GTRModel is well-defined)

### Consolidation 3: ML Level Progression

**Current State**:
- ml_tree.py: Level 1 (basic GTR, placeholder likelihood)
- ml_tree_level2.py: Level 2 (full Felsenstein, branch optimization)
- ml_tree_level3.py: Level 3 (GTR+Gamma, pattern compression)

**Proposed**: Two Options

**Option A: Archive Levels 1 & 2**
- Keep ml_tree_level3.py as the canonical ML implementation
- Move ml_tree.py → archive/ml_level1.py
- Move ml_tree_level2.py → archive/ml_level2.py
- Refactor ml_tree_level3.py to not depend on previous levels
- Rename ml_tree_level3.py → rrna_phylo/methods/ml_nucleotide.py

**Option B: Keep All Levels (Educational)**
- Keep all three levels (for learning/demonstration)
- Move to methods/_ml_level1.py, _ml_level2.py, ml_nucleotide.py
- Prefix with _ to indicate internal/educational
- Add README explaining progression

**Recommendation**: Option A (cleaner, production-ready)

**Benefits**:
- Single ML implementation (easier to maintain)
- No version numbers in filenames
- Clear what users should use

**Risk**: Medium (need to ensure Level 3 is fully standalone)

### Consolidation 4: Duplicate Files (Partial Migration)

**Current State**:
- aligner.py (root) vs rrna_phylo/io/aligner.py - IDENTICAL
- fasta_parser.py (root) vs rrna_phylo/io/fasta_parser.py - IDENTICAL
- sequence_type.py (root) vs rrna_phylo/core/sequence_type.py - ALMOST IDENTICAL (import paths differ)

**Proposed**:
1. Verify rrna_phylo versions are correct
2. Update imports in rrna_phylo versions (if needed)
3. Delete root versions
4. Update all other files to import from rrna_phylo

**Benefits**:
- No duplication
- Clear package structure
- Single source of truth

**Risk**: Low (just need to update imports)

### Consolidation 5: GammaRates and SitePatternCompressor

**Current State**:
- ml_tree_level3.py contains:
  - GammaRates class (lines 25-100+)
  - SitePatternCompressor class
  - LikelihoodCalculator class
  - Tree search code

**Proposed**:
- Extract GammaRates → rrna_phylo/models/gamma.py
- Keep SitePatternCompressor and LikelihoodCalculator in ml_nucleotide.py (they're method-specific)

**Benefits**:
- Gamma model can be reused (protein ML also needs it)
- Cleaner module boundaries
- Easier testing

**Risk**: Low (GammaRates is independent)

---

## Refactoring Phases

### Phase 0: Preparation (30 minutes)

**Objective**: Set up infrastructure and safeguards

#### Step 0.1: Backup Current State
```bash
cd c:\Users\User\Desktop\projects\rrna-phylo
git status  # Ensure clean working directory
git branch refactor/package-structure
git checkout refactor/package-structure
```

#### Step 0.2: Run All Tests (Baseline)
```bash
cd backend
pytest test_*.py -v --tb=short > test_results_baseline.txt
```
- Document which tests pass
- Save baseline for comparison

#### Step 0.3: Create pyproject.toml
```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "rrna-phylo"
version = "0.1.0"
description = "rRNA prediction and phylogenetic tree builder"
readme = "README.md"
requires-python = ">=3.9"
dependencies = [
    "numpy>=1.20",
    "scipy>=1.7",
    "biopython>=1.79",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "pytest-cov>=3.0",
]

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
```

#### Step 0.4: Verify Package Structure
```bash
# Check what exists
ls -la rrna_phylo/
ls -la rrna_phylo/core/
ls -la rrna_phylo/io/
```

**Deliverables**:
- ✅ Git branch created
- ✅ Baseline test results saved
- ✅ pyproject.toml created
- ✅ Current structure documented

---

### Phase 1: Core Module Consolidation (1.5 hours)

**Objective**: Establish canonical core modules and fix TreeNode duplication

#### Step 1.1: Consolidate TreeNode (20 minutes)

**Actions**:
1. Verify rrna_phylo/core/tree.py is complete
2. Remove TreeNode class from upgma.py (lines 12-53)
3. Add import at top of upgma.py: `from rrna_phylo.core.tree import TreeNode`
4. Update upgma.py docstring to note TreeNode is imported

**Test**:
```bash
pytest test_upgma.py -v
```

**Files Modified**:
- upgma.py (remove class, add import)

#### Step 1.2: Update TreeNode Imports in Methods (30 minutes)

**Files to Update**:
1. bionj.py: `from upgma import TreeNode` → `from rrna_phylo.core.tree import TreeNode`
2. ml_tree.py: same
3. ml_tree_level2.py: same
4. ml_tree_level3.py: same
5. protein_ml.py: same

**Test After Each File**:
```bash
pytest test_bionj.py -v
pytest test_ml_tree.py -v
pytest test_ml_level2.py -v
pytest test_ml_level3.py -v
pytest test_protein_phylo.py -v
```

**Files Modified**: 5 files

#### Step 1.3: Consolidate Core Modules (30 minutes)

**Actions**:
1. Compare sequence_type.py (root) vs rrna_phylo/core/sequence_type.py
2. Ensure rrna_phylo version has correct imports
3. Delete root sequence_type.py
4. Update files that import sequence_type

**Files that import sequence_type**:
- phylo_builder.py: `from sequence_type import ...` → `from rrna_phylo.core.sequence_type import ...`
- distance.py: same
- (check with grep)

**Test**:
```bash
pytest test_sequence_type.py -v
pytest test_phylo_builder.py -v
pytest test_distance.py -v
```

**Files Modified**:
- Delete: sequence_type.py (root)
- Update: phylo_builder.py, distance.py, others (check with grep)

#### Step 1.4: Update __init__.py Files (10 minutes)

**rrna_phylo/__init__.py**:
```python
"""
rRNA-Phylo: General-Purpose Phylogenetic Analysis
"""

__version__ = "0.1.0"

# Core exports
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
from rrna_phylo.io.fasta_parser import Sequence, FastaParser
from rrna_phylo.io.aligner import MuscleAligner

__all__ = [
    "TreeNode",
    "Sequence",
    "FastaParser",
    "MuscleAligner",
    "SequenceType",
    "SequenceTypeDetector",
    "__version__",
]
```

**Test**:
```bash
python -c "from rrna_phylo import TreeNode, Sequence, SequenceType; print('OK')"
```

**Phase 1 Checkpoint**:
```bash
pytest test_*.py -v --tb=short
# All tests should still pass
```

---

### Phase 2: IO Module Consolidation (45 minutes)

**Objective**: Remove duplicate IO files, establish rrna_phylo/io as canonical

#### Step 2.1: Consolidate fasta_parser.py (20 minutes)

**Actions**:
1. Verify rrna_phylo/io/fasta_parser.py is complete
2. Compare with root fasta_parser.py (should be identical)
3. Delete root fasta_parser.py
4. Update all imports

**Files that import fasta_parser** (check with grep):
```bash
cd backend
grep -l "from fasta_parser import" *.py
```

Update each:
- `from fasta_parser import ...` → `from rrna_phylo.io.fasta_parser import ...`

**Test After Each Update**:
```bash
pytest test_parser.py -v
pytest test_aligner.py -v
# ... test files that use Sequence
```

**Files Modified**:
- Delete: fasta_parser.py (root)
- Update: ~10-15 files (all that import Sequence)

#### Step 2.2: Consolidate aligner.py (15 minutes)

**Actions**:
1. Verify rrna_phylo/io/aligner.py is complete
2. Delete root aligner.py
3. Update imports in test files

**Test**:
```bash
pytest test_aligner.py -v
```

**Files Modified**:
- Delete: aligner.py (root)
- Update: test_aligner.py, any other files

#### Step 2.3: Update rrna_phylo/io/__init__.py (10 minutes)

```python
"""
Input/Output module for sequence handling.
"""

from rrna_phylo.io.fasta_parser import Sequence, FastaParser
from rrna_phylo.io.aligner import MuscleAligner

__all__ = [
    "Sequence",
    "FastaParser",
    "MuscleAligner",
]
```

**Phase 2 Checkpoint**:
```bash
pytest test_*.py -v --tb=short
# All tests should still pass
```

---

### Phase 3: Distance Module Organization (30 minutes)

**Objective**: Move distance calculations to rrna_phylo/distance/

#### Step 3.1: Move distance.py to rrna_phylo/distance/nucleotide.py

**Actions**:
```bash
cd backend
cp distance.py rrna_phylo/distance/nucleotide.py
```

**Update imports in rrna_phylo/distance/nucleotide.py**:
```python
# Old
from fasta_parser import Sequence
from sequence_type import SequenceType

# New
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.sequence_type import SequenceType
```

**Update files that import distance**:
```bash
grep -l "from distance import" *.py
```

Change to:
```python
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
```

**Test**:
```bash
pytest test_distance.py -v
pytest test_upgma.py -v
pytest test_bionj.py -v
```

**Delete root distance.py** after tests pass.

#### Step 3.2: Move protein_distance.py to rrna_phylo/distance/protein.py

**Similar process**:
1. Copy to rrna_phylo/distance/protein.py
2. Update imports inside file
3. Update files that import protein_distance
4. Test
5. Delete root file

**Test**:
```bash
pytest test_protein_phylo.py -v
```

#### Step 3.3: Create rrna_phylo/distance/__init__.py

```python
"""
Distance calculation methods for phylogenetics.
"""

from rrna_phylo.distance.nucleotide import calculate_distance_matrix
from rrna_phylo.distance.protein import calculate_protein_distance_matrix

__all__ = [
    "calculate_distance_matrix",
    "calculate_protein_distance_matrix",
]
```

**Phase 3 Checkpoint**:
```bash
pytest test_distance.py test_protein_phylo.py -v
```

---

### Phase 4: Models Module Organization (1 hour)

**Objective**: Extract substitution models to rrna_phylo/models/

#### Step 4.1: Extract GTRModel to rrna_phylo/models/gtr.py (25 minutes)

**Actions**:
1. Copy GTRModel class from ml_tree.py (lines 28-189) to new file rrna_phylo/models/gtr.py
2. Update imports in gtr.py:
   ```python
   from rrna_phylo.io.fasta_parser import Sequence
   ```
3. Update ml_tree.py:
   ```python
   from rrna_phylo.models.gtr import GTRModel
   ```
4. Update ml_tree_level2.py: same
5. Update ml_tree_level3.py: same

**Test**:
```bash
pytest test_ml_tree.py test_ml_level2.py test_ml_level3.py -v
```

**Files Modified**:
- Create: rrna_phylo/models/gtr.py
- Update: ml_tree.py (remove GTRModel class, add import)
- Update: ml_tree_level2.py, ml_tree_level3.py (update import)

#### Step 4.2: Extract GammaRates to rrna_phylo/models/gamma.py (20 minutes)

**Actions**:
1. Extract GammaRates class from ml_tree_level3.py to rrna_phylo/models/gamma.py
2. Update ml_tree_level3.py to import: `from rrna_phylo.models.gamma import GammaRates`
3. Update protein_ml.py: same

**Test**:
```bash
pytest test_ml_level3.py test_protein_phylo.py -v
```

#### Step 4.3: Move protein_models.py to rrna_phylo/models/protein_models.py (15 minutes)

**Actions**:
1. Copy protein_models.py to rrna_phylo/models/protein_models.py
2. Update imports inside file
3. Update protein_ml.py to import from new location
4. Test
5. Delete root file

**Test**:
```bash
pytest test_protein_phylo.py -v
```

#### Step 4.4: Create rrna_phylo/models/__init__.py

```python
"""
Substitution models for phylogenetic inference.
"""

from rrna_phylo.models.gtr import GTRModel
from rrna_phylo.models.gamma import GammaRates
from rrna_phylo.models.protein_models import ProteinModel, WAGModel, LGModel

__all__ = [
    "GTRModel",
    "GammaRates",
    "ProteinModel",
    "WAGModel",
    "LGModel",
]
```

**Phase 4 Checkpoint**:
```bash
pytest test_ml_tree.py test_ml_level2.py test_ml_level3.py test_protein_phylo.py -v
```

---

### Phase 5: Methods Module Organization (1 hour)

**Objective**: Move tree building methods to rrna_phylo/methods/

#### Step 5.1: Move upgma.py to rrna_phylo/methods/upgma.py (15 minutes)

**Actions**:
1. Copy upgma.py to rrna_phylo/methods/upgma.py
2. Update imports (TreeNode already imported correctly)
3. Update files that import upgma:
   ```python
   from rrna_phylo.methods.upgma import build_upgma_tree, UPGMABuilder
   ```
4. Test
5. Delete root upgma.py

**Test**:
```bash
pytest test_upgma.py -v
```

#### Step 5.2: Move bionj.py to rrna_phylo/methods/bionj.py (15 minutes)

**Similar process**:
1. Copy to rrna_phylo/methods/bionj.py
2. Update imports inside:
   ```python
   from rrna_phylo.core.tree import TreeNode
   from rrna_phylo.distance.nucleotide import calculate_distance_matrix
   from rrna_phylo.io.fasta_parser import Sequence
   ```
3. Update files that import bionj
4. Test
5. Delete root file

**Test**:
```bash
pytest test_bionj.py -v
```

#### Step 5.3: Move ML implementations (30 minutes)

**Decision Point**: Archive or refactor?

**Option A - Archive Levels 1 & 2** (RECOMMENDED):
1. Create backend/archive/ directory
2. Move ml_tree.py → archive/ml_level1.py
3. Move ml_tree_level2.py → archive/ml_level2.py
4. Move ml_tree_level3.py → rrna_phylo/methods/ml_nucleotide.py
5. Refactor ml_nucleotide.py to be self-contained (no imports from archived files)
6. Update test files

**Option B - Keep All Levels**:
1. Move ml_tree.py → rrna_phylo/methods/_ml_level1.py
2. Move ml_tree_level2.py → rrna_phylo/methods/_ml_level2.py
3. Move ml_tree_level3.py → rrna_phylo/methods/ml_nucleotide.py
4. Update all imports

**Recommended**: Option A

**Actions for Option A**:
```bash
mkdir -p backend/archive
mv ml_tree.py archive/ml_level1.py
mv ml_tree_level2.py archive/ml_level2.py
cp ml_tree_level3.py rrna_phylo/methods/ml_nucleotide.py
```

**Update rrna_phylo/methods/ml_nucleotide.py**:
- Ensure all classes are self-contained
- Update imports to use rrna_phylo paths
- If it imports from ml_tree_level2, copy needed classes

**Test**:
```bash
pytest test_ml_level3.py -v
```

#### Step 5.4: Move protein_ml.py to rrna_phylo/methods/ml_protein.py

**Similar process**:
1. Copy to rrna_phylo/methods/ml_protein.py
2. Update imports
3. Test
4. Delete root file

**Test**:
```bash
pytest test_protein_phylo.py -v
```

#### Step 5.5: Create rrna_phylo/methods/__init__.py

```python
"""
Tree building methods for phylogenetic inference.
"""

from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.ml_nucleotide import build_ml_tree_level3 as build_ml_tree
from rrna_phylo.methods.ml_protein import build_protein_ml_tree

__all__ = [
    "build_upgma_tree",
    "build_bionj_tree",
    "build_ml_tree",
    "build_protein_ml_tree",
]
```

**Phase 5 Checkpoint**:
```bash
pytest test_upgma.py test_bionj.py test_ml_level3.py test_protein_phylo.py -v
```

---

### Phase 6: Utilities and Builder (30 minutes)

**Objective**: Complete package structure

#### Step 6.1: Move utilities to rrna_phylo/utils/

**Actions**:
1. Copy visualize_trees.py → rrna_phylo/utils/visualization.py
2. Copy compare_trees.py → rrna_phylo/utils/comparison.py
3. Update imports in both files
4. Delete root files

**No tests for these** (they're utility scripts, not tested)

#### Step 6.2: Move phylo_builder.py to rrna_phylo/builder.py

**Actions**:
1. Copy phylo_builder.py → rrna_phylo/builder.py
2. Update imports to use rrna_phylo paths:
   ```python
   from rrna_phylo.io.fasta_parser import Sequence
   from rrna_phylo.core.sequence_type import SequenceTypeDetector, SequenceType
   from rrna_phylo.methods.upgma import build_upgma_tree
   from rrna_phylo.methods.bionj import build_bionj_tree
   from rrna_phylo.methods.ml_nucleotide import build_ml_tree_level3
   from rrna_phylo.methods.ml_protein import build_protein_ml_tree
   from rrna_phylo.distance.nucleotide import calculate_distance_matrix
   from rrna_phylo.distance.protein import calculate_protein_distance_matrix
   ```
3. Test
4. Delete root phylo_builder.py

**Test**:
```bash
pytest test_phylo_builder.py -v
```

#### Step 6.3: Update main rrna_phylo/__init__.py

```python
"""
rRNA-Phylo: General-Purpose Phylogenetic Analysis
"""

__version__ = "0.1.0"
__author__ = "rRNA-Phylo Contributors"

# Core exports - main user-facing API
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
from rrna_phylo.builder import PhylogeneticTreeBuilder, build_trees
from rrna_phylo.io.fasta_parser import Sequence, FastaParser

__all__ = [
    # Main functions
    "build_trees",

    # Core classes
    "TreeNode",
    "Sequence",
    "FastaParser",
    "SequenceType",
    "SequenceTypeDetector",
    "PhylogeneticTreeBuilder",

    # Version
    "__version__",
]
```

**Test**:
```bash
python -c "from rrna_phylo import build_trees; print('OK')"
```

**Phase 6 Checkpoint**:
```bash
pytest test_phylo_builder.py -v
```

---

### Phase 7: Test Reorganization (1 hour)

**Objective**: Move tests to tests/ directory, organize by module

#### Step 7.1: Create tests structure

```bash
mkdir -p tests/core
mkdir -p tests/io
mkdir -p tests/distance
mkdir -p tests/methods
mkdir -p tests/integration
```

#### Step 7.2: Move and rename test files

**Core tests**:
```bash
cp test_sequence_type.py tests/core/test_sequence_type.py
```

**IO tests**:
```bash
cp test_parser.py tests/io/test_fasta_parser.py
cp test_aligner.py tests/io/test_aligner.py
```

**Distance tests**:
```bash
cp test_distance.py tests/distance/test_nucleotide.py
```

**Methods tests**:
```bash
cp test_upgma.py tests/methods/test_upgma.py
cp test_bionj.py tests/methods/test_bionj.py
cp test_ml_level3.py tests/methods/test_ml_nucleotide.py
cp test_protein_phylo.py tests/methods/test_protein.py
```

**Integration tests**:
```bash
cp test_phylo_builder.py tests/integration/test_builder.py
```

**Archive old ML tests**:
```bash
mv test_ml_tree.py archive/test_ml_level1.py
mv test_ml_level2.py archive/test_ml_level2.py
```

#### Step 7.3: Update imports in all test files

**Pattern**:
```python
# Old
from fasta_parser import Sequence
from upgma import build_upgma_tree

# New
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.methods.upgma import build_upgma_tree
```

**Apply to all test files in tests/**

#### Step 7.4: Create conftest.py

**tests/conftest.py**:
```python
"""
Shared pytest fixtures for rRNA-Phylo tests.
"""

import pytest
from rrna_phylo.io.fasta_parser import Sequence


@pytest.fixture
def simple_dna_sequences():
    """Simple DNA sequences for testing."""
    return [
        Sequence("seq1", "Species 1", "ATGCAT"),
        Sequence("seq2", "Species 2", "ATGCAC"),
        Sequence("seq3", "Species 3", "ATCCAA"),
    ]


@pytest.fixture
def simple_rna_sequences():
    """Simple RNA sequences for testing."""
    return [
        Sequence("rna1", "16S rRNA 1", "AUGCAU"),
        Sequence("rna2", "16S rRNA 2", "AUGCAC"),
        Sequence("rna3", "16S rRNA 3", "AUCCAA"),
    ]


@pytest.fixture
def simple_protein_sequences():
    """Simple protein sequences for testing."""
    return [
        Sequence("prot1", "Protein 1", "MKTAYIAK"),
        Sequence("prot2", "Protein 2", "MKTAYIAQ"),
        Sequence("prot3", "Protein 3", "MSTAYIAK"),
    ]
```

#### Step 7.5: Test from new location

```bash
cd backend
pytest tests/ -v --tb=short
```

**If all pass, delete old test files from root**:
```bash
rm test_*.py
```

**Phase 7 Checkpoint**:
```bash
pytest tests/ -v --tb=short
# All tests should pass from new location
```

---

### Phase 8: Final Cleanup (30 minutes)

**Objective**: Remove old files, update documentation, verify everything works

#### Step 8.1: Remove old root files

**Verify all these are moved/deleted**:
```bash
# Should be moved to archive/
ls -la archive/

# Should be moved to rrna_phylo/
ls -la rrna_phylo/core/
ls -la rrna_phylo/io/
ls -la rrna_phylo/distance/
ls -la rrna_phylo/models/
ls -la rrna_phylo/methods/
ls -la rrna_phylo/utils/

# Should be moved to tests/
ls -la tests/

# Root should only have:
# - rrna_phylo/
# - tests/
# - examples/
# - docs/
# - archive/
# - README.md
# - pyproject.toml
# - ARCHITECTURE.md (move to docs/)
# - migrate_structure.py (delete or move to archive/)
```

**Delete**:
```bash
rm migrate_structure.py  # No longer needed
```

**Move**:
```bash
mv ARCHITECTURE.md docs/architecture.md
mv README.md docs/README.md  # Or keep in root
```

#### Step 8.2: Create examples directory

**examples/basic_usage.py**:
```python
"""
Basic usage example for rRNA-Phylo.
"""

from rrna_phylo import build_trees, Sequence

# Example sequences (DNA)
sequences = [
    Sequence("ecoli", "E. coli", "ATGCATGCATGC"),
    Sequence("salm", "Salmonella", "ATGCATGCATCC"),
    Sequence("bacil", "B. subtilis", "ATCCATGCATGC"),
]

# Build trees with all three methods
results = build_trees(sequences)

# Access results
upgma_tree = results["upgma"]
bionj_tree = results["bionj"]
ml_tree, log_likelihood = results["ml"]

# Print Newick format
print("UPGMA:", upgma_tree.to_newick())
print("BioNJ:", bionj_tree.to_newick())
print("ML:", ml_tree.to_newick())
print("ML Log-likelihood:", log_likelihood)
```

**examples/16s_rrna_example.py**:
```python
"""
16S rRNA phylogenetics example.
"""

from rrna_phylo import build_trees
from rrna_phylo.io import FastaParser

# Parse 16S rRNA sequences
parser = FastaParser()
sequences = parser.parse("16S_sequences.fasta")

# Build trees (auto-detects RNA)
results = build_trees(sequences)

# Compare methods
print("=" * 70)
print("16S rRNA PHYLOGENETIC ANALYSIS")
print("=" * 70)
print()
print("UPGMA:", results["upgma"].to_newick())
print("BioNJ:", results["bionj"].to_newick())
print("ML:", results["ml"][0].to_newick())
print(f"ML Log-likelihood: {results['ml'][1]:.2f}")
```

#### Step 8.3: Update README.md

**Add installation and usage sections**:
```markdown
## Installation

```bash
cd backend
pip install -e .
```

## Quick Start

```python
from rrna_phylo import build_trees, Sequence

sequences = [
    Sequence("seq1", "Species 1", "ATGCAT"),
    Sequence("seq2", "Species 2", "ATGCAC"),
    Sequence("seq3", "Species 3", "ATCCAA"),
]

results = build_trees(sequences)
print(results["upgma"].to_newick())
```

## Project Structure

```
backend/
├── rrna_phylo/           # Main package
│   ├── core/             # Core data structures
│   ├── io/               # Input/Output
│   ├── distance/         # Distance calculations
│   ├── models/           # Substitution models
│   ├── methods/          # Tree building methods
│   ├── utils/            # Utilities
│   └── builder.py        # High-level interface
├── tests/                # All tests
├── examples/             # Usage examples
└── docs/                 # Documentation
```
```

#### Step 8.4: Verify package works

**Test package import**:
```bash
cd backend
python -c "import rrna_phylo; print(rrna_phylo.__version__)"
python -c "from rrna_phylo import build_trees, Sequence, TreeNode; print('OK')"
```

**Run all tests**:
```bash
pytest tests/ -v --tb=short --cov=rrna_phylo
```

**Test examples**:
```bash
python examples/basic_usage.py
```

#### Step 8.5: Compare with baseline

```bash
# Compare test results
diff test_results_baseline.txt test_results_final.txt

# Should show same number of tests passing
# (maybe different file names, but all pass)
```

**Phase 8 Complete**: Refactoring done!

---

## Risk Assessment

### High-Risk Areas

#### R1. Import Chain Breakage
- **Risk**: Updating imports in wrong order causes cascading failures
- **Mitigation**:
  - Update bottom-up (low-level modules first)
  - Test after each file update
  - Use git to checkpoint after each phase
- **Recovery**: `git checkout HEAD -- <file>` to revert

#### R2. TreeNode Import Circular Dependencies
- **Risk**: Some files might have circular imports we didn't detect
- **Mitigation**:
  - TreeNode has no dependencies (pure data structure)
  - Import analysis shows no circularity
  - Test imports in isolation
- **Recovery**: Keep TreeNode in separate file with no imports

#### R3. Test File Import Paths
- **Risk**: Tests have many imports, easy to miss one
- **Mitigation**:
  - Update imports systematically
  - Test each test file individually
  - Use grep to find all imports before starting
- **Recovery**: Compare with baseline, fix missed imports

### Medium-Risk Areas

#### R4. ML Level Dependencies
- **Risk**: ml_tree_level3.py imports from level2, hard to untangle
- **Mitigation**:
  - Option A: Archive levels 1 & 2, make level 3 standalone
  - Option B: Keep all levels, update imports
  - Review dependencies before deciding
- **Recovery**: Keep all three levels if refactoring is too complex

#### R5. Protein Module Completeness
- **Risk**: Protein support is marked "in development", might be incomplete
- **Mitigation**:
  - Test protein_phylo tests pass
  - If any fail, keep protein modules separate or marked as experimental
- **Recovery**: Mark protein modules as rrna_phylo.experimental.protein

#### R6. External Tool Dependencies (MUSCLE)
- **Risk**: aligner.py calls external MUSCLE, might have path issues
- **Mitigation**:
  - Test aligner tests thoroughly
  - Ensure MUSCLE path is configurable
- **Recovery**: Keep path detection logic robust

### Low-Risk Areas

#### R7. File Moves
- **Risk**: Files get lost or corrupted during moves
- **Mitigation**: Copy first, test, then delete (never move directly)
- **Recovery**: Files in git history

#### R8. Newick Format Changes
- **Risk**: Tree serialization might differ
- **Mitigation**: TreeNode.to_newick() is well-tested, unchanged
- **Recovery**: N/A (very unlikely)

### Risk Mitigation Checklist

Before starting:
- ✅ Git branch created
- ✅ All tests pass (baseline)
- ✅ No uncommitted changes

During refactoring:
- ✅ Copy files, never move directly
- ✅ Test after each file update
- ✅ Commit after each phase
- ✅ Use descriptive commit messages

If something breaks:
1. Identify which phase caused the issue
2. Revert to last working commit
3. Re-attempt that phase more carefully
4. Document what went wrong

---

## Testing Strategy

### Test Levels

#### Unit Tests (rrna_phylo package)
- Core: TreeNode, SequenceType
- IO: FastaParser, MuscleAligner
- Distance: nucleotide, protein calculations
- Models: GTRModel, GammaRates, ProteinModel
- Methods: UPGMA, BioNJ, ML builders

#### Integration Tests (cross-module)
- Builder: PhylogeneticTreeBuilder with all methods
- End-to-end: Parse → Align → Detect → Build → Compare

#### Regression Tests (functionality preserved)
- All existing test cases should pass
- Same tree topologies produced
- Same Newick strings generated
- Same log-likelihoods calculated

### Testing Workflow

**After Each File Modification**:
```bash
pytest tests/<specific_test>.py -v
```

**After Each Phase**:
```bash
pytest tests/ -v --tb=short
```

**Final Verification**:
```bash
# All tests
pytest tests/ -v --tb=short --cov=rrna_phylo --cov-report=html

# Specific test suites
pytest tests/core/ -v
pytest tests/io/ -v
pytest tests/distance/ -v
pytest tests/models/ -v
pytest tests/methods/ -v
pytest tests/integration/ -v

# Test package imports
python -c "from rrna_phylo import build_trees, Sequence, TreeNode"
python -c "from rrna_phylo.methods import build_upgma_tree, build_bionj_tree"

# Run examples
python examples/basic_usage.py
python examples/16s_rrna_example.py
```

### Test Coverage Goals

- Overall: >80% coverage
- Core: >95% coverage (critical code)
- Methods: >80% coverage (complex algorithms)
- Utils: >60% coverage (utility code)

### Continuous Testing

**During refactoring, run after every change**:
```bash
# Quick smoke test (30 seconds)
pytest tests/core/ tests/io/ -v

# Full test suite (2-3 minutes)
pytest tests/ -v
```

**Git hooks** (optional):
```bash
# .git/hooks/pre-commit
#!/bin/bash
cd backend
pytest tests/ --tb=short -q || exit 1
```

---

## Import Update Strategy

### Import Update Patterns

#### Pattern 1: Low-Level Modules (No Dependencies)
**Files**: fasta_parser.py, tree.py
**Action**: These have no imports to update, just verify they work

#### Pattern 2: Core Modules (Depend on IO)
**Files**: sequence_type.py, aligner.py
**Old Import**:
```python
from fasta_parser import Sequence
```
**New Import**:
```python
from rrna_phylo.io.fasta_parser import Sequence
```

#### Pattern 3: Distance Modules
**Files**: distance.py, protein_distance.py
**Old Imports**:
```python
from fasta_parser import Sequence
from sequence_type import SequenceType
```
**New Imports**:
```python
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.sequence_type import SequenceType
```

#### Pattern 4: Model Modules
**Files**: ml_tree.py (GTRModel), protein_models.py
**Old Imports**:
```python
from fasta_parser import Sequence
from scipy.linalg import expm
```
**New Imports**:
```python
from rrna_phylo.io.fasta_parser import Sequence
from scipy.linalg import expm
```

#### Pattern 5: Method Modules (Most Complex)
**Files**: upgma.py, bionj.py, ml_tree_level3.py
**Old Imports**:
```python
from fasta_parser import Sequence
from upgma import TreeNode
from distance import calculate_distance_matrix
from ml_tree import GTRModel
```
**New Imports**:
```python
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
from rrna_phylo.models.gtr import GTRModel
```

#### Pattern 6: Builder (Imports Everything)
**File**: phylo_builder.py
**Old Imports** (12+ imports):
```python
from fasta_parser import Sequence
from sequence_type import SequenceTypeDetector, SequenceType
from upgma import TreeNode, build_upgma_tree
from bionj import build_bionj_tree
from distance import calculate_distance_matrix
from protein_distance import calculate_protein_distance_matrix
from ml_tree_level3 import build_ml_tree_level3
from protein_ml import build_protein_ml_tree
```
**New Imports**:
```python
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.sequence_type import SequenceTypeDetector, SequenceType
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.ml_nucleotide import build_ml_tree_level3
from rrna_phylo.methods.ml_protein import build_protein_ml_tree
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
from rrna_phylo.distance.protein import calculate_protein_distance_matrix
```

#### Pattern 7: Test Files
**Old Imports**:
```python
from fasta_parser import Sequence, FastaParser
from upgma import build_upgma_tree
from distance import calculate_distance_matrix
```
**New Imports**:
```python
from rrna_phylo.io.fasta_parser import Sequence, FastaParser
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
```

### Import Update Order (Critical!)

**Follow this order to avoid breakage**:

1. **Core modules first** (no dependencies):
   - rrna_phylo/core/tree.py
   - rrna_phylo/io/fasta_parser.py

2. **Modules that depend on core**:
   - rrna_phylo/core/sequence_type.py (depends on fasta_parser)
   - rrna_phylo/io/aligner.py (depends on fasta_parser)

3. **Distance and models**:
   - rrna_phylo/distance/nucleotide.py
   - rrna_phylo/distance/protein.py
   - rrna_phylo/models/gtr.py
   - rrna_phylo/models/gamma.py
   - rrna_phylo/models/protein_models.py

4. **Methods** (depend on core, models, distance):
   - rrna_phylo/methods/upgma.py
   - rrna_phylo/methods/bionj.py
   - rrna_phylo/methods/ml_nucleotide.py
   - rrna_phylo/methods/ml_protein.py

5. **Builder** (depends on everything):
   - rrna_phylo/builder.py

6. **Tests** (depend on everything):
   - tests/core/
   - tests/io/
   - tests/distance/
   - tests/models/
   - tests/methods/
   - tests/integration/

### Automated Import Finding

**Before starting each phase**:
```bash
# Find all files that import a specific module
cd backend
grep -rn "from fasta_parser import" *.py
grep -rn "from upgma import" *.py
grep -rn "from distance import" *.py
grep -rn "from ml_tree import" *.py
grep -rn "from sequence_type import" *.py

# Save to file for reference
grep -rn "from [a-z_]* import" *.py > import_audit.txt
```

### Import Verification

**After updating imports in a file**:
```python
# Test the imports work
python -c "import rrna_phylo.core.tree"
python -c "from rrna_phylo.core.tree import TreeNode"
python -c "from rrna_phylo.methods.upgma import build_upgma_tree"

# Test the file itself
python rrna_phylo/methods/upgma.py  # Should run without errors
```

### Common Import Errors

#### Error 1: ModuleNotFoundError
```
ModuleNotFoundError: No module named 'fasta_parser'
```
**Solution**: Update import to `from rrna_phylo.io.fasta_parser import ...`

#### Error 2: Circular Import
```
ImportError: cannot import name 'TreeNode' from partially initialized module
```
**Solution**: Check for circular dependencies, TreeNode should have no imports

#### Error 3: Relative Import Error
```
ImportError: attempted relative import with no known parent package
```
**Solution**: Use absolute imports (from rrna_phylo...), not relative (from .core...)

---

## Success Metrics

### Quantitative Metrics

#### Code Organization
- ✅ 0 files in backend/ root (except package dirs and config)
- ✅ All source in rrna_phylo/ package
- ✅ All tests in tests/ directory
- ✅ Clear module boundaries (core, io, distance, models, methods, utils)

#### Code Quality
- ✅ 0 duplicate classes (TreeNode, GTRModel)
- ✅ 0 duplicate files (fasta_parser, aligner, sequence_type)
- ✅ <3 levels of import nesting
- ✅ 100% tests passing (same number as baseline)

#### Package Health
- ✅ Can install with `pip install -e .`
- ✅ Can import as `from rrna_phylo import ...`
- ✅ pyproject.toml with dependencies
- ✅ __init__.py in all packages with clear exports

### Qualitative Metrics

#### Developer Experience
- ✅ Clear where to add new distance methods (rrna_phylo/distance/)
- ✅ Clear where to add new tree methods (rrna_phylo/methods/)
- ✅ Clear where to add new models (rrna_phylo/models/)
- ✅ Tests organized by functionality

#### Code Maintainability
- ✅ Changes to TreeNode only need updating one file
- ✅ Changes to GTRModel only need updating one file
- ✅ Import paths reflect logical organization
- ✅ Related code grouped together

#### User Experience
- ✅ Simple high-level API: `from rrna_phylo import build_trees`
- ✅ Advanced API: `from rrna_phylo.methods.upgma import UPGMABuilder`
- ✅ Clear examples in examples/
- ✅ Updated documentation

### Acceptance Criteria

**Before declaring refactoring complete, verify**:

1. **All Tests Pass**
   ```bash
   pytest tests/ -v --tb=short
   # Should show 12 tests passing (same as baseline)
   ```

2. **Package Imports Work**
   ```python
   from rrna_phylo import build_trees, Sequence, TreeNode
   from rrna_phylo.methods import build_upgma_tree, build_bionj_tree
   from rrna_phylo.models import GTRModel, GammaRates
   ```

3. **Examples Run**
   ```bash
   python examples/basic_usage.py
   python examples/16s_rrna_example.py
   ```

4. **No Duplicates**
   ```bash
   # Should find only one TreeNode definition
   grep -rn "class TreeNode" rrna_phylo/
   # Output: rrna_phylo/core/tree.py:8:class TreeNode:

   # Should find only one GTRModel definition
   grep -rn "class GTRModel" rrna_phylo/
   # Output: rrna_phylo/models/gtr.py:X:class GTRModel:
   ```

5. **Clean Root Directory**
   ```bash
   ls backend/
   # Should only see:
   # - rrna_phylo/
   # - tests/
   # - examples/
   # - docs/
   # - archive/
   # - pyproject.toml
   # - README.md
   ```

6. **Documentation Updated**
   - ✅ README.md has installation instructions
   - ✅ README.md has updated structure diagram
   - ✅ ARCHITECTURE.md moved to docs/
   - ✅ Examples are up-to-date

### Performance Verification

**Ensure refactoring doesn't slow down code**:
```bash
# Run tests with timing
pytest tests/ -v --durations=10

# Compare with baseline
# Should be similar (±10%)
```

### Regression Verification

**Ensure same functionality**:
```bash
# Run specific validation
pytest tests/integration/test_builder.py -v

# Check tree outputs match
python -c "
from rrna_phylo import build_trees, Sequence
seqs = [
    Sequence('s1', 'Species 1', 'ATGCAT'),
    Sequence('s2', 'Species 2', 'ATGCAC'),
]
results = build_trees(seqs)
print(results['upgma'].to_newick())
"

# Compare with known-good output
```

---

## Rollback Plan

### If Phase Fails

**Immediate rollback**:
```bash
git status  # See what changed
git diff    # Review changes
git checkout HEAD -- <file>  # Revert specific file
git reset --hard  # Revert all changes (nuclear option)
```

**Return to last checkpoint**:
```bash
git log --oneline -10  # See recent commits
git checkout <commit_hash>  # Go back to working state
git checkout -b refactor/retry  # Start fresh
```

### If Tests Fail

**Diagnose**:
```bash
pytest tests/ -v --tb=long  # Full traceback
pytest tests/methods/test_upgma.py -v -s  # Specific test with output
```

**Common failures**:
1. **ImportError**: Update import paths
2. **AttributeError**: Check if class/function moved
3. **NameError**: Check if variable/class renamed

**Fix or rollback**:
- If quick fix: Update imports, retest
- If complex: Rollback phase, plan better

### If Everything Breaks

**Nuclear rollback**:
```bash
git checkout main  # Go back to main branch
git branch -D refactor/package-structure  # Delete failed attempt
git checkout -b refactor/package-structure-v2  # Start over
```

**Lessons learned**:
- Document what went wrong
- Update this plan with mitigations
- Consider smaller phases

---

## Post-Refactoring Tasks

### Documentation

1. **Update README.md**:
   - Installation instructions
   - Package structure diagram
   - Import examples

2. **Update ARCHITECTURE.md** (now docs/architecture.md):
   - New file organization
   - Module responsibilities
   - Import patterns

3. **Create API Documentation**:
   ```bash
   pip install sphinx
   sphinx-quickstart docs/
   sphinx-apidoc -o docs/api rrna_phylo/
   cd docs && make html
   ```

4. **Write Migration Guide** (for existing users):
   - Old imports → New imports mapping
   - Breaking changes
   - Deprecation warnings

### Code Quality

1. **Add Type Hints** (optional):
   ```python
   from typing import List, Tuple, Optional
   def build_tree(sequences: List[Sequence]) -> TreeNode:
       ...
   ```

2. **Add Docstrings** (ensure all public functions have them):
   ```python
   def build_trees(sequences: List[Sequence]) -> dict:
       """
       Build phylogenetic trees using all three methods.

       Args:
           sequences: List of aligned sequences

       Returns:
           Dictionary with keys "upgma", "bionj", "ml"
       """
   ```

3. **Linting**:
   ```bash
   pip install ruff
   ruff check rrna_phylo/
   ruff format rrna_phylo/
   ```

### Testing Enhancements

1. **Add Coverage Badges**:
   ```bash
   pytest tests/ --cov=rrna_phylo --cov-report=html
   # Upload to codecov or similar
   ```

2. **Add Integration Tests**:
   - Real 16S rRNA data
   - Real protein data
   - Large datasets (performance)

3. **Add Benchmark Tests**:
   ```python
   def test_ml_performance(benchmark):
       benchmark(build_ml_tree, sequences)
   ```

### Packaging

1. **Publish to TestPyPI**:
   ```bash
   pip install build twine
   python -m build
   twine upload --repository testpypi dist/*
   ```

2. **Publish to PyPI** (when ready):
   ```bash
   twine upload dist/*
   ```

3. **Create Release**:
   - Tag version: `git tag v0.1.0`
   - Create GitHub release
   - Add changelog

---

## Timeline and Effort Estimate

### Estimated Timeline

| Phase | Description | Time Estimate | Cumulative |
|-------|-------------|---------------|------------|
| Phase 0 | Preparation | 30 min | 0.5 hr |
| Phase 1 | Core Consolidation | 1.5 hr | 2.0 hr |
| Phase 2 | IO Consolidation | 45 min | 2.75 hr |
| Phase 3 | Distance Organization | 30 min | 3.25 hr |
| Phase 4 | Models Organization | 1 hr | 4.25 hr |
| Phase 5 | Methods Organization | 1 hr | 5.25 hr |
| Phase 6 | Utils and Builder | 30 min | 5.75 hr |
| Phase 7 | Test Reorganization | 1 hr | 6.75 hr |
| Phase 8 | Final Cleanup | 30 min | 7.25 hr |
| **Total** | **All Phases** | **~7-8 hours** | |

### Parallel Work Opportunities

Some phases can be done in parallel (if multiple developers):
- Phase 2 & 3 (IO and Distance) - independent
- Phase 4 & 5 (Models and Methods) - models first, then methods

### Suggested Schedule

**Day 1 (4 hours)**:
- Phase 0: Preparation (0.5 hr)
- Phase 1: Core Consolidation (1.5 hr)
- Phase 2: IO Consolidation (0.75 hr)
- Phase 3: Distance Organization (0.5 hr)
- Phase 4: Models Organization (start, 0.75 hr)

**Day 2 (3-4 hours)**:
- Phase 4: Models Organization (finish, 0.25 hr)
- Phase 5: Methods Organization (1 hr)
- Phase 6: Utils and Builder (0.5 hr)
- Phase 7: Test Reorganization (1 hr)
- Phase 8: Final Cleanup (0.5 hr)

### Contingency Time

Add 20-30% buffer for:
- Unexpected import issues
- Test failures requiring investigation
- Merge conflicts (if on team)
- Documentation updates

**Total with buffer**: 8-10 hours

---

## Conclusion

This refactoring plan transforms the rRNA-Phylo backend from a flat structure with 26 files to a well-organized Python package with clear module boundaries. The approach is:

- **Incremental**: Each phase is self-contained and testable
- **Safe**: Test after every change, git checkpoints after every phase
- **Reversible**: Clear rollback strategy if something breaks
- **Comprehensive**: Covers code, tests, documentation, and examples

**Key Benefits**:
1. **Maintainability**: Changes to TreeNode or GTRModel only need one file update
2. **Clarity**: Clear where to add new functionality
3. **Usability**: Clean import paths (`from rrna_phylo import build_trees`)
4. **Professionalism**: Proper package structure, can publish to PyPI
5. **Testability**: Tests separated, organized by module

**Success Criteria**:
- ✅ All 12 tests pass
- ✅ No code duplication
- ✅ Clean package structure
- ✅ Can import as library

**Estimated Effort**: 7-10 hours of focused work

This refactoring sets the foundation for future growth (Phase 2: ML integration, Phase 3: Multi-tree consensus, Phase 4: Web API) by establishing clear architectural boundaries and module responsibilities.

---

## Appendix A: File Mapping

### Complete File Moves

| Old Location | New Location | Status |
|--------------|--------------|--------|
| fasta_parser.py | rrna_phylo/io/fasta_parser.py | Exists (delete old) |
| aligner.py | rrna_phylo/io/aligner.py | Exists (delete old) |
| sequence_type.py | rrna_phylo/core/sequence_type.py | Exists (delete old) |
| upgma.py | rrna_phylo/methods/upgma.py | To move |
| bionj.py | rrna_phylo/methods/bionj.py | To move |
| distance.py | rrna_phylo/distance/nucleotide.py | To move |
| protein_distance.py | rrna_phylo/distance/protein.py | To move |
| ml_tree.py | archive/ml_level1.py | Archive |
| ml_tree_level2.py | archive/ml_level2.py | Archive |
| ml_tree_level3.py | rrna_phylo/methods/ml_nucleotide.py | To move |
| protein_models.py | rrna_phylo/models/protein_models.py | To move |
| protein_ml.py | rrna_phylo/methods/ml_protein.py | To move |
| phylo_builder.py | rrna_phylo/builder.py | To move |
| visualize_trees.py | rrna_phylo/utils/visualization.py | To move |
| compare_trees.py | rrna_phylo/utils/comparison.py | To move |

### Test File Moves

| Old Location | New Location |
|--------------|--------------|
| test_parser.py | tests/io/test_fasta_parser.py |
| test_aligner.py | tests/io/test_aligner.py |
| test_sequence_type.py | tests/core/test_sequence_type.py |
| test_distance.py | tests/distance/test_nucleotide.py |
| test_upgma.py | tests/methods/test_upgma.py |
| test_bionj.py | tests/methods/test_bionj.py |
| test_ml_level3.py | tests/methods/test_ml_nucleotide.py |
| test_protein_phylo.py | tests/methods/test_protein.py |
| test_phylo_builder.py | tests/integration/test_builder.py |
| test_ml_tree.py | archive/test_ml_level1.py |
| test_ml_level2.py | archive/test_ml_level2.py |

---

## Appendix B: Import Cheat Sheet

### Quick Reference

**Sequence and Parser**:
```python
from rrna_phylo.io.fasta_parser import Sequence, FastaParser
```

**Tree Node**:
```python
from rrna_phylo.core.tree import TreeNode
```

**Sequence Type**:
```python
from rrna_phylo.core.sequence_type import SequenceType, SequenceTypeDetector
```

**Distance Calculations**:
```python
from rrna_phylo.distance.nucleotide import calculate_distance_matrix
from rrna_phylo.distance.protein import calculate_protein_distance_matrix
```

**Models**:
```python
from rrna_phylo.models.gtr import GTRModel
from rrna_phylo.models.gamma import GammaRates
from rrna_phylo.models.protein_models import WAGModel, LGModel
```

**Tree Methods**:
```python
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.methods.ml_nucleotide import build_ml_tree_level3
from rrna_phylo.methods.ml_protein import build_protein_ml_tree
```

**High-Level API**:
```python
from rrna_phylo import build_trees, Sequence, TreeNode, SequenceType
```

---

## Appendix C: Git Workflow

### Recommended Commits

**Phase 0**:
```bash
git commit -m "refactor: Add pyproject.toml and test baseline"
```

**Phase 1**:
```bash
git commit -m "refactor: Consolidate TreeNode to rrna_phylo.core.tree"
git commit -m "refactor: Update TreeNode imports in all modules"
git commit -m "refactor: Remove duplicate sequence_type.py"
```

**Phase 2**:
```bash
git commit -m "refactor: Consolidate IO modules (fasta_parser, aligner)"
```

**Phase 3**:
```bash
git commit -m "refactor: Move distance calculations to rrna_phylo.distance"
```

**Phase 4**:
```bash
git commit -m "refactor: Extract GTRModel to rrna_phylo.models.gtr"
git commit -m "refactor: Extract GammaRates to rrna_phylo.models.gamma"
git commit -m "refactor: Move protein models to rrna_phylo.models"
```

**Phase 5**:
```bash
git commit -m "refactor: Move tree methods to rrna_phylo.methods"
git commit -m "refactor: Archive ML levels 1 and 2"
```

**Phase 6**:
```bash
git commit -m "refactor: Move utilities to rrna_phylo.utils"
git commit -m "refactor: Move phylo_builder to rrna_phylo.builder"
```

**Phase 7**:
```bash
git commit -m "refactor: Reorganize tests into tests/ directory"
git commit -m "refactor: Update test imports for new package structure"
```

**Phase 8**:
```bash
git commit -m "refactor: Clean up root directory and archive old files"
git commit -m "docs: Update README and ARCHITECTURE for new structure"
git commit -m "chore: Add examples and finalize package structure"
```

**Final**:
```bash
git tag -a v0.1.0 -m "Version 0.1.0 - Package restructuring complete"
```

---

**End of Refactoring Plan**

This plan is ready to execute. Follow phases in order, test thoroughly, and commit frequently. Good luck!

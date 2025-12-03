# rRNA-Phylo Architecture Refactoring Plan

**Date**: 2025-12-03
**Version**: 1.0
**Status**: Ready for Review

---

## Executive Summary

This refactoring plan addresses architectural organization, code cleanliness, and maintainability improvements for the rRNA-Phylo project after recent cleanup efforts. The project has a solid foundation with 44 Python files (~11,400 lines) organized into 8 packages, but several organizational issues need attention.

**Key Findings**:
- Well-structured core architecture with clear separation of concerns
- Several modules misplaced in `utils/` that belong elsewhere
- Opportunity to consolidate ML-related functionality (5,000+ lines across 13 files in `models/`)
- Some duplicate functionality between `core/builder.py` and `core/builder_smart.py`
- Generally clean code with good documentation, minimal comment cleanup needed

**Impact**: Medium effort, high value. Most changes are file moves and renames (low risk).

---

## Current State Analysis

### Package Structure

```
backend/rrna_phylo/
├── core/               # 4 files, ~1,000 lines  - Tree structures, builders, type detection
├── io/                 # 3 files, ~400 lines    - FASTA parsing, alignment
├── distance/           # 3 files, ~400 lines    - Distance calculations
├── methods/            # 4 files, ~900 lines    - Tree building algorithms
├── models/             # 13 files, ~5,000 lines - ML models, GPU, rate matrices
├── utils/              # 9 files, ~2,500 lines  - Mixed utilities (problematic)
├── consensus/          # 3 files, ~400 lines    - Tree comparison, bipartitions
├── visualization/      # 2 files, ~100 lines    - ETE3 visualization
└── cli.py + config.py  # 2 files, ~700 lines    - Command-line interface
```

### File Size Distribution

**Largest files** (potential refactoring targets):
1. `cli.py` - 686 lines (acceptable, handles many options)
2. `models/ml_tree_level3.py` - 593 lines (core ML implementation)
3. `models/gpu_likelihood_torch.py` - 538 lines (GPU acceleration)
4. `utils/strain_handler.py` - 512 lines (should be moved)
5. `models/model_selection.py` - 491 lines (could be split)

### Dependency Analysis

**Clean imports** - no circular dependencies detected:
- `core/` is foundational (tree structure, no dependencies)
- `io/` depends only on `core/`
- `distance/` depends on `io/` and `core/`
- `methods/` depends on `core/`, `io/`, `distance/`
- `models/` is mostly self-contained with some `methods/` usage
- `utils/` has scattered dependencies (indicates misplaced modules)

---

## Issues Identified

### Critical Issues (Must Fix)

#### 1. **Misplaced Modules in `utils/` Package**

**Problem**: The `utils/` package contains modules that aren't utilities:

| File | Lines | Should Be In | Reason |
|------|-------|--------------|---------|
| `strain_handler.py` | 512 | `preprocessing/` (new) | Core preprocessing logic, not a utility |
| `dataset_analyzer.py` | 294 | `analysis/` or `methods/` | Analyzes datasets for method selection |
| `bootstrap.py` | 398 | `analysis/` (new) | Statistical analysis, not a utility |
| `sampling_strategy.py` | 350 | `preprocessing/` (new) | Dataset sampling, part of preprocessing |
| `outgroup_handler.py` | ~150 | `preprocessing/` (new) | Preprocessing logic |

**Impact**: Confusing organization, hard to find functionality.

#### 2. **`models/` Package Too Large and Unfocused**

**Problem**: 13 files totaling 5,000+ lines with mixed concerns:

**Subgroups identified**:
- **Core ML** (keep in `models/`): `ml_tree.py`, `ml_tree_level2.py`, `ml_tree_level3.py`, `ml_tree_level4.py`
- **Model Selection** (could be separate): `model_selection.py`, `substitution_models.py`, `rate_matrices.py`
- **Optimization** (could be separate): `branch_length_optimizer.py`, `tree_search.py`
- **Acceleration** (could be separate): `numba_likelihood.py`, `gpu_likelihood_torch.py`
- **Protein** (already separate): `protein_models.py`

**Current organization works but could be clearer with subpackages**.

#### 3. **Duplicate Builder Classes**

**Problem**: Two builder classes with overlapping functionality:
- `core/builder.py` - `PhylogeneticTreeBuilder` (451 lines)
- `core/builder_smart.py` - `SmartPhylogeneticTreeBuilder` (258 lines)

The "smart" builder extends the regular builder but adds complexity. Consider:
- Merge into single class with optional "smart" mode
- OR clearly document when to use each

#### 4. **`visualize_trees.py` in Wrong Package**

**Problem**: `utils/visualize_trees.py` contains visualization code but we have a `visualization/` package.

Should be: `visualization/ascii_viz.py` or similar.

### Major Issues (Should Fix)

#### 5. **Unclear ML Progression: Level 1-4**

**Problem**: Four ML implementations with progressive features:
- `ml_tree.py` (Level 1) - Basic GTR
- `ml_tree_level2.py` (Level 2) - Felsenstein algorithm
- `ml_tree_level3.py` (Level 3) - GTR+Gamma, compression
- `ml_tree_level4.py` (Level 4) - Model selection, NNI

**Why all four?** Only Level 4 is used in production. Others are kept for:
- Educational purposes (showing progression)
- Benchmarking
- Fallback if Level 4 breaks

**Recommendation**: Move Levels 1-3 to `models/deprecated/` or `models/reference/` with clear README explaining their purpose.

#### 6. **GPU Module Has "torch" in Name**

**Problem**: `gpu_likelihood_torch.py` - filename is implementation-specific.

Better: `gpu_likelihood.py` (implementation detail is internal).

### Minor Issues (Nice to Have)

#### 7. **`consensus/` Package Underutilized**

Only 3 files, ~400 lines. Consider:
- Move `tree_distance.py` functionality into `bipartitions.py`
- OR expand consensus methods (strict, majority-rule, etc.)

#### 8. **Config Module is Minimal**

`config.py` is very small. Consider merging into `__init__.py` or expanding with:
- Default parameters for all methods
- GPU detection settings
- Path configurations

---

## Proposed Refactoring Plan

### Phase 1: Package Reorganization (Priority: HIGH)

**Goal**: Move misplaced modules to appropriate packages.

#### Step 1.1: Create `preprocessing/` Package

**New structure**:
```
rrna_phylo/preprocessing/
├── __init__.py
├── strain_handler.py      # from utils/
├── outgroup_handler.py    # from utils/
├── sampling_strategy.py   # from utils/
└── deduplication.py       # future: extract from strain_handler
```

**Rationale**:
- Preprocessing is a major feature (40-60% dataset reduction)
- Currently scattered across "utils"
- Should be first-class package

**Files to move**:
1. `utils/strain_handler.py` → `preprocessing/strain_handler.py`
2. `utils/outgroup_handler.py` → `preprocessing/outgroup_handler.py`
3. `utils/sampling_strategy.py` → `preprocessing/sampling_strategy.py`

**Update imports in**:
- `cli.py` (uses all three)
- Any test files

**Risk**: LOW - Just file moves, imports easy to update.

#### Step 1.2: Create `analysis/` Package

**New structure**:
```
rrna_phylo/analysis/
├── __init__.py
├── bootstrap.py           # from utils/
└── dataset_analyzer.py    # from utils/
```

**Rationale**:
- Bootstrap is statistical analysis, not a utility
- Dataset analyzer is method selection analysis
- Both are analytical tools, not utilities

**Files to move**:
1. `utils/bootstrap.py` → `analysis/bootstrap.py`
2. `utils/dataset_analyzer.py` → `analysis/dataset_analyzer.py`

**Update imports in**:
- `cli.py`
- `core/builder_smart.py` (uses dataset_analyzer)

**Risk**: LOW - Clean moves, limited import updates.

#### Step 1.3: Move `visualize_trees.py`

**Action**:
```
utils/visualize_trees.py → visualization/ascii_viz.py
```

**Rationale**: Consolidate visualization in one package.

**Update imports in**:
- `utils/__init__.py`
- `cli.py`

**Risk**: VERY LOW - Single file move.

#### Step 1.4: Update `utils/` Package

**After moves, `utils/` should contain**:
- `console.py` (actual utility - console output formatting)

**Consider removing** `utils/` package entirely if only one file remains:
- Move `console.py` → `io/console.py` (I/O related)
- OR keep `utils/` as minimal utility package

**Risk**: LOW - Final cleanup.

---

### Phase 2: `models/` Package Restructuring (Priority: MEDIUM)

**Goal**: Organize ML code into logical subpackages.

#### Step 2.1: Create `models/ml/` Subpackage

**New structure**:
```
rrna_phylo/models/
├── ml/
│   ├── __init__.py
│   ├── ml_tree.py         # Level 1 (basic GTR)
│   ├── ml_tree_level2.py  # Level 2 (Felsenstein)
│   ├── ml_tree_level3.py  # Level 3 (production)
│   └── ml_tree_level4.py  # Level 4 (full features)
├── selection/
│   ├── __init__.py
│   ├── model_selection.py
│   ├── substitution_models.py
│   └── rate_matrices.py
├── optimization/
│   ├── __init__.py
│   ├── branch_length_optimizer.py
│   └── tree_search.py
├── acceleration/
│   ├── __init__.py
│   ├── numba_likelihood.py
│   └── gpu_likelihood.py  # renamed from gpu_likelihood_torch.py
├── protein_models.py      # stays at top level
└── __init__.py
```

**Rationale**:
- Clear separation of concerns
- Easier to navigate large codebase
- Logical grouping of related functionality

**Alternative (Simpler)**: Keep flat structure but rename files:
- `ml_tree_level3.py` → `ml_tree.py` (production)
- Others → `ml_tree_reference_level{1,2}.py`
- `ml_tree_level4.py` → `ml_tree_full.py` or `ml_tree_advanced.py`

**Recommendation**: Start with **Alternative** (simpler), consider subpackages if `models/` grows further.

**Risk**: MEDIUM - Many import changes, but clear benefits.

#### Step 2.2: Deprecate or Archive Levels 1-3

**Option A: Archive** (Recommended)
```
rrna_phylo/models/
├── reference/
│   ├── README.md          # Explains progression
│   ├── ml_tree_level1.py  # Educational: Basic GTR
│   ├── ml_tree_level2.py  # Educational: Felsenstein
│   └── ml_tree_level3.py  # Reference: Production without model selection
└── ml_tree_level4.py      # MAIN: Use this one
```

**Option B: Keep Current** (If used for testing)
- Keep all four if actively used for benchmarking
- Add clear docstrings explaining purpose of each

**Recommendation**: If Levels 1-3 aren't imported by production code → Archive.

**Risk**: LOW if not used in production.

---

### Phase 3: Builder Class Consolidation (Priority: MEDIUM)

**Goal**: Simplify builder interface.

#### Step 3.1: Merge or Clarify Builders

**Current situation**:
- `PhylogeneticTreeBuilder` - Main builder (451 lines)
- `SmartPhylogeneticTreeBuilder` - Extends with dataset analysis (258 lines)

**Option A: Merge** (Recommended)
```python
class PhylogeneticTreeBuilder:
    def __init__(self, verbose=True, smart_mode=True):
        self.smart_mode = smart_mode
        # ...

    def build_trees(self, sequences, method="all", alpha=1.0):
        if self.smart_mode:
            # Analyze dataset, skip unsuitable methods
            analysis = self.analyze_dataset(sequences)
            # Use recommendations...
        else:
            # Build all requested methods
            # ...
```

**Option B: Keep Separate but Document**
- Add clear docstrings explaining when to use each
- Update README to guide users

**Recommendation**: **Option B** (less risky, maintains backward compatibility).

**Risk**: LOW (documentation change) or MEDIUM (code merge).

---

### Phase 4: Code Quality Improvements (Priority: LOW)

**Goal**: Clean up minor issues.

#### Step 4.1: Standardize File Naming

**Current inconsistencies**:
- `gpu_likelihood_torch.py` - Implementation detail in name
- `ml_tree_level{1-4}.py` - Version numbers in name

**Recommendations**:
- `gpu_likelihood_torch.py` → `gpu_likelihood.py`
- Consider renaming ML levels (see Phase 2)

#### Step 4.2: Consolidate Small Packages

**`consensus/` package** (3 files, 400 lines):
- Option 1: Merge `tree_distance.py` into `bipartitions.py`
- Option 2: Keep as is (reasonable size)

**Recommendation**: Keep as is unless expanding consensus methods.

#### Step 4.3: Comment Cleanup

**Good news**: Code is generally well-documented with minimal unnecessary comments.

**No "# ..." comment patterns found** (excellent!).

**Minor cleanup needed**:
- Some files have overly verbose docstrings explaining obvious things
- Example patterns that could be simplified:
  ```python
  # Good (keep):
  def calculate_distance(self, seq1, seq2):
      """Calculate evolutionary distance between two sequences."""

  # Too verbose (simplify):
  def calculate_distance(self, seq1, seq2):
      """
      Calculate the evolutionary distance between two sequences.

      This function computes the pairwise distance using the specified
      model and applies corrections for saturation effects.

      The distance represents the expected number of substitutions per site.
      """
  ```

**Files with verbose docstrings** (review for simplification):
- `models/ml_tree_level3.py` - Some class docstrings are very long
- `models/substitution_models.py` - Detailed explanations (good for educational purposes)
- `utils/bootstrap.py` - Very detailed explanations

**Recommendation**: Keep detailed docstrings in complex ML code (educational value). Simplify only obvious utilities.

**Risk**: VERY LOW - Cosmetic changes.

---

## Implementation Order

### Sprint 1: High-Impact, Low-Risk Moves (Week 1)

**Focus**: Package reorganization (Phase 1)

1. Create `preprocessing/` package and move 3 files
2. Create `analysis/` package and move 2 files
3. Move `visualize_trees.py` to `visualization/`
4. Update all imports in affected files
5. Run full test suite to verify

**Deliverable**: Cleaner package structure, more intuitive organization.

**Estimated effort**: 4-8 hours

**Risk**: LOW

---

### Sprint 2: Builder Improvements (Week 2)

**Focus**: Builder documentation and minor refactoring (Phase 3)

1. Add comprehensive docstrings to both builders
2. Update README with guidance on when to use each
3. Add examples to docstrings
4. Consider merging if time permits (optional)

**Deliverable**: Clear builder interface, better documentation.

**Estimated effort**: 4-6 hours

**Risk**: LOW (documentation) to MEDIUM (merge)

---

### Sprint 3: Models Package Cleanup (Week 3)

**Focus**: ML code organization (Phase 2)

1. Archive or clearly mark reference implementations (Levels 1-3)
2. Rename `gpu_likelihood_torch.py` → `gpu_likelihood.py`
3. Update imports
4. Add README explaining ML progression
5. Consider subpackages if complexity warrants

**Deliverable**: Clearer ML code organization.

**Estimated effort**: 6-10 hours

**Risk**: MEDIUM (import updates)

---

### Sprint 4: Polish and Documentation (Week 4)

**Focus**: Final cleanup (Phase 4)

1. Review and simplify overly verbose docstrings (optional)
2. Update architecture documentation
3. Create package-level READMEs for new packages
4. Update main README with new structure

**Deliverable**: Polished, well-documented codebase.

**Estimated effort**: 4-6 hours

**Risk**: VERY LOW

---

## Risk Assessment

### Low-Risk Changes
- File moves within packages (import updates only)
- Documentation improvements
- File renames
- Adding new packages

### Medium-Risk Changes
- Moving files between packages (broader import updates)
- Merging builder classes (backward compatibility concerns)
- ML level reorganization (educational materials affected)

### High-Risk Changes
- None identified (all changes are organizational/cosmetic)

### Mitigation Strategies
1. **Test-driven refactoring**: Run test suite after each change
2. **Incremental commits**: One logical change per commit
3. **Backward compatibility**: Add deprecation warnings before removing old imports
4. **Branch protection**: Do refactoring in feature branch, review before merge

---

## Testing Strategy

### Pre-Refactoring
1. Document current test coverage (`pytest --cov`)
2. Create baseline test results
3. Identify critical paths that must not break

### During Refactoring
1. Run tests after each file move
2. Update tests to use new import paths
3. Verify no functionality regression

### Post-Refactoring
1. Full test suite must pass
2. Compare coverage (should be same or better)
3. Manual testing of CLI with real data
4. Verify all examples in README still work

### Test Files to Update
After package moves, these test files need import updates:
- `tests/test_strain_handler.py` (if exists)
- `tests/test_bootstrap.py` (if exists)
- `tests/test_builder.py` (if exists)
- Integration tests using moved modules

---

## Success Metrics

### Organization Metrics
- ✅ No modules in `utils/` that aren't utilities
- ✅ Package structure matches functional boundaries
- ✅ Clear separation between preprocessing, analysis, and core functionality

### Code Quality Metrics
- ✅ Test coverage maintained or improved (current: check with pytest)
- ✅ No broken imports in any file
- ✅ All CLI examples in README work

### Documentation Metrics
- ✅ Each package has clear purpose documented
- ✅ Builder usage is documented with examples
- ✅ Architecture diagram updated to match new structure

### User Experience Metrics
- ✅ Easier to find functionality (new users)
- ✅ Clear naming conventions throughout
- ✅ Logical package organization

---

## Backward Compatibility

### Import Aliases (Recommended)

To maintain backward compatibility during transition:

```python
# rrna_phylo/utils/__init__.py
from rrna_phylo.preprocessing.strain_handler import *
from rrna_phylo.analysis.bootstrap import *
import warnings

warnings.warn(
    "Importing from rrna_phylo.utils is deprecated. "
    "Use rrna_phylo.preprocessing or rrna_phylo.analysis instead.",
    DeprecationWarning,
    stacklevel=2
)
```

**Timeline**:
- Version 1.x: Add deprecation warnings
- Version 2.0: Remove old import paths

---

## Future Considerations

### Beyond This Refactoring

**When the project grows further**, consider:

1. **Further ML Subpackaging**:
   ```
   models/
   ├── ml/          # DNA/RNA ML
   ├── protein/     # Protein ML
   ├── selection/   # Model selection
   └── gpu/         # GPU acceleration
   ```

2. **API Layer**:
   - If adding web API, create `api/` package
   - Keep CLI separate from API

3. **Plugin Architecture**:
   - Allow custom substitution models
   - Allow custom tree building methods
   - Registry pattern for extensibility

4. **Performance Monitoring**:
   - Add `profiling/` package
   - Built-in performance benchmarks

---

## Appendix A: File Move Checklist

Use this checklist for each file move:

- [ ] Create destination directory (if new)
- [ ] Move file to new location
- [ ] Update imports in moved file (if any)
- [ ] Find all files importing the moved file: `grep -r "from.*{old_path}" .`
- [ ] Update each importing file
- [ ] Update package `__init__.py` files
- [ ] Run tests
- [ ] Update documentation
- [ ] Commit with descriptive message

---

## Appendix B: Detailed File Inventory

### `core/` Package (Keep as is - well organized)
- `tree.py` (239 lines) - Core TreeNode class ✅
- `sequence_type.py` (239 lines) - Type detection ✅
- `builder.py` (451 lines) - Main builder ⚠️ (see Phase 3)
- `builder_smart.py` (258 lines) - Smart builder ⚠️ (see Phase 3)

### `io/` Package (Keep as is - well organized)
- `fasta_parser.py` (331 lines) - FASTA parsing ✅
- `aligner.py` - MUSCLE wrapper ✅

### `distance/` Package (Keep as is - well organized)
- `distance.py` - Jukes-Cantor distances ✅
- `protein_distance.py` - Protein distances ✅

### `methods/` Package (Keep as is - well organized)
- `upgma.py` (91 lines) - UPGMA algorithm ✅
- `bionj.py` (~200 lines) - BioNJ algorithm ✅
- `protein_ml.py` (429 lines) - Protein ML ✅

### `models/` Package (Refactor - see Phase 2)
- `ml_tree.py` (322 lines) - Level 1 ⚠️
- `ml_tree_level2.py` (443 lines) - Level 2 ⚠️
- `ml_tree_level3.py` (593 lines) - Level 3 ⚠️
- `ml_tree_level4.py` (409 lines) - Level 4 ✅
- `model_selection.py` (491 lines) - Model selection ⚠️
- `substitution_models.py` (442 lines) - Models ⚠️
- `rate_matrices.py` (291 lines) - Rate matrices ⚠️
- `branch_length_optimizer.py` (260 lines) - Optimization ⚠️
- `tree_search.py` (420 lines) - NNI search ⚠️
- `numba_likelihood.py` (368 lines) - Numba acceleration ⚠️
- `gpu_likelihood_torch.py` (538 lines) - GPU acceleration ⚠️
- `protein_models.py` (356 lines) - Protein models ✅

### `utils/` Package (Refactor - see Phase 1)
- `strain_handler.py` (512 lines) - MOVE to preprocessing/ 🔴
- `dataset_analyzer.py` (294 lines) - MOVE to analysis/ 🔴
- `bootstrap.py` (398 lines) - MOVE to analysis/ 🔴
- `sampling_strategy.py` (350 lines) - MOVE to preprocessing/ 🔴
- `outgroup_handler.py` (~150 lines) - MOVE to preprocessing/ 🔴
- `visualize_trees.py` (~150 lines) - MOVE to visualization/ 🔴
- `console.py` - Keep (actual utility) ✅

### `consensus/` Package (Keep as is or minor consolidation)
- `bipartitions.py` - Bipartition extraction ✅
- `tree_distance.py` - Robinson-Foulds distance ⚠️ (could merge)

### `visualization/` Package (Add moved file)
- `ete3_viz.py` - ETE3 visualization ✅
- ADD: `ascii_viz.py` (from utils/visualize_trees.py) 🟢

### Top-level
- `cli.py` (686 lines) - Command-line interface ✅
- `config.py` - Configuration (minimal) ⚠️ (could expand)

**Legend**:
- ✅ Good as is
- ⚠️ Consider refactoring
- 🔴 Move required
- 🟢 Add/create

---

## Appendix C: Import Path Changes Reference

### After Phase 1 Refactoring

**Old Imports** → **New Imports**

```python
# Preprocessing modules
from rrna_phylo.utils.strain_handler import ...
→ from rrna_phylo.preprocessing.strain_handler import ...

from rrna_phylo.utils.outgroup_handler import ...
→ from rrna_phylo.preprocessing.outgroup_handler import ...

from rrna_phylo.utils.sampling_strategy import ...
→ from rrna_phylo.preprocessing.sampling_strategy import ...

# Analysis modules
from rrna_phylo.utils.bootstrap import ...
→ from rrna_phylo.analysis.bootstrap import ...

from rrna_phylo.utils.dataset_analyzer import ...
→ from rrna_phylo.analysis.dataset_analyzer import ...

# Visualization
from rrna_phylo.utils.visualize_trees import ...
→ from rrna_phylo.visualization.ascii_viz import ...

# Console utility (stays in utils or moves to io)
from rrna_phylo.utils.console import ...
→ from rrna_phylo.io.console import ...  # if moved
```

### Files Requiring Import Updates

**After Phase 1**:
- `cli.py` - Uses strain_handler, bootstrap, sampling, outgroup
- `core/builder_smart.py` - Uses dataset_analyzer
- `utils/__init__.py` - Re-exports moved modules

**After Phase 2** (if subpackaging models):
- All files importing from `models/ml_tree*`
- All files importing from `models/model_selection`
- Test files for ML functionality

---

## Appendix D: Documentation Updates Required

### README.md Updates

**Section: "File Organization"** - Update to reflect new structure:

```markdown
## File Organization

### Core Modules
- `core/` - Tree structures, builders, type detection
  - `tree.py` - TreeNode class
  - `sequence_type.py` - DNA/RNA/Protein detection
  - `builder.py` - Phylogenetic tree builder

### Input/Output
- `io/` - FASTA parsing, alignment
- `visualization/` - Tree visualization (ASCII, ETE3)

### Preprocessing
- `preprocessing/` - Data preparation
  - `strain_handler.py` - Duplicate handling, deduplication
  - `sampling_strategy.py` - Dataset sampling
  - `outgroup_handler.py` - Outgroup management

### Tree Building Methods
- `methods/` - Tree building algorithms
  - `upgma.py` - UPGMA (molecular clock)
  - `bionj.py` - BioNJ (variance-weighted)
  - `protein_ml.py` - Protein ML
- `distance/` - Distance calculations
  - `distance.py` - Jukes-Cantor for DNA/RNA
  - `protein_distance.py` - Protein distances

### Maximum Likelihood
- `models/` - ML models and optimization
  - `ml_tree_level4.py` - Production ML (use this!)
  - `model_selection.py` - AIC/BIC model selection
  - `tree_search.py` - NNI topology search
  - `gpu_likelihood.py` - GPU acceleration

### Analysis
- `analysis/` - Statistical analysis
  - `bootstrap.py` - Bootstrap support
  - `dataset_analyzer.py` - Dataset quality analysis
- `consensus/` - Tree comparison
  - `bipartitions.py` - Split extraction
  - `tree_distance.py` - Robinson-Foulds distance
```

### New Package READMEs

Create `README.md` in each new package:

**`preprocessing/README.md`**:
```markdown
# Preprocessing Module

Handles data preparation before tree building:

- **Strain handling**: Manage multiple rRNA copies from same genome
- **Deduplication**: Remove exact and near-exact duplicates (99.5% threshold)
- **Sampling**: Reduce dataset size while maintaining diversity
- **Outgroup selection**: Identify and manage outgroup sequences

See documentation for detailed usage.
```

**`analysis/README.md`**:
```markdown
# Analysis Module

Statistical analysis tools for phylogenetic trees:

- **Bootstrap**: Non-parametric bootstrap for branch support
- **Dataset analysis**: Assess dataset quality and method suitability

See documentation for detailed usage.
```

---

## Questions for Review

Before implementing this plan, consider:

1. **ML Levels**: Do we actively use Levels 1-3 for testing/education? If yes, keep; if no, archive.

2. **Builder Merge**: Should we merge the two builders now or wait? Merging adds risk but simplifies interface.

3. **Subpackaging**: Is the `models/` package large enough to warrant subpackages now, or wait until it grows more?

4. **Backward Compatibility**: How long should we maintain deprecated import paths? (Recommendation: 1-2 versions)

5. **Testing Coverage**: What's the current test coverage? Should we add tests before refactoring?

---

## Conclusion

This refactoring plan focuses on **organizational improvements** with minimal risk. The codebase is already well-written with good documentation; we're primarily moving files to more logical locations and clarifying architecture.

**Recommended approach**:
1. Start with Phase 1 (package reorganization) - highest impact, lowest risk
2. Gather feedback before proceeding to Phases 2-3
3. Phase 4 (polish) can be done incrementally over time

**Total estimated effort**: 20-30 hours spread over 3-4 weeks.

**Expected outcome**: More intuitive package structure, easier onboarding for new developers, clearer separation of concerns, maintained backward compatibility.

---

**End of Refactoring Plan**

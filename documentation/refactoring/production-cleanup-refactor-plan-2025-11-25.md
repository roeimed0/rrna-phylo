# rRNA-Phylo Production Cleanup Refactoring Plan

**Date**: 2025-11-25
**Branch**: refactor/package-structure
**Status**: Ready for execution
**Risk Level**: Low-Medium (mostly file deletions and cleanup)

---

## Executive Summary

The rRNA-Phylo project has completed rapid development with a fully functional phylogenetic tree building pipeline. The codebase now requires production-ready cleanup to remove temporary files, consolidate duplicate code, and improve maintainability.

**Key Findings**:
- 3 temporary test files in backend root (should move to tests/)
- 12 documentation files in backend/ (should consolidate/reorganize)
- Multiple ML tree implementations (Level 1-4) with duplication
- Broken consensus module (disabled but not removed)
- Unused visualization code (ETE3 vs ASCII)
- Multiple builder classes with overlapping functionality
- Python cache directories (__pycache__) throughout

**Scope**: Aggressive cleanup targeting production-ready state
**Timeline**: 2-3 hours for full execution
**Testing**: Run existing test suite after each phase

---

## Current State Analysis

### Directory Structure

```
backend/
├── rrna_phylo/              # Main package
│   ├── consensus/           # BROKEN - needs removal/fix
│   ├── core/                # Core functionality
│   ├── distance/            # Distance calculations
│   ├── io/                  # I/O operations
│   ├── methods/             # Tree building methods
│   ├── models/              # ML models (4 levels!)
│   ├── utils/               # Utilities
│   └── visualization/       # ETE3 visualization
├── tests/                   # Proper test directory
├── examples/                # Example scripts
├── test_*.py                # TEMPORARY tests in root
├── *.md (12 files)          # Documentation files
└── __pycache__/             # Cache directories
```

### Code Statistics

**Total Python Files**: 41
**Test Files**: 14 (11 in tests/, 3 in root)
**Documentation Files**: 12 .md files
**Builder Classes**: 4 (PhylogeneticTreeBuilder, UPGMABuilder, BioNJBuilder, MLTreeBuilder)
**ML Implementations**: 4 levels (ml_tree.py, ml_tree_level2.py, ml_tree_level3.py, ml_tree_level4.py)

---

## Identified Issues and Opportunities

### CRITICAL Issues

#### 1. Broken Consensus Module (SEVERITY: CRITICAL)
**Location**: `backend/rrna_phylo/consensus/`

**Problem**:
- Algorithm creates invalid bipartitions not present in input trees
- Currently disabled in `builder.py` but still imported
- Taking up package space and mental overhead

**Evidence**:
```python
# From SUMMARY.md:
# Example: Input trees all show {RecA_Ecoli, RecA_Salm} as sisters
# Consensus incorrectly created {RecA_Salm, RecA_Bsub} (not in any tree!)
```

**Files**:
- `rrna_phylo/consensus/__init__.py`
- `rrna_phylo/consensus/bipartitions.py`
- `rrna_phylo/consensus/consensus.py`
- `rrna_phylo/consensus/tree_distance.py`

**Decision**: Remove entire consensus package for production
**Rationale**: Broken code should not ship. Can restore from git when fixed.

---

#### 2. Duplicate ML Tree Implementations (SEVERITY: MAJOR)
**Location**: `backend/rrna_phylo/models/`

**Problem**: 4 ML implementations with increasing sophistication but duplication

**Files**:
1. `ml_tree.py` - Basic GTR model (323 lines)
   - Status: OBSOLETE - placeholder implementation
   - Used by: None (Level 2+ use GTRModel from here)

2. `ml_tree_level2.py` - Full Felsenstein + optimization (444 lines)
   - Status: OBSOLETE - superseded by Level 3
   - Used by: None

3. `ml_tree_level3.py` - GTR+Gamma + site compression (600+ lines)
   - Status: ACTIVE - used by Level 4
   - Used by: `ml_tree_level4.py`, `protein_ml.py`, `visualize_trees.py`

4. `ml_tree_level4.py` - Model selection + tree search (200+ lines)
   - Status: ACTIVE - current production implementation
   - Used by: `builder.py`, `cli.py`

**Consolidation Strategy**:
- Keep: `ml_tree_level3.py`, `ml_tree_level4.py`
- Remove: `ml_tree.py`, `ml_tree_level2.py`
- Extract: GTRModel to `substitution_models.py` (consolidate with existing)

---

#### 3. Temporary Test Files in Root (SEVERITY: MAJOR)
**Location**: `backend/`

**Problem**: Test files outside proper test directory

**Files to Move/Remove**:
1. `test_bias_handling.py` - Demo script, not a test
2. `test_complete_pipeline.py` - Integration test
3. `test_visualization.py` - ETE3 visualization test

**Decision**:
- `test_complete_pipeline.py` → Move to `tests/integration/test_complete_pipeline.py`
- `test_bias_handling.py` → Move to `examples/demo_bias_handling.py` (it's a demo)
- `test_visualization.py` → Move to `examples/demo_ete3_visualization.py` (requires ETE3)

---

### MAJOR Issues

#### 4. Documentation Sprawl (SEVERITY: MAJOR)
**Location**: `backend/`

**Problem**: 12 markdown files in backend root, unclear organization

**Files**:
```
API-USAGE.md              (4.9K)
ARCHITECTURE.md           (12K)
AUTO_BIAS_DETECTION.md    (12K)
BOOTSTRAP_STATUS.md       (8.1K)
CLI_USAGE.md              (8.6K)
DATABASE_BIAS.md          (8.7K)
PERFORMANCE.md            (5.8K)
PROJECT_STATUS.md         (15K)
README.md                 (9.1K)
TESTING_REPORT.md         (21K)
USAGE_GUIDE.md            (13K)
VISUALIZATION_STATUS.md   (12K)
```

**Consolidation Strategy**:

Create organized documentation structure:
```
backend/
├── README.md                    # Keep - main entry point
└── docs/
    ├── user-guide/
    │   ├── cli-usage.md         # CLI_USAGE.md
    │   ├── api-usage.md         # API-USAGE.md
    │   └── usage-guide.md       # USAGE_GUIDE.md
    ├── development/
    │   ├── architecture.md      # ARCHITECTURE.md
    │   ├── testing-report.md    # TESTING_REPORT.md
    │   └── project-status.md    # PROJECT_STATUS.md
    └── features/
        ├── bootstrap.md         # BOOTSTRAP_STATUS.md
        ├── bias-detection.md    # Merge AUTO_BIAS_DETECTION.md + DATABASE_BIAS.md
        ├── performance.md       # PERFORMANCE.md
        └── visualization.md     # VISUALIZATION_STATUS.md
```

---

#### 5. Unused Visualization Code (SEVERITY: MAJOR)
**Location**: `backend/rrna_phylo/visualization/`

**Problem**: ETE3 visualization exists but project uses ASCII only

**Files**:
- `visualization/ete3_viz.py` (197 lines)
- `visualization/__init__.py`

**Evidence from docs**:
```
Tree Visualization - ASCII Only
Decision: Use ASCII-based tree visualization only
Rationale:
- Simple, works everywhere without dependencies
- No external libraries required
- Use external tools (FigTree, iTOL) for publication figures
```

**Decision**: Keep but mark as optional
**Rationale**:
- ETE3 is valuable for users who have it
- ASCII is the default/recommended
- Document clearly that ETE3 is optional

**Action**: Update docs to clarify ETE3 is optional enhancement

---

#### 6. Multiple Builder Classes (SEVERITY: MINOR)
**Location**: Various

**Classes**:
1. `PhylogeneticTreeBuilder` (core/builder.py) - Main unified builder
2. `UPGMABuilder` (methods/upgma.py) - Legacy, still has print_tree method
3. `BioNJBuilder` (methods/bionj.py) - Legacy, still has print_tree method
4. `MLTreeBuilder` (models/ml_tree_level2.py) - OBSOLETE

**Problem**: Duplication of tree printing logic

**Decision**: Consolidate tree printing in utils
**Action**: Remove print_tree methods from UPGMABuilder/BioNJBuilder (use utils.print_tree_ascii)

---

### MINOR Issues

#### 7. Python Cache Directories (SEVERITY: MINOR)
**Location**: Throughout

**Problem**: __pycache__ directories committed or present

**Action**:
- Verify .gitignore has `__pycache__/`
- Clean existing cache: `find . -name __pycache__ -type d -exec rm -rf {} +`

---

#### 8. Unused Imports (SEVERITY: MINOR)

**Files with potential unused imports** (requires detailed analysis):
- `utils/visualize_trees.py` - imports UPGMABuilder, BioNJBuilder classes (only uses functions)
- `consensus/__init__.py` - imports from broken modules
- Various __init__.py files may have outdated imports

**Action**: Run import analysis after file cleanup

---

#### 9. Example Scripts Organization (SEVERITY: MINOR)
**Location**: `backend/examples/`

**Current**: Only `compare_trees.py`
**Needed**: Move demo scripts here

**Action**: Consolidate demo/example scripts in examples/

---

## Proposed Refactoring Plan

### Phase 1: Critical Cleanup (30 min)

**Priority**: Remove broken/obsolete code

#### 1.1 Remove Broken Consensus Module

**Risk**: LOW - already disabled in builder.py

```bash
# Files to delete:
rm -rf backend/rrna_phylo/consensus/
```

**Update imports in**:
- `backend/rrna_phylo/core/builder.py` - Remove consensus import (line 22)
- `backend/rrna_phylo/__init__.py` - Remove consensus exports if any

**Test**:
```bash
cd backend && python -m pytest tests/
```

---

#### 1.2 Remove Obsolete ML Implementations

**Risk**: MEDIUM - ensure Level 3/4 don't depend on Level 1/2

**Before removal, verify imports**:
```bash
grep -r "from.*ml_tree import\|from.*ml_tree_level2 import" backend/rrna_phylo/
```

**Files to delete**:
```bash
rm backend/rrna_phylo/models/ml_tree.py
rm backend/rrna_phylo/models/ml_tree_level2.py
```

**Exception**: Keep GTRModel class - extract first:

1. Move GTRModel from `ml_tree.py` to `substitution_models.py`
2. Update imports in:
   - `ml_tree_level3.py`
   - Any other files importing GTRModel

**Test**:
```bash
cd backend && python -m pytest tests/test_ml*.py
```

---

#### 1.3 Clean Python Cache

```bash
cd backend
find . -name "__pycache__" -type d -exec rm -rf {} +
find . -name "*.pyc" -delete
```

**Update .gitignore**:
```
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
```

---

### Phase 2: Test Organization (20 min)

**Priority**: Organize test files properly

#### 2.1 Move Temporary Tests

**Create structure**:
```bash
mkdir -p backend/tests/integration
mkdir -p backend/examples
```

**Move files**:
```bash
# Integration test
mv backend/test_complete_pipeline.py backend/tests/integration/test_complete_pipeline.py

# Demo scripts (not tests)
mv backend/test_bias_handling.py backend/examples/demo_bias_handling.py
mv backend/test_visualization.py backend/examples/demo_ete3_visualization.py
```

**Update shebang/imports if needed**:
- Ensure relative imports still work
- Update any hardcoded paths

**Test**:
```bash
cd backend && python -m pytest tests/
cd backend && python examples/demo_bias_handling.py  # Manual verification
```

---

### Phase 3: Documentation Consolidation (30 min)

**Priority**: Organize documentation

#### 3.1 Create Documentation Structure

```bash
cd backend
mkdir -p docs/user-guide
mkdir -p docs/development
mkdir -p docs/features
```

#### 3.2 Move and Consolidate Files

**User Guide**:
```bash
mv CLI_USAGE.md docs/user-guide/cli-usage.md
mv API-USAGE.md docs/user-guide/api-usage.md
mv USAGE_GUIDE.md docs/user-guide/usage-guide.md
```

**Development**:
```bash
mv ARCHITECTURE.md docs/development/architecture.md
mv TESTING_REPORT.md docs/development/testing-report.md
mv PROJECT_STATUS.md docs/development/project-status.md
```

**Features**:
```bash
mv BOOTSTRAP_STATUS.md docs/features/bootstrap.md
mv PERFORMANCE.md docs/features/performance.md
mv VISUALIZATION_STATUS.md docs/features/visualization.md
```

**Merge bias documentation**:
```bash
# Combine AUTO_BIAS_DETECTION.md + DATABASE_BIAS.md → docs/features/bias-detection.md
cat AUTO_BIAS_DETECTION.md DATABASE_BIAS.md > docs/features/bias-detection.md
rm AUTO_BIAS_DETECTION.md DATABASE_BIAS.md
```

#### 3.3 Update README.md

Create new backend/README.md with clear structure:

```markdown
# rRNA-Phylo Backend

## Quick Start
See [User Guide](docs/user-guide/)

## Documentation
- **User Guide**: [CLI Usage](docs/user-guide/cli-usage.md), [API Usage](docs/user-guide/api-usage.md)
- **Development**: [Architecture](docs/development/architecture.md), [Testing](docs/development/testing-report.md)
- **Features**: [Bootstrap](docs/features/bootstrap.md), [Bias Detection](docs/features/bias-detection.md)
```

---

### Phase 4: Code Consolidation (40 min)

**Priority**: Remove code duplication

#### 4.1 Consolidate GTRModel

**Current**: GTRModel in `models/ml_tree.py` (to be deleted)
**Target**: Merge into `models/substitution_models.py`

**Action**:
1. Copy GTRModel class from ml_tree.py
2. Add to substitution_models.py (after line 200+)
3. Update imports:
   - `models/ml_tree_level3.py`: `from rrna_phylo.models.substitution_models import GTRModel`
   - Any other files

**Test**:
```bash
cd backend && python -m pytest tests/test_ml*.py
```

---

#### 4.2 Remove Duplicate Tree Printing

**Current**: print_tree methods in UPGMABuilder, BioNJBuilder
**Target**: Use `utils.print_tree_ascii` only

**Action**:
1. Remove `print_tree` method from `methods/upgma.py` (UPGMABuilder)
2. Remove `print_tree` method from `methods/bionj.py` (BioNJBuilder)
3. Update `examples/compare_trees.py`:
   ```python
   from rrna_phylo.utils import print_tree_ascii
   # Replace upgma_builder.print_tree(tree) with print_tree_ascii(tree)
   ```

**Test**:
```bash
cd backend && python examples/compare_trees.py
```

---

#### 4.3 Clean Up Imports in __init__.py Files

**Target files**:
1. `models/__init__.py` - Remove ml_tree_level2 exports
2. `rrna_phylo/__init__.py` - Remove consensus exports

**Action**: Update each __init__.py to only export active code

**Test**: Import check
```python
python -c "from rrna_phylo import *; print('OK')"
```

---

### Phase 5: Final Cleanup (20 min)

#### 5.1 Update Package Metadata

**Files to update**:
1. `backend/rrna_phylo/__init__.py` - Update version/exports
2. Root `README.md` - Update project status

#### 5.2 Run Full Test Suite

```bash
cd backend
python -m pytest tests/ -v
```

#### 5.3 Update SUMMARY.md

Remove references to deleted files:
- Consensus module
- ML Level 1/2
- Old test files

---

## Risk Assessment and Mitigation

### High Risk Changes

#### Removing ML Level 1/2
**Risk**: Breaking imports in unknown locations
**Mitigation**:
- Grep entire codebase before deletion
- Run full test suite immediately after
- Keep in git history (easy rollback)

**Rollback**: `git checkout HEAD -- backend/rrna_phylo/models/ml_tree*.py`

---

#### Moving Test Files
**Risk**: Breaking CI/CD if paths are hardcoded
**Mitigation**:
- Check for any CI config files
- Update test discovery paths
- Manual verification after move

**Rollback**: `git checkout HEAD -- backend/test_*.py`

---

### Medium Risk Changes

#### Removing Consensus Module
**Risk**: Import errors if used elsewhere
**Mitigation**:
- Already disabled in builder.py
- Grep for imports before removal
- Test suite will catch issues

**Rollback**: `git checkout HEAD -- backend/rrna_phylo/consensus/`

---

#### Documentation Reorganization
**Risk**: Breaking external links/references
**Mitigation**:
- Low risk - documentation is internal
- Can add redirects/notes if needed

**Rollback**: Simple git restore

---

### Low Risk Changes

- Python cache cleanup: Zero risk
- Example file moves: Low risk (non-critical)
- Import cleanup: Caught by tests

---

## Testing Strategy

### After Each Phase

**Unit Tests**:
```bash
cd backend
python -m pytest tests/ -v
```

**Import Validation**:
```python
python -c "from rrna_phylo import *"
python -c "from rrna_phylo.models import build_ml_tree_level3"
python -c "from rrna_phylo.methods import build_upgma_tree, build_bionj_tree"
```

**Integration Test**:
```bash
cd backend
python -m pytest tests/integration/test_complete_pipeline.py -v
```

### Final Validation

**Full Pipeline Test**:
```bash
cd backend
python -m rrna_phylo.cli test_real_rrana.fasta --method all -o test_output/
```

**Example Scripts**:
```bash
python examples/compare_trees.py
python examples/demo_bias_handling.py
```

---

## Success Metrics

### Quantitative Metrics

**Before**:
- Python files: 41
- Test files in root: 3
- Documentation files: 12 in backend/
- ML implementations: 4 levels
- Builder classes: 4
- Lines of code: ~8,000

**After**:
- Python files: 37 (-4, -10%)
- Test files in root: 0 (-3, clean)
- Documentation files: 1 README + organized docs/
- ML implementations: 2 levels (Level 3 + 4)
- Builder classes: 3 (remove MLTreeBuilder)
- Lines of code: ~7,200 (-800, -10%)

### Qualitative Metrics

- Clear separation: tests/ vs examples/ vs source
- Organized documentation structure
- Single source of truth for each feature
- Production-ready (no broken/temp code)
- Easier onboarding (clear README → docs/)

---

## Post-Refactoring Checklist

- [ ] All tests passing
- [ ] No broken imports
- [ ] Documentation updated
- [ ] README.md reflects new structure
- [ ] Git history clean (meaningful commits)
- [ ] No __pycache__ directories
- [ ] Examples run successfully
- [ ] CLI works as expected
- [ ] Code review completed

---

## Rollback Plan

### If Issues Arise

**Step 1**: Identify which phase caused the issue

**Step 2**: Rollback that phase
```bash
# Example: Rollback Phase 1.2 (ML cleanup)
git checkout HEAD -- backend/rrna_phylo/models/ml_tree.py
git checkout HEAD -- backend/rrna_phylo/models/ml_tree_level2.py
```

**Step 3**: Run tests to verify rollback
```bash
cd backend && python -m pytest tests/ -v
```

**Step 4**: Debug issue, fix, retry

### Complete Rollback

If all else fails:
```bash
git reset --hard HEAD
```

---

## Future Recommendations

### After This Refactoring

1. **Add Import Linter**
   - Use `pylint` or `flake8` to catch unused imports
   - Integrate into pre-commit hooks

2. **Add Documentation Generator**
   - Use Sphinx or MkDocs for API documentation
   - Auto-generate from docstrings

3. **Add Type Checking**
   - Use `mypy` for static type checking
   - Add type hints to all public APIs

4. **Add Code Coverage**
   - Use `pytest-cov` to track test coverage
   - Target >80% coverage

5. **Fix or Remove Consensus**
   - Either: Implement correct bipartition algorithm
   - Or: Permanently remove from project scope

6. **Consolidate Builder Pattern**
   - Consider single TreeBuilder class with strategy pattern
   - Remove individual builder classes

---

## Execution Timeline

**Total Estimated Time**: 2-3 hours

| Phase | Duration | Dependencies |
|-------|----------|--------------|
| Phase 1: Critical Cleanup | 30 min | None |
| Phase 2: Test Organization | 20 min | Phase 1 |
| Phase 3: Documentation | 30 min | None (parallel with Phase 2) |
| Phase 4: Code Consolidation | 40 min | Phase 1 |
| Phase 5: Final Cleanup | 20 min | All previous |
| **Buffer** | 20 min | Testing/fixes |

**Recommended Approach**: Execute phases sequentially with git commits after each

---

## Git Commit Strategy

### Recommended Commits

```bash
# Phase 1
git commit -m "refactor: remove broken consensus module"
git commit -m "refactor: remove obsolete ML Level 1 and 2 implementations"
git commit -m "chore: clean Python cache directories"

# Phase 2
git commit -m "refactor: reorganize test files into proper structure"
git commit -m "refactor: move demo scripts to examples/"

# Phase 3
git commit -m "docs: reorganize documentation into structured hierarchy"

# Phase 4
git commit -m "refactor: consolidate GTRModel into substitution_models"
git commit -m "refactor: remove duplicate tree printing methods"
git commit -m "refactor: clean up package imports"

# Phase 5
git commit -m "docs: update README and project metadata"
git commit -m "chore: production cleanup complete"
```

---

## Appendix A: File Deletion Summary

### Files to Delete

**Consensus Module** (4 files):
```
backend/rrna_phylo/consensus/__init__.py
backend/rrna_phylo/consensus/bipartitions.py
backend/rrna_phylo/consensus/consensus.py
backend/rrna_phylo/consensus/tree_distance.py
```

**Obsolete ML** (2 files):
```
backend/rrna_phylo/models/ml_tree.py
backend/rrna_phylo/models/ml_tree_level2.py
```

**Documentation** (move, not delete):
```
backend/CLI_USAGE.md → docs/user-guide/cli-usage.md
backend/API-USAGE.md → docs/user-guide/api-usage.md
backend/USAGE_GUIDE.md → docs/user-guide/usage-guide.md
backend/ARCHITECTURE.md → docs/development/architecture.md
backend/TESTING_REPORT.md → docs/development/testing-report.md
backend/PROJECT_STATUS.md → docs/development/project-status.md
backend/BOOTSTRAP_STATUS.md → docs/features/bootstrap.md
backend/PERFORMANCE.md → docs/features/performance.md
backend/VISUALIZATION_STATUS.md → docs/features/visualization.md
backend/AUTO_BIAS_DETECTION.md + DATABASE_BIAS.md → docs/features/bias-detection.md (merge)
```

**Temporary Tests** (move, not delete):
```
backend/test_complete_pipeline.py → tests/integration/test_complete_pipeline.py
backend/test_bias_handling.py → examples/demo_bias_handling.py
backend/test_visualization.py → examples/demo_ete3_visualization.py
```

**Total**: 6 deletions, 13 moves, 2 merges

---

## Appendix B: Import Update Checklist

### Files Requiring Import Updates

After removing consensus module:
- [ ] `backend/rrna_phylo/core/builder.py` (line 22)
- [ ] `backend/rrna_phylo/__init__.py` (check for consensus exports)

After removing ML Level 1/2:
- [ ] `backend/rrna_phylo/models/__init__.py` (remove exports)

After consolidating GTRModel:
- [ ] `backend/rrna_phylo/models/ml_tree_level3.py` (update import)
- [ ] Any files importing from ml_tree.py

After moving test files:
- [ ] `backend/tests/integration/test_complete_pipeline.py` (verify imports)
- [ ] `backend/examples/demo_bias_handling.py` (verify imports)
- [ ] `backend/examples/demo_ete3_visualization.py` (verify imports)

---

## Appendix C: Quick Reference Commands

### Search for imports of removed modules
```bash
cd backend
grep -r "from.*consensus import\|import.*consensus" rrna_phylo/
grep -r "from.*ml_tree import\|import.*ml_tree" rrna_phylo/
grep -r "from.*ml_tree_level2 import" rrna_phylo/
```

### Find all __init__.py files
```bash
find backend/rrna_phylo -name "__init__.py"
```

### Count Python files
```bash
find backend/rrna_phylo -name "*.py" | wc -l
```

### Run specific test suites
```bash
# ML tests
python -m pytest tests/test_ml*.py -v

# Integration tests
python -m pytest tests/integration/ -v

# All tests
python -m pytest tests/ -v
```

### Verify imports
```bash
python -c "from rrna_phylo import *; print('rrna_phylo: OK')"
python -c "from rrna_phylo.models import build_ml_tree_level3, build_ml_tree_level4; print('models: OK')"
python -c "from rrna_phylo.methods import build_upgma_tree, build_bionj_tree; print('methods: OK')"
```

---

**END OF REFACTORING PLAN**

This plan is ready for execution. Proceed phase-by-phase with testing after each step.

# rRNA-Phylo Refactoring Summary

**Quick Reference Guide**

---

## TL;DR

The codebase is **well-structured** overall, but some modules are in the wrong packages. Main issue: `utils/` contains modules that aren't utilities.

**Impact**: Medium effort (~20-30 hours), high value, low risk.

---

## Main Issues & Solutions

### 1. Misplaced Modules in `utils/` (HIGH PRIORITY)

**Problem**: 6 files (~2,200 lines) in `utils/` that aren't utilities.

**Solution**: Create proper packages

```
utils/strain_handler.py       → preprocessing/strain_handler.py      (512 lines)
utils/sampling_strategy.py    → preprocessing/sampling_strategy.py   (350 lines)
utils/outgroup_handler.py     → preprocessing/outgroup_handler.py    (~150 lines)
utils/bootstrap.py            → analysis/bootstrap.py                (398 lines)
utils/dataset_analyzer.py     → analysis/dataset_analyzer.py         (294 lines)
utils/visualize_trees.py      → visualization/ascii_viz.py           (~150 lines)
utils/console.py              → Keep (actual utility) ✅
```

**Effort**: 4-8 hours
**Risk**: LOW (just file moves)

---

### 2. `models/` Package Organization (MEDIUM PRIORITY)

**Problem**: 13 files (5,000+ lines) with mixed concerns.

**Solution**: Either:
- **Option A**: Keep flat, but rename for clarity
  - `ml_tree_level3.py` → `ml_tree.py` (production)
  - `ml_tree_level4.py` → `ml_tree_full.py`
  - Levels 1-2 → `reference/` subdirectory

- **Option B**: Create subpackages (only if continuing to grow)
  ```
  models/
  ├── ml/              # Core ML implementations
  ├── selection/       # Model selection
  ├── optimization/    # Branch length, tree search
  └── acceleration/    # GPU, Numba
  ```

**Recommendation**: Start with Option A (simpler).

**Effort**: 6-10 hours
**Risk**: MEDIUM (import updates)

---

### 3. Duplicate Builder Classes (MEDIUM PRIORITY)

**Problem**: Two builders with overlapping functionality:
- `PhylogeneticTreeBuilder` (451 lines)
- `SmartPhylogeneticTreeBuilder` (258 lines)

**Solution**:
- **Short term**: Document when to use each
- **Long term**: Consider merge with `smart_mode` flag

**Effort**: 4-6 hours (docs) or 8-12 hours (merge)
**Risk**: LOW (docs) or MEDIUM (merge)

---

## Quick Wins (Do First)

1. **Move files from `utils/` to proper packages** (Phase 1)
   - Creates `preprocessing/` and `analysis/` packages
   - Immediate clarity improvement
   - 4-8 hours, LOW risk

2. **Add builder documentation**
   - Clear examples for when to use each
   - 2-3 hours, VERY LOW risk

3. **Rename `gpu_likelihood_torch.py` → `gpu_likelihood.py`**
   - Remove implementation detail from name
   - 1 hour, VERY LOW risk

---

## What NOT to Change

These are **already well-organized**:

- ✅ `core/` - Tree structures, builders, type detection
- ✅ `io/` - FASTA parsing, alignment
- ✅ `distance/` - Distance calculations
- ✅ `methods/` - UPGMA, BioNJ, protein ML
- ✅ `consensus/` - Tree comparison, bipartitions
- ✅ `visualization/` - ETE3 visualization
- ✅ `cli.py` - Command-line interface (686 lines is reasonable)

---

## Implementation Roadmap

### Week 1: Package Reorganization ⭐ START HERE
- Create `preprocessing/` and `analysis/` packages
- Move 6 files from `utils/`
- Update imports in `cli.py` and affected files
- Run tests to verify

**Deliverable**: Clearer package structure

---

### Week 2: Documentation
- Add comprehensive docstrings to builders
- Update README with new structure
- Create package-level READMEs

**Deliverable**: Better documentation

---

### Week 3: Models Cleanup (Optional)
- Archive ML reference implementations (Levels 1-3)
- Rename GPU module
- Add ML progression README

**Deliverable**: Clearer ML code

---

### Week 4: Polish (Optional)
- Review verbose docstrings
- Update architecture docs
- Final testing

**Deliverable**: Polished codebase

---

## File Move Checklist

For each file move:

1. [ ] Create destination directory (if new)
2. [ ] Move file
3. [ ] Update imports in moved file
4. [ ] Find files importing it: `grep -r "from.*{old_path}" .`
5. [ ] Update importing files
6. [ ] Update `__init__.py` files
7. [ ] Run tests
8. [ ] Update docs
9. [ ] Commit

---

## Risk Matrix

| Change | Impact | Effort | Risk | Priority |
|--------|--------|--------|------|----------|
| Create preprocessing/ package | HIGH | 4h | LOW | 🔴 HIGH |
| Create analysis/ package | HIGH | 4h | LOW | 🔴 HIGH |
| Move visualize_trees | MED | 1h | LOW | 🟡 MEDIUM |
| Document builders | MED | 3h | LOW | 🟡 MEDIUM |
| Rename GPU module | LOW | 1h | LOW | 🟢 LOW |
| Archive ML levels | MED | 3h | LOW | 🟢 LOW |
| Subpackage models/ | MED | 8h | MED | 🟢 LOW |
| Merge builders | HIGH | 12h | MED | ⚪ FUTURE |

---

## Before & After

### Current Structure
```
rrna_phylo/
├── core/          # ✅ Good
├── io/            # ✅ Good
├── distance/      # ✅ Good
├── methods/       # ✅ Good
├── models/        # ⚠️ Could be clearer
├── utils/         # 🔴 Misplaced modules
├── consensus/     # ✅ Good
└── visualization/ # ✅ Good
```

### Proposed Structure
```
rrna_phylo/
├── core/            # ✅ Tree structures, builders
├── io/              # ✅ FASTA, alignment
├── distance/        # ✅ Distance calculations
├── methods/         # ✅ UPGMA, BioNJ
├── models/          # ✅ ML (better organized)
├── preprocessing/   # 🟢 NEW: Data prep
├── analysis/        # 🟢 NEW: Statistical analysis
├── consensus/       # ✅ Tree comparison
└── visualization/   # ✅ Tree visualization
```

---

## Success Criteria

After refactoring:

- ✅ No non-utility modules in `utils/`
- ✅ Clear package purposes
- ✅ All tests pass
- ✅ No broken imports
- ✅ Documentation updated
- ✅ Backward compatibility maintained

---

## Questions?

See full plan: `architecture-refactor-plan-2025-12-03.md`

**Recommendation**: Start with Week 1 (package reorganization), then reassess.

---

**Total Effort**: 20-30 hours over 3-4 weeks
**Risk Level**: LOW (mostly file moves)
**Impact Level**: HIGH (much clearer organization)

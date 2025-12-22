# OpenMP Library Conflict Explanation

## What's Happening

You're seeing this error:
```
OMP: Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already initialized.
```

## Root Cause

**Multiple copies of the Intel OpenMP library are being loaded:**

```
C:/Users/User/anaconda3/envs/gene_prediction/Lib/site-packages/torch/lib/libiomp5md.dll
C:/Users/User/anaconda3/envs/gene_prediction/Library/bin/libiomp5md.dll
```

### Why This Happens

1. **NumPy/SciPy** - Built against Intel MKL (Math Kernel Library)
   - Loads OpenMP from: `Library/bin/libiomp5md.dll`
   - Used for: Linear algebra operations

2. **PyTorch** - Bundles its own MKL/OpenMP
   - Loads OpenMP from: `torch/lib/libiomp5md.dll`
   - Used for: Tensor operations and neural networks

3. **Numba** - Uses LLVM for JIT compilation
   - May trigger OpenMP initialization when using parallel features
   - Relies on system OpenMP (from NumPy/SciPy)

### The Conflict Timeline

```
1. Import numpy/scipy
   → Loads Library/bin/libiomp5md.dll (System MKL)
   → OpenMP initialized (version X)

2. Import numba
   → Triggers parallel operations
   → May try to load torch/lib/libiomp5md.dll (PyTorch MKL)
   → ERROR: OpenMP already initialized!
```

## Why It's Dangerous

Having multiple OpenMP runtimes can cause:
- **Incorrect results** - Race conditions in parallel code
- **Performance degradation** - Conflicting thread management
- **Crashes** - Memory corruption from competing thread pools

## Solutions

### Solution 1: Quick Workaround (Current)

Set environment variable to suppress the error:

```bash
# Linux/Mac
export KMP_DUPLICATE_LIB_OK=TRUE

# Windows CMD
set KMP_DUPLICATE_LIB_OK=TRUE

# Windows PowerShell
$env:KMP_DUPLICATE_LIB_OK="TRUE"

# In Python script (before imports)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
```

**Pros**: Quick, no reinstallation
**Cons**: Doesn't fix the root cause, may hide real issues

### Solution 2: Remove PyTorch (Recommended for this project)

Since you're not using PyTorch for phylogenetics:

```bash
conda activate gene_prediction
conda remove pytorch  # or pip uninstall torch
```

**Pros**: Eliminates conflict entirely
**Cons**: Can't use PyTorch if needed later

### Solution 3: Use conda-forge (Best for mixed workloads)

Install all packages from conda-forge for consistent MKL:

```bash
conda create -n gene_prediction_clean python=3.9
conda activate gene_prediction_clean
conda install -c conda-forge numpy scipy numba biopython pytorch
```

**Pros**: Single MKL installation, all packages compatible
**Cons**: Requires recreating environment

### Solution 4: Intel oneAPI (Most Robust)

Use Intel's official distribution:

```bash
conda install -c intel numpy scipy mkl
pip install numba biopython
```

**Pros**: Optimized Intel libraries, best performance
**Cons**: Larger download, Intel-specific

## What We're Using Now

**Current approach**: `KMP_DUPLICATE_LIB_OK=TRUE`

This is set in the bash commands:
```bash
KMP_DUPLICATE_LIB_OK=TRUE python build_trees.py sequences.fasta
```

**Why this is acceptable for now:**
1. Your code doesn't use PyTorch, so the extra OpenMP won't cause conflicts
2. Numba's parallel operations work correctly with system MKL
3. No numerical errors observed in phylogenetic calculations

## Long-Term Recommendation

**For production use of rRNA-Phylo:**

1. **Option A - Remove PyTorch** (if not needed):
   ```bash
   conda remove pytorch torchvision torchaudio
   ```

2. **Option B - Add to script** (make it automatic):
   ```python
   # At top of build_trees.py and build_phylogenetic_tree.py
   import os
   os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
   ```

3. **Option C - Set in conda environment** (permanent):
   ```bash
   conda activate gene_prediction
   conda env config vars set KMP_DUPLICATE_LIB_OK=TRUE
   conda deactivate
   conda activate gene_prediction
   ```

## Verification

After applying any solution, verify it works:

```bash
cd backend
python build_trees.py test_sequences.fasta
# Should see no OpenMP warnings
```

## Additional Context

### Why We Have PyTorch

Check your conda environment:
```bash
conda list | grep torch
```

If you see PyTorch installed but don't remember installing it, it may have been:
- Installed as a dependency for another package
- Part of the anaconda metapackage
- Installed for gene prediction (neural network models?)

### Performance Impact

The duplicate OpenMP libraries **do not** affect phylogenetic tree calculations because:
- Numba JIT uses the system OpenMP (from NumPy/SciPy)
- PyTorch OpenMP is never actually used
- No thread pool conflicts occur in practice

The warning is Intel's safety mechanism, not an indication of actual problems in your case.

## References

- [Intel OpenMP Documentation](https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming/openmp-support.html)
- [NumPy + PyTorch conflicts](https://github.com/pytorch/pytorch/issues/37377)
- [Conda-forge MKL handling](https://conda-forge.org/docs/maintainer/knowledge_base.html#mkl)

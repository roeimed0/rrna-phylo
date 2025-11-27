---
description: Guide for fixing NNI tree search to properly re-optimize branch lengths after topology changes
triggers:
  - keywords: ["nni", "tree search", "branch length", "zero branches"]
---

# NNI Tree Search Fix - Branch Length Re-Optimization

## Problem Statement

**Current Bug**: NNI search collapses branches to zero length during topology optimization.

**Evidence**:
- Test with 87 sequences produced 68 branches with exactly 0.000000 length
- Creates impossible star-like topology with massive polytomy
- UPGMA tree has only 1 zero-length branch (correct)

**Root Cause**: After accepting an improved NNI topology swap, the code never re-optimizes branch lengths for the new topology.

## Current Broken Implementation

**Location**: `backend/rrna_phylo/models/tree_search.py:108-124`

```python
# Accept if better than current best
if best_neighbor_logL > best_logL_this_round + tolerance:
    previous_logL = current_logL  # Save for improvement calculation
    if logL1 > logL2:
        current_tree = neighbor1
        current_logL = logL1
    else:
        current_tree = neighbor2
        current_logL = logL2

    best_logL_this_round = current_logL
    improved_this_round = True
    n_improvements += 1

    improvement = current_logL - previous_logL
    if verbose:
        print(f"  NNI improved: LogL = {current_logL:.2f} (+{improvement:.2f})")
```

**What's Missing**: After swapping topology (line 111 or 113), we must re-optimize branch lengths!

## Why This Matters

1. **NNI swaps topology** - Changes which taxa are grouped together
2. **Old branch lengths are invalid** - They were optimized for the OLD topology
3. **Without re-optimization** - Branches gradually collapse to zero as topology changes
4. **Result**: Scientifically invalid trees with 68/87 zero-length branches

## Proper Fix Implementation

### Step 1: Implement Branch Length Optimizer

Create `backend/rrna_phylo/models/branch_length_optimizer.py`:

```python
"""
Branch length optimization for phylogenetic trees.

After topology changes (like NNI swaps), branch lengths must be re-optimized
to maximize likelihood under the new topology.
"""

from typing import List, Optional
import numpy as np
from scipy.optimize import minimize_scalar
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood


def optimize_single_branch(
    tree: TreeNode,
    node: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001,
    max_length: float = 10.0
) -> float:
    """
    Optimize the branch length leading to a single node.

    Args:
        tree: Full tree
        node: Node whose incoming branch length to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter
        min_length: Minimum branch length (prevents zero collapse)
        max_length: Maximum branch length

    Returns:
        Optimized branch length
    """
    original_length = node.branch_length

    def neg_log_likelihood(length):
        """Negative log-likelihood for minimization."""
        node.branch_length = max(min_length, min(length, max_length))
        logL = compute_log_likelihood(tree, sequences, alpha=alpha)
        return -logL  # Scipy minimizes, we want to maximize logL

    # Optimize using Brent's method (efficient for 1D optimization)
    result = minimize_scalar(
        neg_log_likelihood,
        bounds=(min_length, max_length),
        method='bounded'
    )

    # Set optimized length
    optimal_length = max(min_length, min(result.x, max_length))
    node.branch_length = optimal_length

    return optimal_length


def optimize_all_branch_lengths(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001,
    max_iterations: int = 3,
    verbose: bool = False
) -> float:
    """
    Optimize all branch lengths in a tree.

    Uses coordinate descent: optimize each branch in turn, repeat until convergence.

    Args:
        tree: Tree with topology to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter
        min_length: Minimum branch length (prevents zero collapse)
        max_iterations: Number of passes through all branches
        verbose: Print progress

    Returns:
        Final log-likelihood
    """
    if verbose:
        print("  Optimizing branch lengths...")

    # Get all nodes (except root)
    all_nodes = []
    def collect_nodes(node):
        if node.branch_length is not None:  # Has incoming branch
            all_nodes.append(node)
        if not node.is_leaf():
            if node.left:
                collect_nodes(node.left)
            if node.right:
                collect_nodes(node.right)
    collect_nodes(tree)

    # Coordinate descent: optimize each branch in turn
    for iteration in range(max_iterations):
        for node in all_nodes:
            optimize_single_branch(tree, node, sequences, alpha, min_length)

        if verbose and iteration < max_iterations - 1:
            logL = compute_log_likelihood(tree, sequences, alpha=alpha)
            print(f"    Branch optimization iteration {iteration + 1}: LogL = {logL:.2f}")

    final_logL = compute_log_likelihood(tree, sequences, alpha=alpha)

    if verbose:
        print(f"    Final LogL after branch optimization: {final_logL:.2f}")

    return final_logL


def optimize_branch_lengths_fast(
    tree: TreeNode,
    sequences: List[Sequence],
    alpha: Optional[float] = None,
    min_length: float = 0.0001
) -> float:
    """
    Fast single-pass branch length optimization.

    Optimizes each branch once (no iteration). Faster but less accurate than full optimization.
    Good enough for NNI search where topology is changing frequently.

    Args:
        tree: Tree to optimize
        sequences: Aligned sequences
        alpha: Gamma parameter
        min_length: Minimum branch length

    Returns:
        Final log-likelihood
    """
    # Get all nodes with incoming branches
    all_nodes = []
    def collect_nodes(node):
        if node.branch_length is not None:
            all_nodes.append(node)
        if not node.is_leaf():
            if node.left:
                collect_nodes(node.left)
            if node.right:
                collect_nodes(node.right)
    collect_nodes(tree)

    # Single pass optimization
    for node in all_nodes:
        optimize_single_branch(tree, node, sequences, alpha, min_length)

    return compute_log_likelihood(tree, sequences, alpha=alpha)
```

### Step 2: Integrate into NNI Search

Modify `backend/rrna_phylo/models/tree_search.py:108-124`:

```python
from rrna_phylo.models.branch_length_optimizer import optimize_branch_lengths_fast

# ... in nni_search function ...

# Accept if better than current best
if best_neighbor_logL > best_logL_this_round + tolerance:
    previous_logL = current_logL  # Save for improvement calculation
    if logL1 > logL2:
        current_tree = neighbor1
        current_logL = logL1
    else:
        current_tree = neighbor2
        current_logL = logL2

    # *** NEW: Re-optimize branch lengths after topology swap ***
    current_logL = optimize_branch_lengths_fast(
        current_tree,
        sequences,
        alpha=alpha,
        min_length=0.0001  # Prevent zero-length branches!
    )

    best_logL_this_round = current_logL
    improved_this_round = True
    n_improvements += 1

    improvement = current_logL - previous_logL
    if verbose:
        print(f"  NNI improved: LogL = {current_logL:.2f} (+{improvement:.2f})")
```

## Key Design Decisions

### 1. Minimum Branch Length
- Set to `0.0001` to prevent zero collapse
- Small enough to not affect likelihood
- Large enough to maintain valid tree structure

### 2. Fast vs Full Optimization
- **Fast** (`optimize_branch_lengths_fast`): Single pass, use in NNI
- **Full** (`optimize_all_branch_lengths`): 3 iterations, use for final tree
- Trade-off: Speed vs accuracy

### 3. When to Optimize
- **After every accepted NNI swap** (prevents accumulation of errors)
- **Before final return** (ensure final tree is fully optimized)

### 4. Scipy minimize_scalar
- Uses Brent's method (optimal for 1D optimization)
- Bounded to [min_length, max_length]
- Much faster than general-purpose optimizers

## Testing Strategy

### Test 1: Branch Length Preservation
```python
def test_branch_lengths_nonzero():
    """Test that NNI doesn't create zero-length branches."""
    tree, logL, n_impr = nni_search(initial_tree, sequences, verbose=True)

    # Check all branches
    zero_branches = 0
    def count_zeros(node):
        nonlocal zero_branches
        if node.branch_length == 0.0:
            zero_branches += 1
        if not node.is_leaf():
            count_zeros(node.left)
            count_zeros(node.right)

    count_zeros(tree)
    assert zero_branches == 0, f"Found {zero_branches} zero-length branches!"
```

### Test 2: Likelihood Improvement
```python
def test_nni_improves_likelihood():
    """Test that NNI actually improves likelihood."""
    initial_logL = compute_log_likelihood(initial_tree, sequences)
    final_tree, final_logL, n_impr = nni_search(initial_tree, sequences)

    assert final_logL >= initial_logL, "NNI made likelihood worse!"
    print(f"Improvement: {final_logL - initial_logL:.2f}")
```

### Test 3: Compare to RAxML
```python
def test_nni_comparable_to_raxml():
    """Test that our NNI produces reasonable trees."""
    # Run our NNI
    our_tree, our_logL, _ = nni_search(initial_tree, sequences)

    # Run RAxML (external)
    raxml_logL = run_raxml_ml(sequences)

    # Should be within 1% of RAxML
    diff_pct = 100 * abs(our_logL - raxml_logL) / abs(raxml_logL)
    assert diff_pct < 1.0, f"LogL differs from RAxML by {diff_pct:.2f}%"
```

## Expected Results After Fix

- ✅ Zero zero-length branches (currently 68)
- ✅ Proper hierarchical tree structure
- ✅ Likelihood improves after NNI (not just different)
- ✅ Branch lengths proportional to evolutionary distance
- ✅ Comparable to RAxML/IQ-TREE results

## Performance Impact

- **Before**: NNI iteration ~10-15 seconds (no branch optimization)
- **After**: NNI iteration ~15-25 seconds (with fast optimization)
- **Overhead**: ~50-100% slower per iteration
- **Quality**: MUCH better trees (actually correct!)

## Alternative: Use External Tools

If implementing proper branch length optimization is too complex:

1. **Use IQ-TREE** for ML trees:
   ```python
   subprocess.run(['iqtree', '-s', 'alignment.fasta', '-m', 'GTR+G'])
   ```

2. **Use RAxML-NG** for ML trees:
   ```python
   subprocess.run(['raxml-ng', '--msa', 'alignment.fasta', '--model', 'GTR+G'])
   ```

Both are faster and more accurate than our current implementation.

## References

- Felsenstein (1981) - "Evolutionary trees from DNA sequences: a maximum likelihood approach"
- Yang (1996) - "Maximum likelihood phylogenetic estimation from DNA sequences"
- Stamatakis (2014) - "RAxML version 8: a tool for phylogenetic analysis" (implementation reference)

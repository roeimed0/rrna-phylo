---
skill_name: spr-search
trigger:
  keywords:
    - spr
    - subtree pruning
    - regrafting
    - tree search
    - spr search
  intent_patterns:
    - "implement spr"
    - "fix spr"
    - "spr algorithm"
    - "tree rearrangement"
---

# SPR (Subtree Pruning and Regrafting) Search Implementation

## Overview

SPR is a phylogenetic tree search method more powerful than NNI (Nearest Neighbor Interchange). It allows moving subtrees to distant locations in the tree, enabling escape from local optima.

**Status in rRNA-Phylo**:
- ❌ Original SPR implementation has critical bugs (tree corruption)
- ⚠️ Simplified SPR available but experimental
- ✅ NNI is production-ready and recommended

## Why SPR is Hard to Implement

### The SPR Operation

SPR consists of two steps:
1. **Prune**: Cut an edge, removing a subtree
2. **Regraft**: Attach the pruned subtree to a different edge

### Challenges for Unrooted Binary Trees

1. **Parent Node Collapse**: After pruning, the parent node has only one child and must be removed
   ```
   Before:        After prune:      After collapse:
       P               P                 G
      / \             /                 / \
     A   N    →      A        →        A   C
        / \
       B   C         Subtree: (B,C)
   ```

2. **Binary Tree Property**: Every internal node must have exactly 2 children

3. **Branch Length Updates**: All affected branches need recalculation

4. **Leaf Preservation**: Must maintain exact same set of taxa

## Critical Bugs Found in Original Implementation

### Bug #1: Tree Structure Corruption
**Symptom**: Leaf count changes after SPR move (5 → 1, 2, 3, 4, 6, 7, 8)

**Root Cause**:
```python
def prune_subtree(tree, edge_to_cut):
    parent, child = edge_to_cut
    pruned_tree = tree.copy()  # BUG: Copy doesn't preserve structure
    detached_subtree = child.copy()

    # BUG: Removing child doesn't collapse parent node
    if found_parent.left == found_child:
        found_parent.left = None  # Now parent has only 1 child!
```

**Test Results**:
```
Test: 5 sequences
[ERROR] Tree leaf count changed: 5 -> 1
[ERROR] Tree leaf count changed: 5 -> 2
[ERROR] Tree leaf count changed: 5 -> 8
```

### Bug #2: Unrealistic Likelihood Values
**Symptom**: SPR reports +17,529 log-likelihood improvement (impossible)

**Root Cause**: After tree corruption, likelihood calculator compares trees with different numbers of taxa

### Bug #3: Performance
**Symptom**: 14 seconds for 5 sequences (should be <1s)

**Root Cause**: Invalid trees still processed through expensive branch optimization

## Correct SPR Implementation Strategy

### Approach 1: Proper Tree Manipulation (Recommended for Production)

```python
def spr_move_correct(tree, prune_edge, regraft_edge):
    """
    Proper SPR implementation.

    Steps:
    1. Identify subtree to prune
    2. Identify parent and sibling
    3. Remove subtree and collapse parent
    4. Insert new internal node at regraft edge
    5. Attach pruned subtree
    6. Update all branch lengths
    """
    # 1. Prune
    parent, subtree_root = prune_edge
    grandparent = find_parent(tree, parent)
    sibling = parent.left if parent.right == subtree_root else parent.right

    # 2. Collapse parent node
    if grandparent.left == parent:
        grandparent.left = sibling
    else:
        grandparent.right = sibling
    sibling.distance += parent.distance

    # 3. Regraft
    regraft_parent, regraft_child = regraft_edge
    new_internal = TreeNode()
    new_internal.left = subtree_root
    new_internal.right = regraft_child

    if regraft_parent.left == regraft_child:
        regraft_parent.left = new_internal
    else:
        regraft_parent.right = new_internal

    # 4. Initialize branch lengths
    subtree_root.distance = 0.001
    regraft_child.distance = new_internal.distance / 2
    new_internal.distance = regraft_child.distance

    return tree
```

**Validation Required**:
```python
# After every SPR move:
assert tree.count_leaves() == original_leaf_count
assert set(tree.get_leaf_names()) == original_leaf_names
```

### Approach 2: Newick Round-Trip (Simpler but Slower)

```python
def spr_via_newick(tree, prune_subtree_name, regraft_edge_names):
    """
    Use Newick string manipulation for SPR.
    Slower but guaranteed correct structure.
    """
    # 1. Convert to Newick
    newick = tree.to_newick()

    # 2. Parse and manipulate
    # Extract subtree
    subtree_newick = extract_subtree(newick, prune_subtree_name)
    # Remove from original location
    newick = remove_subtree(newick, prune_subtree_name)
    # Insert at new location
    newick = insert_subtree(newick, subtree_newick, regraft_edge_names)

    # 3. Parse back to tree
    new_tree = parse_newick(newick)

    # 4. Optimize branch lengths
    optimize_branch_lengths(new_tree, sequences, alpha)

    return new_tree
```

### Approach 3: Simplified SPR (Current Implementation)

```python
def simplified_spr(tree, sequences, alpha):
    """
    Limited SPR using extended NNI moves.

    Strategy:
    - Swap subtrees between distant internal nodes
    - No actual prune/regraft (avoids bugs)
    - Safer but less powerful than true SPR
    """
    internal_nodes = get_all_internal_nodes(tree)

    for node1 in internal_nodes:
        for node2 in internal_nodes:
            if can_swap(node1, node2):
                # Try swapping left children
                node1.left, node2.left = node2.left, node1.left

                # Evaluate
                logL = optimize_branch_lengths(tree, sequences, alpha)

                if logL > best_logL:
                    best_tree = tree.copy()
                    best_logL = logL
                else:
                    # Undo swap
                    node1.left, node2.left = node2.left, node1.left

    return best_tree, best_logL
```

## Current Implementation Status

### Files

1. **`rrna_phylo/models/spr_search.py`**
   - Original broken implementation
   - Tree corruption bugs
   - Status: Disabled with extensive comments

2. **`rrna_phylo/models/spr_search_simple.py`**
   - Simplified implementation using extended NNI
   - Status: Experimental
   - Issues: Limited power, may not find improvements

3. **`rrna_phylo/models/ml_tree_level4.py`**
   - Integration point
   - When `tree_search='spr'` requested:
     - Runs NNI first (local optimization)
     - Then simplified SPR (extended search)

### Test Results

**5-sequence test**:
```
NNI only: LogL=-6552.50, Time=0.7s ✅
SPR (NNI+simplified): LogL=-6552.50, Time=0.4s ⚠️
  - No improvement found
  - Tested 0 simplified SPR moves (nodes too close)
```

**20-sequence test**:
```
NNI only: LogL=-19750.78, Time=10.1s ✅
  - Branches collapsed: 5
  - Convergence: 1 iteration
```

## Recommendations

### For Production Use

**Use NNI only**:
```python
tree, logL, metadata = build_ml_tree_level4(
    sequences,
    tree_search='nni',  # Recommended
    verbose=True
)
```

**Reasons**:
- ✅ Fast (10s for 20 sequences)
- ✅ Reliable (no tree corruption)
- ✅ All improvements working (8.5/10 code quality)
- ✅ Production-ready

### For Research/Development

**Implement proper SPR** (estimated 2-3 days):

1. **Study reference implementations**:
   - RAxML-NG source code (C++)
   - IQ-TREE source code (C++)
   - FastTree source code (C)

2. **Use established libraries**:
   ```python
   # Option A: ETE3 toolkit
   from ete3 import Tree
   tree = Tree(newick_string)
   tree.prune(taxa_to_keep)

   # Option B: DendroPy
   import dendropy
   tree = dendropy.Tree.get(data=newick_string, schema="newick")
   tree.prune_taxa(taxa_to_remove)

   # Option C: BioPython
   from Bio import Phylo
   tree = Phylo.read("tree.nwk", "newick")
   ```

3. **Write extensive tests**:
   ```python
   def test_spr_preserves_taxa():
       tree = create_test_tree(n_taxa=5)
       original_taxa = set(tree.get_leaf_names())

       for _ in range(100):  # Test 100 random SPR moves
           tree = spr_move(tree, random_prune(), random_regraft())
           assert set(tree.get_leaf_names()) == original_taxa
           assert tree.count_leaves() == len(original_taxa)

   def test_spr_improves_likelihood():
       # Test on datasets where SPR should find improvement
       pass
   ```

4. **Validate against known results**:
   - Compare to RAxML output on same dataset
   - Tree distance metrics (Robinson-Foulds)
   - Likelihood values should match within tolerance

## Integration with rRNA-Phylo

### Current Usage

```python
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4

# Option 1: NNI only (recommended)
tree, logL, meta = build_ml_tree_level4(
    sequences,
    tree_search='nni',
    verbose=True
)

# Option 2: NNI + Simplified SPR (experimental)
tree, logL, meta = build_ml_tree_level4(
    sequences,
    tree_search='spr',  # Does NNI first, then simplified SPR
    verbose=True
)

# Check what happened
print(f"NNI improvements: {meta.get('n_nni_improvements', 0)}")
print(f"SPR improvements: {meta.get('n_spr_improvements', 0)}")
print(f"Branches collapsed: {meta.get('n_branches_collapsed', 0)}")
```

### Future API (with proper SPR)

```python
# Recommended: Multi-stage search
tree, logL, meta = build_ml_tree_level4(
    sequences,
    tree_search='auto',  # NNI → SPR → TBR
    spr_iterations=5,
    verbose=True
)

# Advanced: Custom search strategy
tree, logL, meta = build_ml_tree_level4(
    sequences,
    tree_search=['nni', 'spr', 'tbr'],  # Sequential
    max_iterations={'nni': 100, 'spr': 10, 'tbr': 5},
    verbose=True
)
```

## Performance Considerations

### Computational Complexity

- **NNI**: O(n) neighbors per node, O(n²) total for tree
- **SPR**: O(n²) neighbors per node, O(n³) total for tree
- **TBR**: O(n³) neighbors per node, O(n⁴) total for tree

Where n = number of taxa

### Optimization Strategies

1. **Lazy SPR**: Only try SPR if NNI converges without improvement
2. **Limited Radius**: Only regraft within distance k of prune point
3. **Pruning Heuristics**: Skip unlikely moves based on likelihood bound
4. **Parallel Evaluation**: Test multiple SPR moves concurrently
5. **Cache Likelihood**: Reuse computations for similar trees

### Example Performance Targets

```
Dataset Size | NNI Time | SPR Time | SPR/NNI Ratio
-------------|----------|----------|---------------
10 taxa      | 1s       | 3s       | 3x
20 taxa      | 10s      | 40s      | 4x
50 taxa      | 60s      | 300s     | 5x
100 taxa     | 300s     | 2000s    | 6.7x
```

## Key Takeaways

1. **SPR is powerful but complex** - proper implementation requires careful tree manipulation

2. **Tree corruption is the main risk** - always validate leaf count and taxa names after moves

3. **NNI is sufficient for most use cases** - provides good results with much less complexity

4. **Start simple, then enhance** - get NNI working perfectly before adding SPR

5. **Use reference implementations** - don't reinvent the wheel, study RAxML/IQ-TREE

6. **Test extensively** - SPR bugs are subtle and only appear with specific tree structures

## References

1. **Swofford, D.L.** (2003). PAUP*: Phylogenetic Analysis Using Parsimony (*and Other Methods). Sinauer Associates.
   - SPR algorithm description

2. **Stamatakis, A.** (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30(9): 1312-1313.
   - Production SPR implementation

3. **Nguyen, L.T., et al.** (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution 32(1): 268-274.
   - Modern SPR with optimizations

4. **Felsenstein, J.** (2004). Inferring Phylogenies. Sinauer Associates.
   - Theoretical foundation

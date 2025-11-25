---
skill_name: consensus-tree-methods
description: Consensus tree methods for combining multiple phylogenetic trees with confidence scoring - includes Robinson-Foulds distance, majority-rule consensus, strict consensus, and bootstrap support calculation
triggers:
  - type: intent_pattern
    patterns:
      - consensus tree
      - combine trees
      - tree comparison
      - robinson-foulds
      - majority rule
      - strict consensus
      - bootstrap support
      - tree distance
      - merge trees
---

# Consensus Tree Methods for rRNA-Phylo

This skill provides comprehensive guidance for implementing consensus tree methods to combine multiple phylogenetic trees (UPGMA, BioNJ, ML) into a single consensus tree with confidence scores.

## Overview

**Goal**: Combine multiple phylogenetic trees to create a consensus tree that represents the most reliable relationships, with confidence scores for each branch.

**Why Consensus Trees?**
- Multiple tree-building methods can produce different topologies
- Consensus provides "forensics-grade" reliability by combining evidence
- Support values show which parts of the tree are well-supported
- Critical for publication-quality phylogenetic analysis

## Core Components

### 1. Tree Comparison Metrics

#### Robinson-Foulds (RF) Distance
**Purpose**: Measure topological difference between two trees

**Algorithm**:
```python
def robinson_foulds_distance(tree1: TreeNode, tree2: TreeNode) -> int:
    """
    Calculate RF distance = (splits in tree1 not in tree2) + (splits in tree2 not in tree1)

    Steps:
    1. Extract all bipartitions (splits) from tree1
    2. Extract all bipartitions (splits) from tree2
    3. Count symmetric difference

    Returns:
        int: RF distance (0 = identical topologies, max = 2*(n-3) for n taxa)
    """
```

**Key Concepts**:
- **Bipartition/Split**: A division of taxa into two groups by removing an edge
- **Normalized RF**: RF distance / max possible RF distance (0-1 scale)
- Works on unrooted trees (convert rooted trees first)

**Implementation Tips**:
- Use frozenset for bipartitions (hashable, unordered)
- Sort taxa names for consistent comparison
- Handle edge cases: identical trees, single taxon, two taxa

#### Branch Score Distance
**Purpose**: Weighted distance accounting for branch lengths

```python
def branch_score_distance(tree1: TreeNode, tree2: TreeNode) -> float:
    """
    Similar to RF but weights by branch length differences

    For each shared bipartition:
        distance += |branch_length1 - branch_length2|

    For non-shared bipartitions:
        distance += branch_length (treat missing as length 0)
    """
```

### 2. Bipartition Extraction

**Critical Function**:
```python
def get_bipartitions(tree: TreeNode) -> Set[frozenset]:
    """
    Extract all bipartitions from a tree.

    For each internal edge:
        1. Get all descendant taxa on one side
        2. Get all other taxa on other side
        3. Create frozenset of smaller group (canonical form)

    Returns:
        Set of bipartitions, each as frozenset of taxa
    """
```

**Example**:
```
Tree: ((A,B),(C,D))

Bipartitions:
- {A,B} | {C,D}  -> store as frozenset({A,B})
- {A} | {B,C,D}  -> store as frozenset({A})
- {B} | {A,C,D}  -> store as frozenset({B})
- {C} | {A,B,D}  -> store as frozenset({C})
- {D} | {A,B,C}  -> store as frozenset({D})
```

### 3. Consensus Algorithms

#### Majority-Rule Consensus (RECOMMENDED)
**Purpose**: Include splits that appear in >50% of input trees

**Algorithm**:
```python
def majority_rule_consensus(trees: List[TreeNode]) -> Tuple[TreeNode, Dict]:
    """
    Build consensus tree from multiple trees.

    Steps:
    1. Get all leaf names (taxa)
    2. Extract bipartitions from each tree
    3. Count frequency of each bipartition
    4. Keep bipartitions appearing in >50% of trees
    5. Build tree from compatible bipartitions (descending frequency)
    6. Assign support values (frequency %)

    Returns:
        (consensus_tree, support_values)
    """
```

**Support Values**:
- For each bipartition: support = (count / num_trees) * 100
- Example: If split appears in 2 out of 3 trees, support = 66.7%

**Compatibility Check**:
- Two bipartitions are compatible if they don't conflict
- Conflict: Neither is a subset of the other AND they overlap
- Build tree by adding bipartitions in order of support (greedy approach)

#### Strict Consensus
**Purpose**: Only include splits in ALL trees (100% agreement)

```python
def strict_consensus(trees: List[TreeNode]) -> TreeNode:
    """
    Like majority-rule but threshold = 100%

    Very conservative - results in polytomies (multifurcations)
    Use when you need absolute certainty
    """
```

#### Extended Majority-Rule
**Purpose**: Add compatible minority splits after majority splits

```python
def extended_majority_rule(trees: List[TreeNode]) -> TreeNode:
    """
    1. Add all majority-rule splits (>50%)
    2. Add compatible splits from 50% down to some threshold (e.g., 10%)
    3. Results in more resolved tree than strict consensus
    """
```

### 4. Building Tree from Bipartitions

**Algorithm**:
```python
def build_tree_from_bipartitions(
    taxa: Set[str],
    bipartitions: List[Tuple[frozenset, float]]  # (split, support)
) -> TreeNode:
    """
    Construct tree from set of compatible bipartitions.

    Steps:
    1. Sort bipartitions by size (smallest first) or support (highest first)
    2. Start with star tree (all taxa connected to root)
    3. For each bipartition:
        a. Find node where this split can be added
        b. Create new internal node
        c. Move appropriate children to new node
        d. Attach support value to new node

    Returns:
        TreeNode with support values
    """
```

### 5. Support Value Integration

**Add to TreeNode class**:
```python
@dataclass
class TreeNode:
    name: Optional[str] = None
    children: List['TreeNode'] = field(default_factory=list)
    dist: float = 0.0
    support: Optional[float] = None  # NEW: Bootstrap/consensus support
```

**Visualization**:
```
`-- Internal (dist: 0.05, support: 95.0)
    |-- seq1 (dist: 0.01)
    `-- seq2 (dist: 0.02)
```

## Implementation Structure

### File Organization
```
rrna_phylo/consensus/
├── __init__.py              # Export main functions
├── tree_distance.py         # RF distance, branch score
├── bipartitions.py          # Bipartition extraction and utils
├── consensus.py             # Majority-rule, strict consensus
└── tree_builder.py          # Build tree from bipartitions
```

### Key Functions to Export
```python
# rrna_phylo/consensus/__init__.py
from rrna_phylo.consensus.tree_distance import (
    robinson_foulds_distance,
    normalized_rf_distance,
    branch_score_distance
)
from rrna_phylo.consensus.consensus import (
    majority_rule_consensus,
    strict_consensus,
    extended_majority_rule
)

__all__ = [
    "robinson_foulds_distance",
    "normalized_rf_distance",
    "branch_score_distance",
    "majority_rule_consensus",
    "strict_consensus",
    "extended_majority_rule"
]
```

## Integration with build_trees()

### Enhanced Return Value
```python
def build_trees(sequences: List[Sequence], **kwargs) -> dict:
    """
    Returns:
        {
            'upgma': TreeNode,
            'bionj': TreeNode,
            'ml': (TreeNode, log_likelihood),
            'consensus': TreeNode,  # NEW
            'support_values': Dict[str, float],  # NEW
            'tree_distances': {  # NEW
                'upgma_vs_bionj': float,
                'upgma_vs_ml': float,
                'bionj_vs_ml': float
            },
            'type': SequenceType,
            'model': str
        }
    """
```

## Testing Strategy

### Test Cases

**1. Simple Identical Trees**
```python
def test_identical_trees():
    # Two identical trees should have RF distance = 0
    tree1 = parse_newick("((A,B),(C,D));")
    tree2 = parse_newick("((A,B),(C,D));")
    assert robinson_foulds_distance(tree1, tree2) == 0
```

**2. Different Topologies**
```python
def test_different_topologies():
    tree1 = parse_newick("((A,B),(C,D));")
    tree2 = parse_newick("((A,C),(B,D));")
    # Should have non-zero RF distance
    assert robinson_foulds_distance(tree1, tree2) > 0
```

**3. Majority-Rule Consensus**
```python
def test_majority_consensus():
    # 3 trees: 2 agree on (A,B), 1 disagrees
    trees = [
        parse_newick("((A,B),(C,D));"),
        parse_newick("((A,B),(C,D));"),
        parse_newick("((A,C),(B,D));")
    ]
    consensus, support = majority_rule_consensus(trees)
    # (A,B) split should have 66.7% support
    assert get_split_support(consensus, frozenset(['A','B'])) == pytest.approx(66.7)
```

**4. Real Tree Comparison**
```python
def test_upgma_vs_bionj():
    sequences = load_test_sequences()
    trees = build_trees(sequences)

    # Calculate distance between methods
    rf_dist = robinson_foulds_distance(trees['upgma'], trees['bionj'])

    # Trees should be similar but not necessarily identical
    max_rf = calculate_max_rf(len(sequences))
    normalized_dist = rf_dist / max_rf
    assert normalized_dist < 0.5  # Less than 50% different
```

## Common Pitfalls

### 1. **Rooting Issues**
- RF distance requires consistent rooting
- Convert to unrooted or root at same location
- Midpoint root or outgroup root

### 2. **Bipartition Representation**
- Always use canonical form (smaller set)
- Use frozenset for hashability
- Sort taxa names for consistency

### 3. **Trivial Splits**
- Single-taxon splits are trivial (always present)
- May want to exclude from RF calculation
- Focus on informative (internal) splits

### 4. **Support Value Interpretation**
- >95% = strong support
- 70-95% = moderate support
- <70% = weak support (polytomy recommended)

### 5. **Compatibility Checking**
- Must ensure added splits don't create conflicts
- Greedy addition (highest support first) usually works
- May need backtracking for complex cases

## Performance Considerations

- **Bipartition extraction**: O(n²) for n taxa
- **RF distance**: O(n) after bipartition extraction
- **Consensus building**: O(n² × m) for m trees, n taxa
- **Large trees**: Consider pruning or sampling

## Example Usage

```python
from rrna_phylo import build_trees, Sequence
from rrna_phylo.consensus import (
    majority_rule_consensus,
    robinson_foulds_distance
)

# Build trees with all methods
sequences = [...]
trees = build_trees(sequences)

# Get individual trees
upgma = trees['upgma']
bionj = trees['bionj']
ml_tree, logL = trees['ml']

# Compare trees
rf_upgma_bionj = robinson_foulds_distance(upgma, bionj)
print(f"RF distance UPGMA vs BioNJ: {rf_upgma_bionj}")

# Build consensus
all_trees = [upgma, bionj, ml_tree]
consensus, support = majority_rule_consensus(all_trees)

# Visualize with support values
print_tree_ascii(consensus, show_support=True)
```

## References

**Key Papers**:
- Robinson & Foulds (1981) - Original RF distance paper
- Felsenstein (1985) - Bootstrap confidence intervals
- Swofford (1991) - PAUP and consensus methods
- Bryant (2003) - Compatibility and consensus algorithms

**Recommended Reading**:
- Felsenstein "Inferring Phylogenies" (Chapter 28: Consensus Trees)
- Semple & Steel "Phylogenetics" (Chapter 5: Tree Metrics)

## Next Steps After Implementation

1. **Bootstrap Analysis** (Future enhancement)
   - Resample alignment columns
   - Rebuild trees from resampled data
   - Calculate bootstrap support from 100-1000 replicates

2. **Bayesian Support** (Advanced)
   - Posterior probabilities from Bayesian inference
   - More accurate than bootstrap for some cases
   - Requires MCMC sampling

3. **Quartet Distance** (Alternative metric)
   - More fine-grained than RF distance
   - Better for large trees
   - More computationally intensive

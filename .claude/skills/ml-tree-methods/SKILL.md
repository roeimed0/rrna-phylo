---
description: Maximum Likelihood phylogenetic tree inference - from simple GTR implementation to RAxML-NG on Windows
triggers:
  keywords:
    - maximum likelihood
    - ML tree
    - GTR model
    - RAxML
    - likelihood calculation
    - Felsenstein
  intent_patterns:
    - "build ML tree"
    - "maximum likelihood phylogenetics"
    - "GTR substitution model"
  file_patterns:
    - "**/ml_tree.py"
    - "**/raxml*.py"
enforcement: suggest
priority: high
---

# Maximum Likelihood Tree Methods

Guide to Maximum Likelihood (ML) phylogenetic inference from basic implementation to RAxML-NG on Windows.

## Overview

**Maximum Likelihood finds the tree that makes your observed data most probable.**

### Core Idea
```
For each possible tree:
  Calculate: P(alignment | tree, model)

Return tree with highest probability
```

### Why ML is Better than Distance Methods
- **Statistical framework**: Uses full alignment, not just pairwise distances
- **Model-based**: Accounts for rate variation, base frequencies
- **Confidence values**: Bootstrap support, likelihood scores
- **More accurate**: Especially for divergent sequences

---

## Part 1: GTR Substitution Model

### What is GTR?

**General Time Reversible** - the most general reversible DNA substitution model.

### The Math (High Level)

**DNA changes over time:**
```
A ↔ C, A ↔ G, A ↔ T, C ↔ G, C ↔ T, G ↔ T
```

**GTR has 9 parameters:**
- 6 rate parameters: r_AC, r_AG, r_AT, r_CG, r_CT, r_GT
- 4 base frequencies: π_A, π_C, π_G, π_T (only 3 free, since they sum to 1)

**Rate Matrix Q (4×4):**
```
        A           C           G           T
A    -μ_A      r_AC*π_C    r_AG*π_G    r_AT*π_T
C   r_AC*π_A     -μ_C      r_CG*π_G    r_CT*π_T
G   r_AG*π_A   r_CG*π_C      -μ_G      r_GT*π_T
T   r_AT*π_A   r_CT*π_C    r_GT*π_G      -μ_T
```

Where: μ_X = -(sum of other entries in row X) to make rows sum to 0

**Probability Matrix P(t):**
```python
P(t) = e^(Qt)  # Matrix exponential
```

P[i,j] = probability that nucleotide i becomes j after time t

### Implementation Pattern

```python
import numpy as np
from scipy.linalg import expm

class GTRModel:
    def __init__(self):
        # Nucleotide indices: A=0, C=1, G=2, T=3
        self.nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        # Base frequencies (estimate from data)
        self.base_freq = np.array([0.25, 0.25, 0.25, 0.25])

        # Rate parameters (estimate or use defaults)
        # Order: AC, AG, AT, CG, CT, GT
        self.rates = np.array([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
        # Note: Transitions (AG, CT) usually faster than transversions

    def estimate_base_frequencies(self, sequences):
        """Estimate π_A, π_C, π_G, π_T from alignment."""
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total = 0

        for seq in sequences:
            for base in seq.sequence.upper():
                if base in counts:
                    counts[base] += 1
                    total += 1

        self.base_freq = np.array([
            counts['A'] / total,
            counts['C'] / total,
            counts['G'] / total,
            counts['T'] / total
        ])

    def build_rate_matrix(self):
        """Build Q matrix from rates and base frequencies."""
        Q = np.zeros((4, 4))

        # Map rate indices
        rate_map = {
            (0,1): 0,  # A-C
            (0,2): 1,  # A-G
            (0,3): 2,  # A-T
            (1,2): 3,  # C-G
            (1,3): 4,  # C-T
            (2,3): 5   # G-T
        }

        # Fill off-diagonal
        for i in range(4):
            for j in range(4):
                if i != j:
                    key = (min(i,j), max(i,j))
                    rate_idx = rate_map[key]
                    Q[i,j] = self.rates[rate_idx] * self.base_freq[j]

        # Fill diagonal (make rows sum to 0)
        for i in range(4):
            Q[i,i] = -np.sum(Q[i,:])

        # Normalize (average rate = 1)
        mu = -np.sum(np.diagonal(Q) * self.base_freq)
        Q = Q / mu

        return Q

    def probability_matrix(self, branch_length, Q):
        """Calculate P(t) = e^(Qt) using matrix exponential."""
        return expm(Q * branch_length)
```

### Key Insights

1. **Reversibility**: Q[i,j] * π_i = Q[j,i] * π_j (detailed balance)
2. **Rows sum to 0**: Probability is conserved
3. **Transitions faster**: AG and CT usually have higher rates
4. **Time dependency**: P(t) changes with branch length

---

## Part 2: Likelihood Calculation

### Felsenstein's Pruning Algorithm (1981)

**The Problem**: Calculate P(alignment | tree)

**Naive approach**: O(4^n) where n = number of sequences (impossible!)

**Felsenstein's solution**: O(n) using dynamic programming

### Algorithm (High Level)

For each site (column) in alignment:

```
1. Start at leaves (we know their nucleotides)

2. Work up the tree:
   For each internal node:
     L[node][state] = P(data below | node has state)

   Calculate from children:
     L[node][A] = Σ_i P(A→i, t_left) * L[left_child][i] *
                  Σ_j P(A→j, t_right) * L[right_child][j]

3. At root:
   L[site] = Σ_state π[state] * L[root][state]

4. Total likelihood:
   L(tree) = Π_sites L[site]

   In practice, use log-likelihood:
   log L(tree) = Σ_sites log L[site]
```

### Implementation Pattern

```python
def calculate_site_likelihood(tree, site_data, model):
    """Calculate likelihood for one alignment column."""

    # Recursive function
    def conditional_likelihood(node, state):
        """L[node][state] = P(data below | node = state)"""

        if node.is_leaf():
            # Known state at leaf
            observed = site_data[node.id]
            return 1.0 if state == observed else 0.0

        # Internal node - recurse to children
        left_prob = 0.0
        right_prob = 0.0

        P_left = model.probability_matrix(node.left.distance)
        P_right = model.probability_matrix(node.right.distance)

        for child_state in range(4):  # A, C, G, T
            left_prob += P_left[state][child_state] * \
                        conditional_likelihood(node.left, child_state)
            right_prob += P_right[state][child_state] * \
                         conditional_likelihood(node.right, child_state)

        return left_prob * right_prob

    # Sum over all possible root states
    likelihood = 0.0
    for root_state in range(4):
        likelihood += model.base_freq[root_state] * \
                     conditional_likelihood(tree, root_state)

    return likelihood

def calculate_tree_likelihood(tree, alignment, model):
    """Calculate total likelihood across all sites."""
    log_likelihood = 0.0

    for site in range(alignment.shape[1]):
        site_data = alignment[:, site]
        L = calculate_site_likelihood(tree, site_data, model)
        log_likelihood += np.log(L)

    return log_likelihood
```

### Optimization: Log Space

Use log-likelihoods to avoid underflow:
```python
# Instead of: L = L1 * L2 * L3 ... (underflows!)
# Use: log(L) = log(L1) + log(L2) + log(L3) ...
```

---

## Part 3: Tree Search

### The Problem

**Tree space is huge!**
- For n sequences: (2n-3)!! possible unrooted trees
- n=10: 2 million trees
- n=20: 2×10^20 trees (can't try them all!)

### Solution: Heuristic Search

**Start with good initial tree** (from BioNJ or UPGMA)

**Improve iteratively** using local rearrangements

### Nearest Neighbor Interchange (NNI)

```
Original tree:      NNI move:

    ----A               ----C
   |                   |
---+                ---+
   |                   |
   |----B              |----B
   |                   |
   +----C              +----A
   |                   |
   +----D              +----D
```

**Algorithm:**
```
1. Start with initial tree T
2. Calculate likelihood L(T)
3. For each internal branch:
     Try both NNI rearrangements
     If likelihood improves, accept move
4. Repeat until no improvement
5. Return best tree
```

### Implementation Pattern

```python
def nni_search(tree, alignment, model):
    """Find ML tree using NNI search."""

    current_tree = tree
    current_likelihood = calculate_tree_likelihood(tree, alignment, model)
    improved = True

    while improved:
        improved = False

        # Try all NNI moves
        for branch in get_internal_branches(current_tree):
            for nni_tree in get_nni_neighbors(current_tree, branch):
                nni_likelihood = calculate_tree_likelihood(nni_tree, alignment, model)

                if nni_likelihood > current_likelihood:
                    current_tree = nni_tree
                    current_likelihood = nni_likelihood
                    improved = True
                    break

            if improved:
                break

    return current_tree, current_likelihood
```

---

## Part 4: RAxML-NG on Windows

### Why RAxML?

**State-of-the-art ML implementation:**
- Highly optimized (C++)
- Multiple tree search strategies
- Bootstrap support
- Parallelized

### Windows Installation

**Option 1: Pre-compiled Binary** (Recommended)

```bash
# Download from GitHub releases
https://github.com/amkozlov/raxml-ng/releases

# Get: raxml-ng_v1.X.X_windows_x86_64.zip
# Extract raxml-ng.exe to project root
```

**Option 2: Conda** (if available)

```bash
conda install -c bioconda raxml-ng
```

### Python Wrapper Pattern

```python
import subprocess
import os

class RAxMLRunner:
    def __init__(self, raxml_executable="raxml-ng.exe"):
        # Check for local executable first
        project_root = os.path.dirname(__file__)
        local_raxml = os.path.join(project_root, "raxml-ng.exe")

        if os.path.exists(local_raxml):
            self.raxml_exe = local_raxml
        else:
            self.raxml_exe = raxml_executable

    def build_tree(self, alignment_file, output_prefix, model="GTR+G"):
        """Run RAxML-NG tree search."""

        cmd = [
            self.raxml_exe,
            "--search",              # Tree search mode
            "--msa", alignment_file, # Input alignment
            "--model", model,        # Substitution model
            "--prefix", output_prefix,
            "--threads", "auto",     # Use available cores
            "--redo"                 # Overwrite existing
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"RAxML failed: {result.stderr}")

        # Parse output tree
        tree_file = f"{output_prefix}.raxml.bestTree"
        return parse_newick(tree_file)
```

### Common Models

**For DNA:**
- `GTR+G` - GTR with gamma rate variation (recommended)
- `GTR+I+G` - GTR with invariant sites + gamma
- `HKY+G` - Hasegawa-Kishino-Yano (simpler than GTR)

**For Protein:**
- `LG+G` - Le-Gascuel model
- `WAG+G` - Whelan and Goldman

### Bootstrap Support

```python
def build_tree_with_bootstrap(alignment_file, n_bootstrap=100):
    """Build tree with bootstrap support values."""

    cmd = [
        "raxml-ng.exe",
        "--all",                    # Search + bootstrap
        "--msa", alignment_file,
        "--model", "GTR+G",
        "--bs-trees", str(n_bootstrap),  # Bootstrap replicates
        "--prefix", "output"
    ]

    subprocess.run(cmd)

    # Output: output.raxml.support (tree with support values)
```

---

## Part 5: Complete ML Pipeline

### Simple Implementation (Our Code)

```python
from ml_tree import build_ml_tree

# Build ML tree (uses GTR + BioNJ starting tree)
tree = build_ml_tree(aligned_sequences)
```

**Pros:**
- Pure Python, works on Windows
- Understand every step
- Good for learning

**Cons:**
- Slower than RAxML
- Less optimized
- No full NNI search yet

### RAxML-NG (Production)

```python
from raxml import RAxMLRunner

runner = RAxMLRunner("raxml-ng.exe")
tree = runner.build_tree_from_sequences(
    aligned_sequences,
    output_prefix="my_tree",
    model="GTR+G",
    bootstrap=100
)
```

**Pros:**
- State-of-the-art accuracy
- Fast (optimized C++)
- Bootstrap support
- Widely accepted

**Cons:**
- External dependency
- Black box (harder to understand)

### When to Use Which

**Use simple implementation:**
- Learning ML methods
- Understanding the math
- Small datasets (<20 sequences)
- Quick prototyping

**Use RAxML-NG:**
- Publication-quality trees
- Large datasets (>20 sequences)
- Need bootstrap support
- Production/forensics use

---

## Part 6: Validation & Diagnostics

### Check Likelihood Values

```python
# Higher (less negative) = better fit
tree1_logL = -1500.32
tree2_logL = -1498.15  # Better!
```

### Bootstrap Support

```
     95
    /  \
   /    \
  72    100
 / \    / \
A   B  C   D
```

- **>95%**: Very strong support
- **75-95%**: Strong support
- **50-75%**: Moderate support
- **<50%**: Weak support

### Compare Models

Use AIC (Akaike Information Criterion):
```
AIC = -2*logL + 2*k

where k = number of parameters
```

Lower AIC = better model

---

## Part 7: Common Patterns

### Pattern 1: Quick ML Tree

```python
from fasta_parser import parse_fasta
from ml_tree import build_ml_tree

# Parse alignment
sequences = parse_fasta("aligned.fasta")

# Build ML tree (simple implementation)
tree = build_ml_tree(sequences)

print(tree.to_newick())
```

### Pattern 2: Production ML with RAxML

```python
from raxml import RAxMLRunner

runner = RAxMLRunner("raxml-ng.exe")

tree = runner.build_tree(
    "aligned.fasta",
    output_prefix="results/my_analysis",
    model="GTR+G",
    bootstrap=100
)
```

### Pattern 3: Compare Multiple Models

```python
models = ["JC", "HKY+G", "GTR+G"]
results = {}

for model in models:
    tree = runner.build_tree("aligned.fasta", f"test_{model}", model)
    results[model] = get_likelihood(f"test_{model}.raxml.log")

# Choose best model by AIC
best_model = min(results, key=lambda m: calculate_aic(results[m]))
```

---

## Key Takeaways

1. **GTR is the standard**: Most general reversible DNA model
2. **Felsenstein's algorithm is clever**: Makes likelihood calculation feasible
3. **Tree search is heuristic**: Can't try all trees, use NNI/SPR
4. **RAxML is industry standard**: Use for publication
5. **Our implementation teaches fundamentals**: Understand before using black box
6. **Bootstrap for confidence**: Need 100+ replicates for reliable support

## References

- Felsenstein (1981) "Evolutionary trees from DNA sequences: a maximum likelihood approach"
- Kozlov et al. (2019) "RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference"
- Yang (2006) "Computational Molecular Evolution" (textbook)

## See Also

- **distance-tree-methods** skill - UPGMA, BioNJ
- **tree-consensus-patterns** skill - Combining multiple trees
- **phylogenetic-methods** skill - Overview of all methods

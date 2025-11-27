---
description: Guide for fixing ML model selection to properly test different substitution models (JC69, K80, F81, HKY85, GTR) with their specific rate matrices
triggers:
  - keywords: ["model selection", "substitution model", "likelihood", "jc69", "k80", "f81", "hky85", "gtr"]
---

# ML Model Selection Fix - Proper Substitution Model Implementation

## Problem Statement

**Current Bug**: All substitution models return identical log-likelihood values.

**Evidence**:
```
Model                   LogL   Params          BIC      Delta
----------------------------------------------------------------------
JC69              -319712.30      172    640818.04       0.00 *
K80               -319712.30      173    640826.14       8.10
F81               -319712.30      175    640842.34      24.30
HKY85             -319712.30      176    640850.45      32.41
GTR               -319712.30      181    640890.95      72.91
```

**Root Cause**: Placeholder implementation always uses GTR likelihood regardless of requested model.

## Current Broken Implementation

**Location**: `backend/rrna_phylo/models/model_selection.py:119-156`

```python
def _compute_likelihood_with_rate_matrix(...):
    """
    Compute likelihood using a specific rate matrix.
    ...
    """
    # Lines 144-148: Comment admits this is wrong!
    # "For now, use a simplified approach"
    # "This is a placeholder - proper implementation would..."

    from rrna_phylo.models.ml_tree_level3 import compute_log_likelihood
    return compute_log_likelihood(tree, sequences, alpha=alpha)  # ALWAYS GTR!
```

## Why This Matters

Model selection is scientifically critical:

1. **Tests evolutionary assumptions** - Different models have different assumptions
2. **Prevents overfitting** - Simpler models preferred if they fit as well
3. **Standard practice** - All phylogenetic papers report model selection results
4. **AIC/BIC only works if models differ** - Currently just counting parameters

## Substitution Models Overview

### JC69 (Jukes-Cantor 1969)
- **Simplest model** - All substitutions equally likely
- **Parameters**: 0 (besides branch lengths)
- **Rate matrix**: All off-diagonal elements equal

```
Q_JC = μ * [[-3   1   1   1  ]
            [ 1  -3   1   1  ]
            [ 1   1  -3   1  ]
            [ 1   1   1  -3 ]]
```

### K80 (Kimura 1980)
- **Transition/transversion** - Different rates for ti/tv
- **Parameters**: 1 (kappa = ti/tv ratio)
- **Rate matrix**: Transitions (A↔G, C↔T) differ from transversions

```
Q_K80 = μ * [[-a-2b    a      b      b   ]
             [  a    -a-2b    b      b   ]
             [  b      b    -a-2b    a   ]
             [  b      b      a    -a-2b ]]
where a = kappa, b = 1
```

### F81 (Felsenstein 1981)
- **Empirical frequencies** - Uses observed base frequencies
- **Parameters**: 3 (πA, πC, πG, πT with constraint Σπ = 1)
- **Rate matrix**: Weighted by equilibrium frequencies

```
Q_F81 = μ * [[-       πC     πG     πT   ]
             [ πA     -      πG     πT   ]
             [ πA     πC     -      πT   ]
             [ πA     πC     πG     -    ]]
where diagonal = -(sum of row)
```

### HKY85 (Hasegawa-Kishino-Yano 1985)
- **Combines K80 + F81** - Unequal frequencies + ti/tv ratio
- **Parameters**: 4 (kappa + 3 frequencies)
- **Rate matrix**: K80 weighted by frequencies

```
Q_HKY = μ * [[-              κ*πC           πG           πT     ]
             [ κ*πA          -              πG           πT     ]
             [ πA            πC             -            κ*πT   ]
             [ πA            πC             κ*πG         -      ]]
```

### GTR (General Time Reversible)
- **Most complex** - Six exchangeability rates
- **Parameters**: 9 (6 rates + 3 frequencies)
- **Rate matrix**: Most general reversible model

```
Q_GTR = [[-              a*πC           b*πG           c*πT   ]
         [ a*πA          -              d*πG           e*πT   ]
         [ b*πA          d*πC           -              f*πT   ]
         [ c*πA          e*πC           f*πG           -      ]]
```

## Proper Fix Implementation

### Step 1: Create Rate Matrix Builder

Create `backend/rrna_phylo/models/rate_matrices.py`:

```python
"""
Substitution model rate matrices for phylogenetic likelihood.

Implements proper Q matrices for JC69, K80, F81, HKY85, and GTR models.
"""

import numpy as np
from typing import Dict, Optional


def get_jc69_rate_matrix() -> np.ndarray:
    """
    Jukes-Cantor 1969 model - all substitutions equally likely.

    Returns:
        4x4 rate matrix Q
    """
    Q = np.array([
        [-3.0,  1.0,  1.0,  1.0],
        [ 1.0, -3.0,  1.0,  1.0],
        [ 1.0,  1.0, -3.0,  1.0],
        [ 1.0,  1.0,  1.0, -3.0]
    ])
    return Q


def get_k80_rate_matrix(kappa: float) -> np.ndarray:
    """
    Kimura 1980 model - different rates for transitions vs transversions.

    Args:
        kappa: Transition/transversion rate ratio (typically 2-4)

    Returns:
        4x4 rate matrix Q
    """
    # Transitions (A<->G, C<->T) get rate kappa
    # Transversions get rate 1.0
    Q = np.array([
        [-(kappa + 2.0),  kappa,        1.0,         1.0        ],
        [ kappa,         -(kappa + 2.0), 1.0,         1.0        ],
        [ 1.0,            1.0,          -(kappa + 2.0), kappa    ],
        [ 1.0,            1.0,           kappa,       -(kappa + 2.0)]
    ])
    return Q


def get_f81_rate_matrix(freqs: np.ndarray) -> np.ndarray:
    """
    Felsenstein 1981 model - equal rates but empirical base frequencies.

    Args:
        freqs: Base frequencies [πA, πC, πG, πT]

    Returns:
        4x4 rate matrix Q
    """
    Q = np.zeros((4, 4))

    # Off-diagonal: rate to base j is πj
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i, j] = freqs[j]

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def get_hky85_rate_matrix(kappa: float, freqs: np.ndarray) -> np.ndarray:
    """
    Hasegawa-Kishino-Yano 1985 model - K80 with empirical frequencies.

    Args:
        kappa: Transition/transversion ratio
        freqs: Base frequencies [πA, πC, πG, πT]

    Returns:
        4x4 rate matrix Q
    """
    Q = np.zeros((4, 4))

    # Define transitions (A=0, C=1, G=2, T=3)
    transitions = [(0, 2), (2, 0), (1, 3), (3, 1)]  # A<->G, C<->T

    for i in range(4):
        for j in range(4):
            if i != j:
                if (i, j) in transitions:
                    # Transition: rate = kappa * πj
                    Q[i, j] = kappa * freqs[j]
                else:
                    # Transversion: rate = πj
                    Q[i, j] = freqs[j]

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def get_gtr_rate_matrix(
    exchangeabilities: np.ndarray,
    freqs: np.ndarray
) -> np.ndarray:
    """
    General Time Reversible model - most complex reversible model.

    Args:
        exchangeabilities: Six rates [rAC, rAG, rAT, rCG, rCT, rGT]
        freqs: Base frequencies [πA, πC, πG, πT]

    Returns:
        4x4 rate matrix Q
    """
    Q = np.zeros((4, 4))

    # Map exchangeabilities to matrix positions
    # Order: AC, AG, AT, CG, CT, GT
    r = exchangeabilities

    # Fill upper triangle (off-diagonal)
    Q[0, 1] = r[0] * freqs[1]  # A->C
    Q[0, 2] = r[1] * freqs[2]  # A->G
    Q[0, 3] = r[2] * freqs[3]  # A->T
    Q[1, 2] = r[3] * freqs[2]  # C->G
    Q[1, 3] = r[4] * freqs[3]  # C->T
    Q[2, 3] = r[5] * freqs[3]  # G->T

    # Fill lower triangle (reversibility: Qij * πi = Qji * πj)
    Q[1, 0] = r[0] * freqs[0]  # C->A
    Q[2, 0] = r[1] * freqs[0]  # G->A
    Q[3, 0] = r[2] * freqs[0]  # T->A
    Q[2, 1] = r[3] * freqs[1]  # G->C
    Q[3, 1] = r[4] * freqs[1]  # T->C
    Q[3, 2] = r[5] * freqs[2]  # T->G

    # Diagonal: -(sum of row)
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])

    return Q


def normalize_rate_matrix(Q: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    """
    Normalize rate matrix so mean rate = 1.0 substitution/site/time.

    Args:
        Q: Unnormalized rate matrix
        freqs: Equilibrium frequencies

    Returns:
        Normalized rate matrix
    """
    # Mean rate = -Σ(πi * Qii)
    mean_rate = -np.sum(freqs * np.diag(Q))

    if mean_rate > 0:
        Q = Q / mean_rate

    return Q
```

### Step 2: Integrate into Likelihood Calculator

Modify `backend/rrna_phylo/models/model_selection.py`:

```python
from rrna_phylo.models.rate_matrices import (
    get_jc69_rate_matrix,
    get_k80_rate_matrix,
    get_f81_rate_matrix,
    get_hky85_rate_matrix,
    get_gtr_rate_matrix,
    normalize_rate_matrix
)

def compute_likelihood_for_model(
    tree: TreeNode,
    sequences: List[Sequence],
    model_name: str,
    alpha: Optional[float] = None,
    verbose: bool = False
) -> Tuple[float, Dict]:
    """
    Compute likelihood for a specific substitution model.
    """
    # Get empirical base frequencies
    freqs = compute_empirical_frequencies(sequences)

    # Get model-specific rate matrix
    if model_name == 'JC69':
        Q = get_jc69_rate_matrix()
        Q = normalize_rate_matrix(Q, freqs)
        params = None

    elif model_name == 'K80':
        kappa = 2.0  # Initial guess
        Q = get_k80_rate_matrix(kappa)
        Q = normalize_rate_matrix(Q, freqs)
        params = np.array([kappa])

    elif model_name == 'F81':
        Q = get_f81_rate_matrix(freqs)
        Q = normalize_rate_matrix(Q, freqs)
        params = freqs

    elif model_name == 'HKY85':
        kappa = 2.0  # Initial guess
        Q = get_hky85_rate_matrix(kappa, freqs)
        Q = normalize_rate_matrix(Q, freqs)
        params = np.array([kappa])

    elif model_name == 'GTR':
        # For GTR, use existing Level 3 implementation
        log_likelihood = compute_log_likelihood(tree, sequences, alpha=alpha)
        model_params = {'model': 'GTR', 'frequencies': freqs, 'alpha': alpha}
        return log_likelihood, model_params

    else:
        raise ValueError(f"Unknown model: {model_name}")

    # Compute likelihood with this rate matrix
    log_likelihood = compute_likelihood_with_Q_matrix(
        tree, sequences, Q, freqs, alpha
    )

    model_params = {
        'model': model_name,
        'params': params,
        'frequencies': freqs,
        'alpha': alpha,
        'Q': Q
    }

    if verbose:
        print(f"{model_name}: LogL = {log_likelihood:.2f}")

    return log_likelihood, model_params
```

### Step 3: Implement Likelihood with Custom Q Matrix

Add to `backend/rrna_phylo/models/model_selection.py`:

```python
def compute_likelihood_with_Q_matrix(
    tree: TreeNode,
    sequences: List[Sequence],
    Q: np.ndarray,
    freqs: np.ndarray,
    alpha: Optional[float] = None
) -> float:
    """
    Compute log-likelihood using a specific Q matrix.

    This is similar to Level 3 compute_log_likelihood but uses
    a custom Q matrix instead of always using GTR.

    Args:
        tree: Phylogenetic tree
        sequences: Aligned sequences
        Q: 4x4 rate matrix
        freqs: Base frequencies
        alpha: Gamma shape parameter

    Returns:
        Log-likelihood
    """
    # Use site pattern compression
    from rrna_phylo.models.ml_tree_level3 import SitePatternCompressor
    compressor = SitePatternCompressor(sequences)

    # Calculate likelihood for each pattern
    log_likelihood = 0.0

    for pattern_idx in range(compressor.n_patterns):
        pattern = compressor.patterns[pattern_idx]
        count = compressor.pattern_counts[pattern_idx]

        # Calculate likelihood for this pattern using Q matrix
        L = calculate_pattern_likelihood_with_Q(
            tree, pattern, sequences, Q, freqs, alpha
        )

        if L > 0:
            log_likelihood += count * np.log(L)
        else:
            log_likelihood += count * (-1000.0)  # Penalty for impossible

    return log_likelihood
```

## Testing Strategy

### Test 1: Models Return Different Likelihoods
```python
def test_models_differ():
    """Test that different models return different likelihoods."""
    models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']
    likelihoods = {}

    for model in models:
        logL, _ = compute_likelihood_for_model(tree, sequences, model)
        likelihoods[model] = logL

    # Check all different
    unique_logLs = set(likelihoods.values())
    assert len(unique_logLs) == 5, f"Only {len(unique_logLs)} unique likelihoods!"

    print("Model likelihoods:")
    for model, logL in likelihoods.items():
        print(f"  {model}: {logL:.2f}")
```

### Test 2: GTR Should Fit Best
```python
def test_gtr_best_likelihood():
    """Test that GTR has highest likelihood (most parameters)."""
    models = ['JC69', 'K80', 'F81', 'HKY85', 'GTR']
    likelihoods = {m: compute_likelihood_for_model(tree, sequences, m)[0]
                   for m in models}

    # GTR should have highest logL (most flexible)
    assert likelihoods['GTR'] >= max(likelihoods.values())
```

### Test 3: Model Selection Works
```python
def test_model_selection_chooses_reasonable():
    """Test that model selection picks a reasonable model."""
    best_model, alpha, score, results = select_best_model_with_gamma(
        tree, sequences, verbose=True
    )

    # Should pick one of the middle models (not JC69, not GTR)
    # BIC typically favors HKY85 or F81+G for real data
    assert best_model in ['HKY85', 'HKY85+G', 'F81+G', 'GTR'], \
           f"Unexpected best model: {best_model}"

    # Likelihoods should differ
    logLs = [r['logL'] for r in results.values()]
    assert len(set(logLs)) > 1, "All models still have same likelihood!"
```

## Expected Results After Fix

- ✅ Each model returns different log-likelihood
- ✅ GTR has highest likelihood (most parameters)
- ✅ JC69 has lowest likelihood (fewest parameters)
- ✅ BIC selects intermediate model (e.g., HKY85+G)
- ✅ Model selection is scientifically meaningful

## Performance Impact

- **Minimal**: All models use site pattern compression
- **Each model test**: Still ~30-60 seconds
- **Total model selection**: Still ~20 minutes (20 models)
- **But**: Now scientifically valid!

## Alternative: Use External Tools

If implementing all rate matrices is too complex:

**Use IQ-TREE for model selection**:
```bash
iqtree -s alignment.fasta -m MFP  # ModelFinder Plus
```

IQ-TREE tests 88+ models and is much faster than our implementation.

## References

- Jukes & Cantor (1969) - "Evolution of protein molecules"
- Kimura (1980) - "A simple method for estimating evolutionary rates"
- Felsenstein (1981) - "Evolutionary trees from DNA sequences"
- Hasegawa et al. (1985) - "Dating of human-ape splitting by molecular clock"
- Tavaré (1986) - "Some probabilistic and statistical problems in the analysis of DNA sequences"

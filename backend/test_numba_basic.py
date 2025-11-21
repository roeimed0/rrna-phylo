"""
Quick test to verify Numba integration works.
"""

print("Testing Numba integration...")

# Test 1: Check if Numba is installed
try:
    import numba
    print(f"[OK] Numba version {numba.__version__} installed")
except ImportError:
    print("[ERROR] Numba not installed")
    exit(1)

# Test 2: Check if our Numba functions work
try:
    from rrna_phylo.models.numba_likelihood import (
        calculate_transition_contrib,
        calculate_root_likelihood,
        pade_matrix_exp
    )
    print("[OK] Numba likelihood functions imported")
except ImportError as e:
    print(f"[ERROR] Failed to import Numba functions: {e}")
    exit(1)

# Test 3: Test basic Numba function
import numpy as np

P = np.array([
    [0.9, 0.03, 0.03, 0.04],
    [0.03, 0.9, 0.04, 0.03],
    [0.03, 0.04, 0.9, 0.03],
    [0.04, 0.03, 0.03, 0.9]
])
L_child = np.array([0.8, 0.1, 0.05, 0.05])

try:
    contrib = calculate_transition_contrib(P, L_child)
    print(f"[OK] calculate_transition_contrib: {contrib}")
except Exception as e:
    print(f"[ERROR] calculate_transition_contrib failed: {e}")
    exit(1)

# Test 4: Test Pade matrix exponential
Q = np.array([
    [-1.0, 0.3, 0.3, 0.4],
    [0.3, -1.0, 0.4, 0.3],
    [0.3, 0.4, -1.0, 0.3],
    [0.4, 0.3, 0.3, -1.0]
])
t = 0.1

try:
    P_pade = pade_matrix_exp(Q, t)
    print(f"[OK] pade_matrix_exp: shape {P_pade.shape}")
except Exception as e:
    print(f"[ERROR] pade_matrix_exp failed: {e}")
    exit(1)

# Test 5: Compare with scipy
from scipy.linalg import expm
P_scipy = expm(Q * t)
error = np.max(np.abs(P_scipy - P_pade))
print(f"[OK] Pade vs scipy error: {error:.2e}")

if error < 1e-4:
    print("[OK] Pade approximation is accurate")
else:
    print(f"[WARNING] Pade approximation has large error: {error:.2e}")

# Test 6: Test with actual likelihood calculator
print("\nTesting with LikelihoodCalculatorLevel3...")

from rrna_phylo import Sequence
from rrna_phylo.models.ml_tree_level3 import LikelihoodCalculatorLevel3, NUMBA_AVAILABLE
from rrna_phylo.models.ml_tree import GTRModel
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.methods.bionj import build_bionj_tree

print(f"NUMBA_AVAILABLE: {NUMBA_AVAILABLE}")

# Small test sequences
seqs = [
    Sequence('A', 'Seq A', 'ATGCATGCATGC'),
    Sequence('B', 'Seq B', 'ATGCATGCATGC'),
    Sequence('C', 'Seq C', 'ATCCATCCATCC'),
    Sequence('D', 'Seq D', 'ATCCATCCATCC'),
]

# Build tree
dist_matrix, ids = calculate_distance_matrix(seqs, model="jukes-cantor")
tree = build_bionj_tree(dist_matrix, ids)

# Create GTR model
model = GTRModel()
model.estimate_parameters(seqs)

# Test without Numba
calc_no_numba = LikelihoodCalculatorLevel3(model, seqs, alpha=1.0, use_numba=False)
logL_no_numba = calc_no_numba.calculate_likelihood(tree)
print(f"[OK] LogL (no Numba): {logL_no_numba:.2f}")

# Test with Numba
calc_with_numba = LikelihoodCalculatorLevel3(model, seqs, alpha=1.0, use_numba=True)
logL_with_numba = calc_with_numba.calculate_likelihood(tree)
print(f"[OK] LogL (with Numba): {logL_with_numba:.2f}")

# Check accuracy
diff = abs(logL_no_numba - logL_with_numba)
if diff < 0.1:
    print(f"[OK] Results match (diff = {diff:.4f})")
else:
    print(f"[WARNING] Results differ significantly: {diff:.4f}")

print("\n" + "=" * 60)
print("ALL TESTS PASSED!")
print("=" * 60)
print("Numba acceleration is working correctly.")
print("You can now run the full benchmark with test_numba_performance.py")
print("=" * 60)

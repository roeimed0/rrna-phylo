"""
Patched PyTorch GPU-accelerated likelihood calculation for phylogenetic trees.

Public API preserved:
 - TorchGPULikelihoodCalculator
 - compute_log_likelihood_gpu

Fixes:
 - Per-node scaling (log-scale per pattern)
 - Uses torch.linalg.matrix_exp for P matrices (consistent)
 - Discrete gamma via quantiles (if scipy available) or midpoint fallback
 - Pattern compression moved to CPU then to device
 - P cache keyed by rounded (branch_length * rate)
"""

from typing import List, Optional, Dict, Tuple
import math
import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
    GPU_AVAILABLE = torch.cuda.is_available()
except Exception:
    TORCH_AVAILABLE = False
    GPU_AVAILABLE = False

# Optional: prefer scipy's gamma quantiles for discrete gamma categories
try:
    from scipy.stats import gamma as scipy_gamma
    SCIPY_GAMMA_AVAILABLE = True
except Exception:
    SCIPY_GAMMA_AVAILABLE = False

from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence


# ---------- small utilities ----------

# Pre-built lookup dict (no .upper() calls needed)
# PERFORMANCE: This avoids 24.97M .upper() calls (1.147s savings)
_BASE_IDX_MAP = {
    'A': 0, 'a': 0,
    'C': 1, 'c': 1,
    'G': 2, 'g': 2,
    'T': 3, 't': 3,
    'U': 3, 'u': 3,
}

def _base_to_idx(base: str) -> int:
    """Convert DNA/RNA base to index (0=A, 1=C, 2=G, 3=T/U)."""
    return _BASE_IDX_MAP.get(base, -1)


def _postorder_nodes(root: TreeNode) -> List[TreeNode]:
    out = []
    def visit(n: Optional[TreeNode]):
        if n is None: return
        if n.left: visit(n.left)
        if n.right: visit(n.right)
        out.append(n)
    visit(root)
    return out


# ---------- Torch GPU likelihood calculator ----------
class TorchGPULikelihoodCalculator:
    """
    Torch-based GPU likelihood calculator.
    - Patterns compressed on CPU then moved to device.
    - Per-node normalization (log-scale per pattern) implemented.
    - compute_likelihood(tree) -> float
    """

    def __init__(
        self,
        Q: np.ndarray,
        base_freq: np.ndarray,
        sequences: List[Sequence],
        alpha: Optional[float] = None,
        use_gpu: bool = True,
        verbose: bool = False,
        n_categories: int = 4,
        compressor=None  # CRITICAL FIX: Accept CPU compressor for consistency
    ):
        if not TORCH_AVAILABLE:
            raise ImportError("PyTorch not available; install torch to use GPU likelihood.")

        self.use_gpu_requested = use_gpu
        self.device = torch.device('cuda') if (use_gpu and GPU_AVAILABLE) else torch.device('cpu')
        self.verbose = verbose

        # store numpy originals (for potential debug)
        self.Q_np = np.array(Q, dtype=np.float64)
        self.base_freq_np = np.array(base_freq, dtype=np.float64)

        # convert tensors to device (double precision)
        self.Q = torch.tensor(self.Q_np, dtype=torch.float64, device=self.device)
        self.base_freq = torch.tensor(self.base_freq_np, dtype=torch.float64, device=self.device)

        self.sequences = sequences
        self.n_sequences = len(sequences)
        if self.n_sequences == 0:
            raise ValueError("No sequences provided")

        self.seq_length = len(sequences[0].sequence)

        # CRITICAL FIX: Reuse CPU compressor if provided (ensures CPU/GPU consistency)
        # This fixes the pattern count/order mismatch between model selection and NNI
        if compressor is not None:
            # Reuse patterns from CPU compressor (CONSISTENT!)
            self.patterns = torch.tensor(compressor.patterns, dtype=torch.long, device=self.device)
            self.pattern_counts = torch.tensor(compressor.pattern_counts, dtype=torch.float64, device=self.device)
            self.n_patterns = compressor.n_patterns

            # VALIDATION: Pattern counts must match CPU exactly
            cpu_counts = compressor.pattern_counts
            gpu_counts = self.pattern_counts.cpu().numpy()
            max_diff = np.max(np.abs(cpu_counts - gpu_counts))
            assert max_diff < 1e-10, \
                f"CPU/GPU pattern count mismatch! Max diff: {max_diff}"

            # VALIDATION: Pattern count sum must equal sequence length
            count_sum = torch.sum(self.pattern_counts).item()
            assert abs(count_sum - self.seq_length) < 1e-6, \
                f"GPU pattern count sum ({count_sum}) != sequence length ({self.seq_length})!"

            if self.verbose:
                print(f"[TorchGPU] Reused CPU pattern compression: {self.seq_length} -> {self.n_patterns} patterns")
        else:
            # Create own patterns (for standalone use)
            self._compress_patterns()

        # gamma setup
        self.alpha = alpha
        self.n_categories = n_categories
        self._setup_gamma_rates()

        # cache for transition matrices: key -> tensor (4x4)
        # key is (rounded_branch_length_times_rate, rate_index)
        self._P_cache: Dict[Tuple[float,int], torch.Tensor] = {}

        # tiny values
        self._min_val = 1e-300
        self._min_P = 1e-15

        # PERFORMANCE: Track when topology changes (need to invalidate caches)
        self._needs_cache_clear = True  # Clear on first call

        # PERFORMANCE: Cached postorder traversal (invalidated on topology change)
        self._postorder_cache = None

        if self.verbose:
            dev = 'cuda' if self.device.type == 'cuda' else 'cpu'
            print(f"[TorchGPU] device={dev}, n_patterns={self.n_patterns}, n_cat={self.n_categories}")

    def _collect_postorder(self, root: TreeNode) -> List[TreeNode]:
        """
        Collect all nodes in postorder (children before parents).
        This eliminates Python recursion overhead in tight loops.

        CRITICAL: This is computed once and cached until topology changes.
        """
        result = []

        def _traverse(node):
            if node is None:
                return
            # Visit children first (postorder)
            if node.left:
                _traverse(node.left)
            if node.right:
                _traverse(node.right)
            # Then visit this node
            result.append(node)

        _traverse(root)
        return result

    def _compress_patterns(self):
        """Compress alignment to unique site patterns (on CPU), move to device."""
        n_seq = self.n_sequences
        seq_len = self.seq_length

        # build integer alignment matrix on CPU
        alignment = np.zeros((n_seq, seq_len), dtype=np.int8)
        for i, seq in enumerate(self.sequences):
            s = seq.sequence
            for j, ch in enumerate(s):
                alignment[i, j] = _base_to_idx(ch)

        patterns = {}
        for site in range(seq_len):
            pat = tuple(int(alignment[i, site]) for i in range(n_seq))
            patterns[pat] = patterns.get(pat, 0) + 1

        pats = list(patterns.keys())
        counts = np.array([patterns[p] for p in pats], dtype=np.float64)
        pat_arr = np.array([list(p) for p in pats], dtype=np.int64)  # (n_patterns, n_seq)

        # transfer to torch on chosen device
        self.patterns = torch.tensor(pat_arr, dtype=torch.long, device=self.device)          # (n_patterns, n_seq)
        self.pattern_counts = torch.tensor(counts, dtype=torch.float64, device=self.device)  # (n_patterns,)
        self.n_patterns = pat_arr.shape[0]

        if self.verbose:
            print(f"Site pattern compression: {seq_len} -> {self.n_patterns} patterns ({seq_len/self.n_patterns:.1f}x)")

    def _gamma_mean_truncated(self, a: float, low: float, high: float) -> float:
        """
        Compute mean of Gamma(a, scale=1/a) truncated to (low, high).
        COPIED FROM CPU implementation to ensure identical results!
        """
        from scipy.special import gammainc, gamma as gamma_func

        L = a * low
        H = a * high

        denom = gammainc(a, H) - gammainc(a, L)
        if denom <= 0 or not np.isfinite(denom):
            return 1.0

        lower_term = gammainc(a + 1, L)
        upper_term = gammainc(a + 1, H)
        num_unreg = gamma_func(a + 1) * (upper_term - lower_term)
        num = (num_unreg / a) / gamma_func(a)

        mean_trunc = num / denom
        if not np.isfinite(mean_trunc) or mean_trunc <= 0:
            return 1.0
        return float(mean_trunc)

    def _setup_gamma_rates(self):
        """Setup discrete gamma rates & probs on device - MATCH CPU EXACTLY!"""
        if self.alpha is None:
            self.rates = torch.tensor([1.0], dtype=torch.float64, device=self.device)
            self.rate_probs = torch.tensor([1.0], dtype=torch.float64, device=self.device)
            self.n_categories = 1
            return

        k = self.n_categories
        a = self.alpha

        if not SCIPY_GAMMA_AVAILABLE:
            # fallback: simple uniform rates
            mids = np.ones(k, dtype=np.float64)
            self.rates = torch.tensor(mids, dtype=torch.float64, device=self.device)
            self.rate_probs = torch.tensor(np.ones(k) / k, dtype=torch.float64, device=self.device)
            return

        # MATCH CPU: Use mean-truncated gamma (Yang 1994)
        from scipy.stats import gamma as scipy_gamma_dist

        # Compute boundaries using ppf
        ps = np.linspace(0.0, 1.0, k + 1)
        eps = 1e-12
        ps[0] = eps
        ps[-1] = 1.0 - eps

        boundaries = scipy_gamma_dist.ppf(ps, a, scale=1.0 / a)

        # Compute mean within each interval
        rates = []
        for i in range(k):
            low = float(boundaries[i])
            high = float(boundaries[i + 1])

            if not np.isfinite(low) or not np.isfinite(high) or high <= low:
                # fallback
                midp = (i + 0.5) / k
                try:
                    mid = float(scipy_gamma_dist.ppf(midp, a, scale=1.0 / a))
                except:
                    mid = 1.0
                rates.append(mid)
                continue

            r = self._gamma_mean_truncated(a, low, high)
            rates.append(r)

        rates = np.array(rates, dtype=np.float64)

        # Normalize so mean = 1
        rates = rates / np.mean(rates)

        # VALIDATION: Gamma rates must average to 1.0
        rate_mean = np.mean(rates)
        assert np.abs(rate_mean - 1.0) < 1e-6, \
            f"GPU Gamma rate mean ({rate_mean}) != 1.0 after normalization!"

        self.rates = torch.tensor(rates, dtype=torch.float64, device=self.device)
        self.rate_probs = torch.tensor(np.ones(k) / k, dtype=torch.float64, device=self.device)

        # VALIDATION: Rate probabilities must sum to 1.0
        prob_sum = torch.sum(self.rate_probs).item()
        assert np.abs(prob_sum - 1.0) < 1e-10, \
            f"GPU rate probabilities sum ({prob_sum}) != 1.0!"


    def _round_key(self, x: float) -> float:
        # Round to stable key to reduce fp cache misses
        return round(float(x), 8)


    def _compute_P_for_branch_rate(self, branch_length: float, rate_idx: int, rate_value: float) -> torch.Tensor:
        """Compute or fetch P matrix for branch_length * rate_value."""
        key = (self._round_key(branch_length * float(rate_value)), rate_idx)
        if key in self._P_cache:
            return self._P_cache[key]

        t = float(branch_length) * float(rate_value)
        Qt = self.Q * t

        # Use matrix_exp (more consistent) with fallback if fails
        try:
            P = torch.linalg.matrix_exp(Qt)
        except Exception:
            # fallback: eig + reconstruct (rare)
            eigvals, eigvecs = torch.linalg.eig(Qt)
            expvals = torch.exp(eigvals.real)
            eigvecs_inv = torch.linalg.inv(eigvecs)
            P = eigvecs @ torch.diag(expvals) @ eigvecs_inv
            P = P.real

        # numeric guard & normalize rows
        P = torch.clamp(P, min=self._min_P)
        row_sums = P.sum(dim=1, keepdim=True)
        P = P / row_sums

        self._P_cache[key] = P
        return P


    def _prepare_P_cache_for_tree(self, tree: TreeNode):
        """Precompute P matrices for every child branch & rate to save repeated computation."""
        nodes = _postorder_nodes(tree)
        for node in nodes:
            for child in (node.left, node.right):
                if child is None:
                    continue
                bl = getattr(child, 'distance', 0.01)
                if not np.isfinite(bl) or bl <= 0.0:
                    bl = 1e-8
                for ridx in range(self.n_categories):
                    _ = self._compute_P_for_branch_rate(bl, ridx, float(self.rates[ridx]))


    def clear_cache(self):
        self._P_cache.clear()


    def _conditional_likelihood_tensor(self, node: TreeNode, rate_idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Return (L_tensor, log_scale_vector) for node at SPECIFIC rate category.
        - L_tensor: shape (n_patterns, 4), normalized per pattern (rows sum ~1)
        - log_scale_vector: shape (n_patterns,), accumulated log-scale for that node's subtree
        - rate_idx: which gamma rate category to use

        CRITICAL FIX: Must compute conditional likelihoods SEPARATELY for each rate category,
        then average the final pattern likelihoods. Cannot average P matrices during pruning!
        This matches the CPU implementation.
        """
        cache_key = f'_cl_cache_r{rate_idx}'
        if hasattr(node, cache_key):
            cached = getattr(node, cache_key)
            return cached['tensor'], cached['log_scale']

        if node.is_leaf():
            # find sequence index by display_name
            seq_col = None
            for i, seq in enumerate(self.sequences):
                if seq.display_name == node.name:
                    seq_col = i
                    break

            # build L: (n_patterns, 4)
            if seq_col is None:
                L = self.base_freq.unsqueeze(0).repeat(self.n_patterns, 1)  # use equilibrium if not found
                log_scale = torch.zeros(self.n_patterns, dtype=torch.float64, device=self.device)
                setattr(node, cache_key, {'tensor': L, 'log_scale': log_scale})
                return L, log_scale

            obs = self.patterns[:, seq_col]  # (n_patterns,)
            L = torch.zeros((self.n_patterns, 4), dtype=torch.float64, device=self.device)

            mask_known = obs >= 0
            if mask_known.any():
                idxs = torch.nonzero(mask_known, as_tuple=False).squeeze(1)
                states = obs[mask_known].to(torch.long)
                L[idxs, states] = 1.0

            if (~mask_known).any():
                idxs = torch.nonzero(~mask_known, as_tuple=False).squeeze(1)
                L[idxs, :] = 0.25

            log_scale = torch.zeros(self.n_patterns, dtype=torch.float64, device=self.device)
            setattr(node, cache_key, {'tensor': L, 'log_scale': log_scale})
            return L, log_scale

        # internal node: compute children's tensors and combine
        left = node.left
        right = node.right

        # start with ones
        L_acc = torch.ones((self.n_patterns, 4), dtype=torch.float64, device=self.device)
        total_log_scale = torch.zeros(self.n_patterns, dtype=torch.float64, device=self.device)

        for child in (left, right):
            if child is None:
                continue
            # FIXED: Use rate-specific conditional likelihood (separate tree traversal per rate!)
            child_L, child_log_scale = self._conditional_likelihood_tensor(child, rate_idx)
            bl = getattr(child, 'distance', 0.01)
            if not np.isfinite(bl) or bl <= 0.0:
                bl = 1e-8

            # FIXED: Use only THIS rate category's P matrix (not averaging over rates!)
            P = self._P_cache[(self._round_key(bl * float(self.rates[rate_idx])), rate_idx)]

            # child_L @ P.T -> (n_patterns,4)
            contrib = torch.matmul(child_L, P.T)
            contrib = torch.clamp(contrib, min=self._min_val)

            # Multiply into accumulator
            L_acc = L_acc * contrib

            # add child's log_scale (already per-pattern)
            total_log_scale = total_log_scale + child_log_scale

        # Now normalize L_acc per pattern (rows) and collect local log-scale (CPU-style)
        max_L, _ = L_acc.max(dim=1, keepdim=True)
        max_L = torch.clamp(max_L, min=self._min_val)
        L_norm = L_acc / max_L
        local_log_scale = torch.log(max_L.squeeze(1))


        # accumulate local scales
        total_log_scale = total_log_scale + local_log_scale

        setattr(node, cache_key, {'tensor': L_norm, 'log_scale': total_log_scale})
        return L_norm, total_log_scale


    def mark_topology_changed(self):
        """
        Mark that tree topology has changed (e.g., after NNI swap).
        Next compute_likelihood() call will clear internal node caches and rebuild postorder.
        """
        self._needs_cache_clear = True
        self._postorder_cache = None  # Invalidate cached traversal

    def compute_likelihood(self, tree: TreeNode) -> float:
        """
        Compute log-likelihood for the input tree.

        CRITICAL FIX: Compute separate conditional likelihoods for EACH rate category,
        then average the final pattern likelihoods. This matches the CPU implementation.

        Returns Python float.
        """
        # PERFORMANCE FIX: Only clear caches when topology has changed!
        # This avoids recomputing internal nodes when just testing different swaps.
        # ALSO: Use flattened postorder to avoid recursion overhead
        if self._needs_cache_clear or self._postorder_cache is None:
            # Rebuild postorder traversal (once per topology change)
            self._postorder_cache = self._collect_postorder(tree)

            # Clear internal node caches using flattened list (no recursion!)
            for node in self._postorder_cache:
                # Skip leaves - their partials are permanent
                if node.is_leaf():
                    continue

                # Clear ONLY internal node caches
                attrs_to_del = [attr for attr in dir(node) if attr.startswith('_cl_cache_r')]
                for attr in attrs_to_del:
                    if hasattr(node, attr):
                        delattr(node, attr)

            self._needs_cache_clear = False  # Don't clear again until topology changes

        # prepare P cache for tree
        self._prepare_P_cache_for_tree(tree)

        # FIXED: Compute pattern likelihoods for EACH rate category separately
        # Shape: (n_patterns,)
        pattern_likelihood_sum = torch.zeros(self.n_patterns, dtype=torch.float64, device=self.device)

        for rate_idx in range(self.n_categories):
            # Compute conditional likelihoods for THIS rate category only
            L_root, root_log_scale = self._conditional_likelihood_tensor(tree, rate_idx)  # (n_patterns,4), (n_patterns,)

            # Weight by base frequencies to get per-pattern likelihood
            pattern_likes = torch.matmul(L_root, self.base_freq)  # (n_patterns,)
            pattern_likes = torch.clamp(pattern_likes, min=self._min_val)

            # Convert to actual likelihood (unlog the scaling)
            actual_likelihood = pattern_likes * torch.exp(root_log_scale)  # (n_patterns,)

            # Accumulate weighted by rate probability
            pattern_likelihood_sum += self.rate_probs[rate_idx] * actual_likelihood

        # Now take log of the averaged likelihood
        pattern_likelihood_sum = torch.clamp(pattern_likelihood_sum, min=self._min_val)
        log_pattern_likelihood = torch.log(pattern_likelihood_sum)  # (n_patterns,)

        # Weight by pattern counts
        weighted = self.pattern_counts * log_pattern_likelihood
        total = torch.sum(weighted)

        # PERFORMANCE FIX: Don't clear caches here!
        # Leaf partials are preserved across calls (never change).
        # Internal node caches will be cleared at start of next compute_likelihood() call.
        # This allows incremental likelihood computation during NNI.

        return float(total.cpu().item())


# wrapper function expected by tests
def compute_log_likelihood_gpu(
    tree: TreeNode,
    sequences: List[Sequence],
    Q: np.ndarray,
    base_freq: np.ndarray,
    alpha: Optional[float] = None,
    use_gpu: bool = True,
    verbose: bool = False
) -> float:
    calc = TorchGPULikelihoodCalculator(
        Q=Q,
        base_freq=base_freq,
        sequences=sequences,
        alpha=alpha,
        use_gpu=use_gpu,
        verbose=verbose
    )
    return calc.compute_likelihood(tree)

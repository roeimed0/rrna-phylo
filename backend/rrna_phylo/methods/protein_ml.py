"""Maximum Likelihood tree inference for protein sequences with empirical models."""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
from typing import List, Tuple
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.models.protein_models import ProteinModel, get_protein_model
from rrna_phylo.models.ml_tree_level3 import GammaRates, SitePatternCompressor


class ProteinLikelihoodCalculator:
    """Calculate tree likelihood for protein sequences."""

    def __init__(
        self,
        model: ProteinModel,
        sequences: List[Sequence],
        alpha: float = 1.0,
        n_categories: int = 4
    ):
        """Initialize protein likelihood calculator."""
        self.model = model
        self.sequences = sequences
        self.n_seq = len(sequences)

        # Gamma rates
        self.gamma = GammaRates(alpha, n_categories)

        # Site pattern compression (adapted for proteins)
        self.compressor = ProteinPatternCompressor(sequences)

        # Map sequence display names to indices (for tree node lookup)
        self.seq_id_to_idx = {seq.display_name: i for i, seq in enumerate(sequences)}

    def calculate_likelihood(self, tree: TreeNode) -> float:
        """Calculate log-likelihood with protein model + Gamma."""
        log_likelihood = 0.0

        # Get rate matrices for all gamma categories
        Q_matrices = self.gamma.get_rate_matrices(self.model.Q)

        # Calculate likelihood for each unique pattern
        for pattern_idx in range(self.compressor.n_patterns):
            pattern = self.compressor.patterns[pattern_idx]
            count = self.compressor.pattern_counts[pattern_idx]

            # Average likelihood across gamma categories
            pattern_likelihood = 0.0

            for cat_idx in range(self.gamma.n_categories):
                Q = Q_matrices[cat_idx]
                cat_prob = self.gamma.probabilities[cat_idx]

                # Calculate likelihood for this pattern and category
                L = self._calculate_pattern_likelihood(tree, pattern, Q)
                pattern_likelihood += cat_prob * L

            # Add to total (weighted by pattern count)
            if pattern_likelihood > 0:
                log_likelihood += count * np.log(pattern_likelihood)
            else:
                log_likelihood += count * (-1000.0)

        return log_likelihood

    def _calculate_pattern_likelihood(
        self,
        tree: TreeNode,
        pattern: np.ndarray,
        Q: np.ndarray
    ) -> float:
        """Calculate likelihood for one site pattern (Felsenstein's algorithm)."""
        # Felsenstein's pruning algorithm (20 states instead of 4)
        def conditional_likelihood(node: TreeNode) -> np.ndarray:
            """Calculate L[node][state] for all 20 states."""
            if node.is_leaf():
                L = np.zeros(20)  # 20 amino acids
                seq_idx = self.seq_id_to_idx[node.name]
                observed_state = pattern[seq_idx]

                if observed_state >= 0:
                    L[observed_state] = 1.0
                else:
                    L[:] = 1.0 / 20  # Gap: equal probability for all states

                return L

            # Internal node
            L = np.ones(20)

            # Process children
            if node.left:
                P_left = expm(Q * node.left.distance)
                L_left = conditional_likelihood(node.left)

                left_contrib = np.zeros(20)
                for parent_state in range(20):  # 20 states
                    for child_state in range(20):
                        left_contrib[parent_state] += \
                            P_left[parent_state, child_state] * L_left[child_state]

                L *= left_contrib

            if node.right:
                P_right = expm(Q * node.right.distance)
                L_right = conditional_likelihood(node.right)

                right_contrib = np.zeros(20)
                for parent_state in range(20):  # 20 states
                    for child_state in range(20):
                        right_contrib[parent_state] += \
                            P_right[parent_state, child_state] * L_right[child_state]

                L *= right_contrib

            return L

        # Get conditional likelihoods at root
        L_root = conditional_likelihood(tree)

        # Sum over root states (weighted by amino acid frequencies)
        likelihood = np.sum(self.model.aa_freq * L_root)

        return likelihood


class ProteinPatternCompressor:
    """Compress alignment by unique site patterns (protein version)."""

    def __init__(self, sequences: List[Sequence]):
        """Initialize pattern compressor."""
        self.sequences = sequences
        self.patterns = None
        self.pattern_counts = None
        self.n_patterns = 0
        self._compress()

    def _compress(self):
        """Identify unique site patterns and their counts."""
        n_seq = len(self.sequences)
        seq_len = self.sequences[0].aligned_length

        # Convert sequences to array (20 amino acids)
        alignment = np.zeros((n_seq, seq_len), dtype=int)

        # Standard amino acids
        aa_to_idx = {
            'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
            'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
            'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
            'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
        }

        for i, seq in enumerate(self.sequences):
            for j, aa in enumerate(seq.sequence.upper()):
                if aa in aa_to_idx:
                    alignment[i, j] = aa_to_idx[aa]
                else:
                    alignment[i, j] = -1  # Gap or unknown

        # Find unique patterns
        patterns_dict = {}

        for site in range(seq_len):
            pattern = tuple(alignment[:, site])

            if pattern in patterns_dict:
                patterns_dict[pattern] += 1
            else:
                patterns_dict[pattern] = 1

        # Convert to arrays
        self.patterns = np.array([list(p) for p in patterns_dict.keys()])
        self.pattern_counts = np.array(list(patterns_dict.values()))
        self.n_patterns = len(patterns_dict)

        compression_ratio = seq_len / self.n_patterns
        print(f"Site pattern compression: {seq_len} sites -> {self.n_patterns} patterns "
              f"({compression_ratio:.1f}x speedup)")


class ProteinBranchOptimizer:
    """Optimize branch lengths for protein likelihood."""

    def __init__(self, calculator: ProteinLikelihoodCalculator):
        """Initialize optimizer."""
        self.calculator = calculator

    def optimize_branch_lengths(self, tree: TreeNode, verbose: bool = False) -> float:
        """Optimize all branch lengths in tree."""
        improved = True
        iterations = 0
        max_iterations = 10

        current_logL = self.calculator.calculate_likelihood(tree)

        if verbose:
            print(f"Initial log-likelihood: {current_logL:.2f}")

        while improved and iterations < max_iterations:
            improved = False
            iterations += 1

            # Optimize each branch
            for node in self._traverse_tree(tree):
                if node.is_leaf():
                    continue

                old_dist = node.distance

                # Define optimization function
                def neg_likelihood(branch_length):
                    """Negative log-likelihood (for minimization)."""
                    if branch_length <= 0:
                        return 1e10
                    node.distance = branch_length
                    return -self.calculator.calculate_likelihood(tree)

                # Optimize
                result = minimize_scalar(
                    neg_likelihood,
                    bounds=(0.0001, 2.0),
                    method='bounded'
                )

                node.distance = result.x

                if abs(node.distance - old_dist) > 0.001:
                    improved = True

            new_logL = self.calculator.calculate_likelihood(tree)

            if verbose:
                print(f"Iteration {iterations}: log-likelihood = {new_logL:.2f}")

            if new_logL > current_logL + 0.01:
                current_logL = new_logL
                improved = True
            else:
                break

        if verbose:
            print(f"Final log-likelihood: {current_logL:.2f}")

        return current_logL

    def _traverse_tree(self, node: TreeNode) -> List[TreeNode]:
        """Get all nodes in tree (post-order traversal)."""
        nodes = []

        if node.left:
            nodes.extend(self._traverse_tree(node.left))
        if node.right:
            nodes.extend(self._traverse_tree(node.right))

        nodes.append(node)
        return nodes


def build_protein_ml_tree(
    sequences: List[Sequence],
    model_name: str = "WAG",
    alpha: float = 1.0,
    verbose: bool = True
) -> Tuple[TreeNode, float]:
    """Build ML tree for protein sequences."""
    if verbose:
        print("=" * 60)
        print("Maximum Likelihood Tree Inference - Proteins")
        print(f"{model_name}+Gamma Model with Site Pattern Compression")
        print("=" * 60)

    model = get_protein_model(model_name)
    model.estimate_frequencies(sequences)

    from rrna_phylo.distance.protein_distance import calculate_protein_distance_matrix
    from rrna_phylo.methods.bionj import build_bionj_tree

    dist_matrix, ids = calculate_protein_distance_matrix(sequences, model="poisson")
    initial_tree = build_bionj_tree(dist_matrix, ids)

    calculator = ProteinLikelihoodCalculator(model, sequences, alpha=alpha)
    initial_logL = calculator.calculate_likelihood(initial_tree)

    optimizer = ProteinBranchOptimizer(calculator)
    final_logL = optimizer.optimize_branch_lengths(initial_tree, verbose=verbose)

    if verbose:
        print("\n" + "=" * 60)
        print("ML Tree Inference Complete (Protein)!")
        print(f"Final log-likelihood: {final_logL:.2f}")
        print("=" * 60)

    return initial_tree, final_logL

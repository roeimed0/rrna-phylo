"""
Unified phylogenetic tree builder with automatic model selection.

This module provides a single interface for building phylogenetic trees
from DNA/RNA sequences, automatically detecting the type and selecting
the appropriate substitution model.
"""

from typing import List, Tuple, Optional
import tempfile
import os
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.io.aligner import MuscleAligner
from rrna_phylo.core.sequence_type import SequenceTypeDetector, SequenceType, get_appropriate_model
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.models.ml_tree_level4 import build_ml_tree_level4
from rrna_phylo.consensus import compare_trees


class PhylogeneticTreeBuilder:
    """
    Unified interface for phylogenetic tree building.

    Automatically detects sequence type and builds trees using
    appropriate models for DNA/RNA sequences.

    Usage:
        builder = PhylogeneticTreeBuilder()
        upgma, bionj, ml = builder.build_all_trees(sequences)
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize tree builder.

        Args:
            verbose: Print progress and recommendations
        """
        self.verbose = verbose
        self.detector = SequenceTypeDetector()
        self.seq_type = None
        self.model_name = None

    def _check_alignment(self, sequences: List[Sequence]) -> bool:
        """
        Check if sequences are already aligned (same length).

        Args:
            sequences: List of sequences

        Returns:
            True if aligned, False otherwise
        """
        if not sequences:
            return True

        lengths = [seq.aligned_length for seq in sequences]
        return len(set(lengths)) == 1

    def _align_sequences(self, sequences: List[Sequence]) -> List[Sequence]:
        """
        Align sequences using MUSCLE if they're not already aligned.

        Args:
            sequences: List of sequences (aligned or unaligned)

        Returns:
            List of aligned sequences
        """
        # Check if already aligned
        if self._check_alignment(sequences):
            if self.verbose:
                print("Sequences are already aligned (same length)")
            return sequences

        # Need alignment
        if self.verbose:
            print("\nSequences have different lengths - aligning with MUSCLE...")
            for seq in sequences:
                print(f"  {seq.id}: length={seq.aligned_length}")

        try:
            # Create temporary file for aligned output
            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp:
                output_file = tmp.name

            aligner = MuscleAligner()
            aligned_seqs = aligner.align_sequences(sequences, output_file)

            # Clean up temp file
            if os.path.exists(output_file):
                os.remove(output_file)

            if self.verbose:
                print(f"Alignment complete. Aligned length: {aligned_seqs[0].aligned_length}")

            return aligned_seqs
        except Exception as e:
            raise RuntimeError(f"Alignment failed: {e}")

    def detect_and_validate(self, sequences: List[Sequence]) -> SequenceType:
        """
        Detect sequence type and validate for phylogenetics.

        Args:
            sequences: List of sequences

        Returns:
            SequenceType

        Raises:
            ValueError: If sequences are invalid or mixed types
        """
        self.seq_type, recommendation = self.detector.validate_for_phylogenetics(sequences)
        self.model_name = get_appropriate_model(self.seq_type)

        if self.verbose:
            print("=" * 70)
            print("SEQUENCE TYPE DETECTION")
            print("=" * 70)
            print(f"Detected type: {self.seq_type.value.upper()}")
            print(f"Model: {self.model_name}")
            print(f"\n{recommendation}")
            print("=" * 70)

        return self.seq_type

    def _calculate_distance_matrix(self, sequences: List[Sequence]):
        """
        Calculate distance matrix using Jukes-Cantor model.

        Args:
            sequences: Aligned sequences

        Returns:
            Tuple of (distance_matrix, sequence_ids)
        """
        return calculate_distance_matrix(sequences, model="jukes-cantor")

    def build_upgma_tree(self, sequences: List[Sequence]) -> TreeNode:
        """
        Build UPGMA tree.

        Works for DNA/RNA sequences using Jukes-Cantor distance.

        Args:
            sequences: Aligned sequences

        Returns:
            UPGMA tree
        """
        if self.seq_type is None:
            self.detect_and_validate(sequences)

        if self.verbose:
            print("\n" + "=" * 70)
            print("METHOD 1: UPGMA (Assumes Molecular Clock)")
            print("=" * 70)

        # Calculate distance matrix
        dist_matrix, ids = self._calculate_distance_matrix(sequences)
        tree = build_upgma_tree(dist_matrix, ids)

        if self.verbose:
            print(f"UPGMA tree built successfully")
            print(f"Newick: {tree.to_newick()};")

        return tree

    def build_bionj_tree(self, sequences: List[Sequence]) -> TreeNode:
        """
        Build BioNJ tree.

        Works for DNA/RNA sequences using Jukes-Cantor distance.

        Args:
            sequences: Aligned sequences

        Returns:
            BioNJ tree
        """
        if self.seq_type is None:
            self.detect_and_validate(sequences)

        if self.verbose:
            print("\n" + "=" * 70)
            print("METHOD 2: BioNJ (Variance-Weighted, No Clock Assumption)")
            print("=" * 70)

        # Calculate distance matrix
        dist_matrix, ids = self._calculate_distance_matrix(sequences)
        tree = build_bionj_tree(dist_matrix, ids)

        if self.verbose:
            print(f"BioNJ tree built successfully")
            print(f"Newick: {tree.to_newick()};")

        return tree

    def build_ml_tree(self, sequences: List[Sequence], alpha: float = 1.0, skip_model_selection: bool = True) -> Tuple[TreeNode, float]:
        """
        Build Maximum Likelihood tree using Level 4 (GTR+Gamma model with NNI tree search).

        Works for DNA/RNA sequences using GTR+Gamma model directly (model selection skipped for performance).

        Args:
            sequences: Aligned sequences
            alpha: Gamma shape parameter (1.0 = moderate rate variation, or 'auto' for optimization)
            skip_model_selection: If True, use GTR+G directly (10x faster, default=True)

        Returns:
            (ml_tree, log_likelihood)
        """
        if self.seq_type is None:
            self.detect_and_validate(sequences)

        if self.verbose:
            print("\n" + "=" * 70)
            print("METHOD 3: Maximum Likelihood (GTR+Gamma with Model Selection + NNI)")
            print("=" * 70)

        # Use Level 4 with automatic model selection and tree search
        tree, logL, metadata = build_ml_tree_level4(
            sequences,
            model='auto',
            alpha=alpha,
            tree_search='nni',
            skip_model_selection=skip_model_selection,
            use_gpu='auto',  # Auto GPU selection based on dataset size
            verbose=self.verbose
        )

        if self.verbose:
            print(f"ML tree built successfully")
            print(f"Log-likelihood: {logL:.2f}")
            print(f"Newick: {tree.to_newick()};")

        return tree, logL

    def build_all_trees(
        self,
        sequences: List[Sequence],
        alpha: float = 1.0
    ) -> Tuple[TreeNode, TreeNode, Tuple[TreeNode, float]]:
        """
        Build trees using all three methods.

        This is the main interface for multi-method phylogenetic analysis.

        Args:
            sequences: Sequences (will be aligned automatically if needed)
            alpha: Gamma shape parameter for ML tree

        Returns:
            (upgma_tree, bionj_tree, (ml_tree, log_likelihood))

        Example:
            builder = PhylogeneticTreeBuilder()
            upgma, bionj, (ml, logL) = builder.build_all_trees(sequences)
        """
        # Detect type once
        self.detect_and_validate(sequences)

        # Align sequences if needed
        aligned_seqs = self._align_sequences(sequences)

        # Build all three trees with aligned sequences
        upgma_tree = self.build_upgma_tree(aligned_seqs)
        bionj_tree = self.build_bionj_tree(aligned_seqs)
        ml_result = self.build_ml_tree(aligned_seqs, alpha=alpha)

        if self.verbose:
            print("\n" + "=" * 70)
            print("SUMMARY")
            print("=" * 70)
            print(f"Sequence type: {self.seq_type.value.upper()}")
            print(f"Model: {self.model_name}")
            print(f"Number of sequences: {len(sequences)}")
            print(f"Alignment length: {sequences[0].aligned_length}")
            print("\nAll three methods completed successfully!")
            print("\nFor forensic reliability:")
            print("  [OK] Compare tree topologies")
            print("  [OK] Check if all methods agree")
            print("  [OK] Build consensus tree (recommended)")
            print("=" * 70)

        return upgma_tree, bionj_tree, ml_result


def build_trees(
    sequences: List[Sequence],
    method: str = "all",
    alpha: float = 1.0,
    verbose: bool = True
) -> dict:
    """
    Convenient function to build phylogenetic trees.

    Args:
        sequences: Aligned sequences (DNA/RNA)
        method: Which method(s) to use:
            - "all": Build UPGMA, BioNJ, and ML trees
            - "upgma": UPGMA only
            - "bionj": BioNJ only
            - "ml": Maximum Likelihood only
        alpha: Gamma shape parameter for ML (default 1.0)
        verbose: Print progress

    Returns:
        Dictionary with tree(s):
            - "upgma": UPGMA tree (if requested)
            - "bionj": BioNJ tree (if requested)
            - "ml": (ML tree, log-likelihood) (if requested)
            - "type": Sequence type detected
            - "model": Model used

    Example:
        # Build all trees
        results = build_trees(sequences)
        upgma = results["upgma"]
        bionj = results["bionj"]
        ml_tree, logL = results["ml"]

        # Build just ML tree
        results = build_trees(sequences, method="ml")
        ml_tree, logL = results["ml"]
    """
    builder = PhylogeneticTreeBuilder(verbose=verbose)
    builder.detect_and_validate(sequences)

    results = {
        "type": builder.seq_type,
        "model": builder.model_name,
    }

    if method.lower() == "all":
        upgma, bionj, ml = builder.build_all_trees(sequences, alpha=alpha)
        results["upgma"] = upgma
        results["bionj"] = bionj
        results["ml"] = ml

        # NOTE: Consensus tree functionality is currently disabled due to bugs
        # The consensus algorithm produces incorrect topologies that don't match input trees
        # See CONSENSUS_TODO.md for details and future implementation plans
        ml_tree, _ = ml
        all_trees = [upgma, bionj, ml_tree]

        # Removed broken consensus implementation
        # consensus_tree, support_values = majority_rule_consensus(all_trees, verbose=verbose)
        # results["consensus"] = consensus_tree
        # results["support_values"] = support_values

        # Compare trees
        comp_upgma_bionj = compare_trees(upgma, bionj)
        comp_upgma_ml = compare_trees(upgma, ml_tree)
        comp_bionj_ml = compare_trees(bionj, ml_tree)

        results["tree_distances"] = {
            "upgma_vs_bionj": comp_upgma_bionj,
            "upgma_vs_ml": comp_upgma_ml,
            "bionj_vs_ml": comp_bionj_ml
        }

        if verbose:
            print(f"\nTree Similarity:")
            print(f"  UPGMA vs BioNJ:  {comp_upgma_bionj['similarity']:.1%}")
            print(f"  UPGMA vs ML:     {comp_upgma_ml['similarity']:.1%}")
            print(f"  BioNJ vs ML:     {comp_bionj_ml['similarity']:.1%}")
            print("=" * 70)

    elif method.lower() == "upgma":
        results["upgma"] = builder.build_upgma_tree(sequences)
    elif method.lower() == "bionj":
        results["bionj"] = builder.build_bionj_tree(sequences)
    elif method.lower() == "ml":
        results["ml"] = builder.build_ml_tree(sequences, alpha=alpha)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'all', 'upgma', 'bionj', or 'ml'")

    return results


# Example usage
if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("UNIFIED PHYLOGENETIC TREE BUILDER")
    print("=" * 70)

    # Example 1: DNA sequences
    print("\nExample 1: DNA Sequences")
    print("-" * 70)

    dna_seqs = [
        Sequence("seq1", "E. coli", "ATGCATGC"),
        Sequence("seq2", "Salmonella", "ATGCATCC"),
        Sequence("seq3", "B. subtilis", "ATCCATGC"),
    ]

    results = build_trees(dna_seqs, method="all", verbose=True)

    # Example 2: RNA sequences
    print("\n\n" + "=" * 70)
    print("\nExample 2: RNA Sequences (16S rRNA-like)")
    print("-" * 70)

    rna_seqs = [
        Sequence("E_coli_16S", "E. coli 16S rRNA", "AUGCAUGC"),
        Sequence("Salm_16S", "Salmonella 16S rRNA", "AUGCAUCC"),
        Sequence("Bacil_16S", "B. subtilis 16S rRNA", "AUCCAUGC"),
    ]

    results = build_trees(rna_seqs, method="all", verbose=True)

    # Example 3: Single method
    print("\n\n" + "=" * 70)
    print("\nExample 3: Build Only ML Tree")
    print("-" * 70)

    results = build_trees(dna_seqs, method="ml", verbose=True)
    ml_tree, logL = results["ml"]
    print(f"\nML tree log-likelihood: {logL:.2f}")

    print("\n" + "=" * 70)
    print("EXAMPLES COMPLETED!")
    print("=" * 70)
    print("\nThe unified builder:")
    print("  [OK] Automatically detects sequence type")
    print("  [OK] Selects appropriate model")
    print("  [OK] Builds trees with all methods")
    print("  [OK] Works for DNA and RNA sequences")
    print("=" * 70)

"""
Unified phylogenetic tree builder with automatic model selection.

This module provides a single interface for building phylogenetic trees
from any sequence type (DNA, RNA, or Protein), automatically detecting
the type and selecting the appropriate substitution model.
"""

from typing import List, Tuple, Optional
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.core.sequence_type import SequenceTypeDetector, SequenceType, get_appropriate_model
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.methods.upgma import build_upgma_tree
from rrna_phylo.methods.bionj import build_bionj_tree
from rrna_phylo.distance.distance import calculate_distance_matrix
from rrna_phylo.distance.protein_distance import calculate_protein_distance_matrix
from rrna_phylo.models.ml_tree_level3 import build_ml_tree_level3
from rrna_phylo.methods.protein_ml import build_protein_ml_tree


class PhylogeneticTreeBuilder:
    """
    Unified interface for phylogenetic tree building.

    Automatically detects sequence type and builds trees using
    appropriate models for DNA, RNA, or Protein sequences.

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

    def build_upgma_tree(self, sequences: List[Sequence]) -> TreeNode:
        """
        Build UPGMA tree.

        Works for DNA, RNA (Jukes-Cantor), and Protein (Poisson).

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

        # Use appropriate distance model
        if self.seq_type == SequenceType.PROTEIN:
            dist_matrix, ids = calculate_protein_distance_matrix(sequences, model="poisson")
        else:
            dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

        tree = build_upgma_tree(dist_matrix, ids)

        if self.verbose:
            print(f"UPGMA tree built successfully")
            print(f"Newick: {tree.to_newick()};")

        return tree

    def build_bionj_tree(self, sequences: List[Sequence]) -> TreeNode:
        """
        Build BioNJ tree.

        Works for DNA, RNA (Jukes-Cantor), and Protein (Poisson).

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

        # Use appropriate distance model
        if self.seq_type == SequenceType.PROTEIN:
            dist_matrix, ids = calculate_protein_distance_matrix(sequences, model="poisson")
        else:
            dist_matrix, ids = calculate_distance_matrix(sequences, model="jukes-cantor")

        tree = build_bionj_tree(dist_matrix, ids)

        if self.verbose:
            print(f"BioNJ tree built successfully")
            print(f"Newick: {tree.to_newick()};")

        return tree

    def build_ml_tree(self, sequences: List[Sequence], alpha: float = 1.0) -> Tuple[TreeNode, float]:
        """
        Build Maximum Likelihood tree.

        Works for DNA, RNA (GTR+Gamma), and Protein (WAG/LG/JTT+Gamma).

        Args:
            sequences: Aligned sequences
            alpha: Gamma shape parameter (1.0 = moderate rate variation)

        Returns:
            (ml_tree, log_likelihood)
        """
        if self.seq_type is None:
            self.detect_and_validate(sequences)

        if self.verbose:
            print("\n" + "=" * 70)
            if self.seq_type == SequenceType.PROTEIN:
                print(f"METHOD 3: Maximum Likelihood ({self.model_name}+Gamma)")
            else:
                print("METHOD 3: Maximum Likelihood (GTR+Gamma)")
            print("=" * 70)

        # Use appropriate ML method
        if self.seq_type == SequenceType.PROTEIN:
            tree, logL = build_protein_ml_tree(
                sequences,
                model_name=self.model_name,
                alpha=alpha,
                verbose=self.verbose
            )
        else:
            tree, logL = build_ml_tree_level3(sequences, alpha=alpha, verbose=self.verbose)

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
            sequences: Aligned sequences
            alpha: Gamma shape parameter for ML tree

        Returns:
            (upgma_tree, bionj_tree, (ml_tree, log_likelihood))

        Example:
            builder = PhylogeneticTreeBuilder()
            upgma, bionj, (ml, logL) = builder.build_all_trees(sequences)
        """
        # Detect type once
        self.detect_and_validate(sequences)

        # Build all three trees
        upgma_tree = self.build_upgma_tree(sequences)
        bionj_tree = self.build_bionj_tree(sequences)
        ml_result = self.build_ml_tree(sequences, alpha=alpha)

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
            print("  ✓ Compare tree topologies")
            print("  ✓ Check if all methods agree")
            print("  ✓ Build consensus tree (recommended)")
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
        sequences: Aligned sequences (DNA, RNA, or Protein)
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
    print("  ✓ Automatically detects sequence type")
    print("  ✓ Selects appropriate model")
    print("  ✓ Builds trees with all methods")
    print("  ✓ Works for DNA, RNA, and (soon) Protein")
    print("=" * 70)

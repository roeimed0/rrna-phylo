"""
SMART Phylogenetic Tree Builder with Intelligent Method Selection.

This module extends the PhylogeneticTreeBuilder to analyze datasets and
automatically recommend which methods are suitable based on:
- Dataset size
- Phylogenetic divergence
- Distance matrix saturation

Prevents timeouts by skipping BioNJ when inappropriate.
"""

from typing import List, Tuple, Optional, Dict
from rrna_phylo.core.builder import PhylogeneticTreeBuilder
from rrna_phylo.core.tree import TreeNode
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.utils.dataset_analyzer import recommend_methods, DatasetAnalysis


class SmartPhylogeneticTreeBuilder(PhylogeneticTreeBuilder):
    """
    Intelligent tree builder that analyzes datasets before building trees.

    Automatically skips BioNJ if the dataset is:
    - Too large (>100 sequences)
    - Too divergent (>15% saturated distances)
    - Has too many invalid comparisons (>10%)

    Usage:
        builder = SmartPhylogeneticTreeBuilder()
        results = builder.build_recommended_trees(sequences)
        # Only builds trees that are suitable for the dataset
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize smart tree builder.

        Args:
            verbose: Print progress and recommendations
        """
        super().__init__(verbose=verbose)
        self.dataset_analysis: Optional[DatasetAnalysis] = None

    def analyze_and_recommend(self, sequences: List[Sequence]) -> DatasetAnalysis:
        """
        Analyze dataset and recommend appropriate methods.

        Args:
            sequences: Aligned sequences

        Returns:
            DatasetAnalysis object with recommendations
        """
        # Detect sequence type first
        if self.seq_type is None:
            self.detect_and_validate(sequences)

        # Calculate distance matrix for analysis
        dist_matrix, ids = self._calculate_distance_matrix(sequences)

        # Analyze and get recommendations
        self.dataset_analysis = recommend_methods(
            sequences,
            dist_matrix=dist_matrix,
            sequence_ids=ids,
            verbose=self.verbose
        )

        return self.dataset_analysis

    def build_recommended_trees(
        self,
        sequences: List[Sequence],
        alpha: float = 1.0,
        force_all: bool = False
    ) -> Dict[str, any]:
        """
        Build only the recommended trees for this dataset.

        Args:
            sequences: Sequences (will be aligned automatically if needed)
            alpha: Gamma shape parameter for ML tree
            force_all: Force building all methods even if not recommended

        Returns:
            Dictionary with built trees:
                - "upgma": UPGMA tree (always built)
                - "bionj": BioNJ tree (only if recommended)
                - "ml": (ML tree, log-likelihood) (always built)
                - "analysis": DatasetAnalysis object
                - "methods_built": List of methods that were built
                - "methods_skipped": List of methods that were skipped

        Example:
            builder = SmartPhylogeneticTreeBuilder()
            results = builder.build_recommended_trees(sequences)

            if "bionj" in results:
                bionj = results["bionj"]
            else:
                print(f"BioNJ skipped: {results['methods_skipped']}")
        """
        # Detect type and align if needed
        self.detect_and_validate(sequences)
        aligned_seqs = self._align_sequences(sequences)

        # Analyze dataset
        analysis = self.analyze_and_recommend(aligned_seqs)

        results = {
            "analysis": analysis,
            "methods_built": [],
            "methods_skipped": []
        }

        # UPGMA: Always build (fast and robust)
        if force_all or "upgma" in analysis.recommended_methods:
            results["upgma"] = self.build_upgma_tree(aligned_seqs)
            results["methods_built"].append("upgma")
        else:
            results["methods_skipped"].append("upgma")

        # BioNJ: Only build if recommended
        if force_all or "bionj" in analysis.recommended_methods:
            if self.verbose and "bionj" not in analysis.recommended_methods:
                print("\n" + "!" * 80)
                print("WARNING: Building BioNJ despite recommendation to skip")
                print("This may take a very long time or timeout!")
                print("!" * 80)

            results["bionj"] = self.build_bionj_tree(aligned_seqs)
            results["methods_built"].append("bionj")
        else:
            if self.verbose:
                print("\n" + "=" * 80)
                print("SKIPPING BIONJ: Not suitable for this dataset")
                print("  Reason: See dataset analysis above")
                print("=" * 80)
            results["methods_skipped"].append("bionj")

        # ML: Always build (most accurate)
        # Use skip_model_selection=True for 10x speedup (model selection currently broken anyway)
        if force_all or "ml" in analysis.recommended_methods:
            results["ml"] = self.build_ml_tree(aligned_seqs, alpha=alpha, skip_model_selection=True)
            results["methods_built"].append("ml")
        else:
            results["methods_skipped"].append("ml")

        if self.verbose:
            print("\n" + "=" * 80)
            print("SMART TREE BUILDING COMPLETE")
            print("=" * 80)
            print(f"Methods built: {', '.join(m.upper() for m in results['methods_built'])}")
            if results['methods_skipped']:
                print(f"Methods skipped: {', '.join(m.upper() for m in results['methods_skipped'])}")
            print("=" * 80)

        return results


# Convenience function
def build_trees_smart(
    sequences: List[Sequence],
    alpha: float = 1.0,
    verbose: bool = True,
    force_all: bool = False
) -> Dict[str, any]:
    """
    Build phylogenetic trees with intelligent method selection.

    This function automatically analyzes the dataset and only builds
    trees using methods appropriate for the data.

    Args:
        sequences: Aligned sequences (DNA, RNA, or Protein)
        alpha: Gamma shape parameter for ML (default 1.0)
        verbose: Print progress and analysis
        force_all: Force building all methods even if not recommended

    Returns:
        Dictionary with tree(s) and analysis:
            - "upgma": UPGMA tree (if built)
            - "bionj": BioNJ tree (if built)
            - "ml": (ML tree, log-likelihood) (if built)
            - "analysis": DatasetAnalysis object
            - "methods_built": List of methods built
            - "methods_skipped": List of methods skipped

    Example:
        # Automatic method selection
        results = build_trees_smart(sequences)

        # Access trees safely
        upgma = results.get("upgma")
        bionj = results.get("bionj")  # May be None if skipped
        ml_tree, logL = results.get("ml", (None, None))

        # Check what was built
        print(f"Built: {results['methods_built']}")
        print(f"Skipped: {results['methods_skipped']}")
    """
    builder = SmartPhylogeneticTreeBuilder(verbose=verbose)
    return builder.build_recommended_trees(sequences, alpha=alpha, force_all=force_all)


# Example usage
if __name__ == "__main__":
    from rrna_phylo.io.fasta_parser import Sequence

    print("=" * 80)
    print("SMART PHYLOGENETIC TREE BUILDER - EXAMPLES")
    print("=" * 80)

    # Example 1: Small closely-related dataset (all methods suitable)
    print("\nExample 1: Small dataset (all methods recommended)")
    print("-" * 80)

    seqs1 = [
        Sequence("s1", "Ecoli1", "ATGCATGCATGC"),
        Sequence("s2", "Ecoli2", "ATGCATGCATCC"),
        Sequence("s3", "Salm1", "ATGCATCCATGC"),
        Sequence("s4", "Bacil1", "ATCCATGCATGC"),
    ]

    results1 = build_trees_smart(seqs1, verbose=True)
    print(f"\nBuilt: {results1['methods_built']}")
    print(f"Skipped: {results1['methods_skipped']}")

    # Example 2: Large divergent dataset (BioNJ skipped)
    print("\n\nExample 2: Large divergent dataset (BioNJ skipped)")
    print("-" * 80)

    # Create 100 highly divergent sequences
    import random
    bases = "ATGC"
    seqs2 = []
    for i in range(100):
        # Create random divergent sequences
        seq = "".join(random.choice(bases) for _ in range(1000))
        seqs2.append(Sequence(f"sp{i}", f"Species{i}", seq))

    results2 = build_trees_smart(seqs2, verbose=True)
    print(f"\nBuilt: {results2['methods_built']}")
    print(f"Skipped: {results2['methods_skipped']}")

    # Example 3: Force all methods
    print("\n\nExample 3: Force all methods (override recommendations)")
    print("-" * 80)

    results3 = build_trees_smart(seqs2, verbose=True, force_all=True)
    print(f"\nBuilt: {results3['methods_built']}")
    print(f"Skipped: {results3['methods_skipped']}")

    print("\n" + "=" * 80)
    print("EXAMPLES COMPLETE!")
    print("=" * 80)

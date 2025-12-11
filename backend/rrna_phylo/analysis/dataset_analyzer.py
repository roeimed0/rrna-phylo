"""Analyze datasets to recommend phylogenetic methods."""

from typing import List, Tuple, Dict, Optional
import numpy as np
from rrna_phylo.io.fasta_parser import Sequence


class DatasetAnalysis:
    """Results of dataset analysis."""

    def __init__(
        self,
        n_sequences: int,
        alignment_length: int,
        max_distance: float,
        mean_distance: float,
        median_distance: float,
        pct_saturated: float,
        pct_invalid: float,
        recommended_methods: List[str],
        warnings: List[str]
    ):
        self.n_sequences = n_sequences
        self.alignment_length = alignment_length
        self.max_distance = max_distance
        self.mean_distance = mean_distance
        self.median_distance = median_distance
        self.pct_saturated = pct_saturated
        self.pct_invalid = pct_invalid
        self.recommended_methods = recommended_methods
        self.warnings = warnings


def analyze_distance_matrix(
    dist_matrix: np.ndarray,
    sequence_ids: List[str],
    saturation_threshold: float = 0.75,
    high_divergence_threshold: float = 0.5
) -> Tuple[Dict[str, float], List[str]]:
    """Analyze distance matrix for quality and divergence patterns."""
    n = len(dist_matrix)

    # Extract upper triangle (excluding diagonal)
    triu_indices = np.triu_indices(n, k=1)
    distances = dist_matrix[triu_indices]

    # Handle NaN/inf values
    valid_distances = distances[np.isfinite(distances)]
    invalid_distances = distances[~np.isfinite(distances)]

    # Statistics
    stats = {
        'n_sequences': n,
        'n_comparisons': len(distances),
        'n_valid': len(valid_distances),
        'n_invalid': len(invalid_distances),
        'pct_invalid': 100 * len(invalid_distances) / len(distances) if len(distances) > 0 else 0,
    }

    if len(valid_distances) > 0:
        stats['min_distance'] = float(np.min(valid_distances))
        stats['max_distance'] = float(np.max(valid_distances))
        stats['mean_distance'] = float(np.mean(valid_distances))
        stats['median_distance'] = float(np.median(valid_distances))
        stats['std_distance'] = float(np.std(valid_distances))

        # Count saturated distances (too divergent for Jukes-Cantor)
        saturated = valid_distances >= saturation_threshold
        stats['n_saturated'] = int(np.sum(saturated))
        stats['pct_saturated'] = 100 * stats['n_saturated'] / len(valid_distances)

        # Count highly divergent (but not saturated)
        high_divergence = (valid_distances >= high_divergence_threshold) & (valid_distances < saturation_threshold)
        stats['n_high_divergence'] = int(np.sum(high_divergence))
        stats['pct_high_divergence'] = 100 * stats['n_high_divergence'] / len(valid_distances)
    else:
        stats['min_distance'] = np.nan
        stats['max_distance'] = np.nan
        stats['mean_distance'] = np.nan
        stats['median_distance'] = np.nan
        stats['std_distance'] = np.nan
        stats['n_saturated'] = len(distances)
        stats['pct_saturated'] = 100.0
        stats['n_high_divergence'] = 0
        stats['pct_high_divergence'] = 0.0

    # Generate warnings
    warnings = []

    if stats['pct_invalid'] > 5:
        warnings.append(
            f"WARNING: {stats['pct_invalid']:.1f}% of pairwise distances are invalid (NaN/inf). "
            f"This indicates sequences with no comparable positions."
        )

    if stats['pct_saturated'] > 20:
        warnings.append(
            f"WARNING: {stats['pct_saturated']:.1f}% of pairwise distances are saturated (p >= {saturation_threshold}). "
            f"Jukes-Cantor correction fails for such divergent sequences."
        )

    if stats['pct_high_divergence'] > 50:
        warnings.append(
            f"INFO: {stats['pct_high_divergence']:.1f}% of distances show high divergence "
            f"(p >= {high_divergence_threshold}). Consider using Maximum Likelihood instead of distance methods."
        )

    if stats.get('mean_distance', 0) > high_divergence_threshold:
        warnings.append(
            f"INFO: Mean distance is {stats['mean_distance']:.3f}, indicating phylogenetically diverse dataset. "
            f"Distance-based methods (UPGMA, BioNJ) may be less reliable than ML."
        )

    return stats, warnings


def recommend_methods(
    sequences: List[Sequence],
    dist_matrix: Optional[np.ndarray] = None,
    sequence_ids: Optional[List[str]] = None,
    max_bionj_sequences: int = 50,
    max_bionj_saturation_pct: float = 3.0,
    max_bionj_mean_distance: float = 0.4,
    verbose: bool = True
) -> DatasetAnalysis:
    """Analyze dataset and recommend appropriate phylogenetic methods."""
    n_sequences = len(sequences)
    alignment_length = sequences[0].aligned_length if sequences else 0

    # If no distance matrix provided, create one
    if dist_matrix is None:
        from rrna_phylo.distance.distance import calculate_distance_matrix
        dist_matrix, sequence_ids = calculate_distance_matrix(sequences, model="jukes-cantor")

    # Analyze distance matrix
    stats, warnings = analyze_distance_matrix(dist_matrix, sequence_ids)

    # Determine recommended methods
    recommended = []

    # UPGMA: Always recommended (fast and robust)
    recommended.append("upgma")

    # BioNJ: Only if dataset is not too large or divergent
    bionj_suitable = True
    bionj_reasons = []

    if n_sequences > max_bionj_sequences:
        bionj_suitable = False
        bionj_reasons.append(f"dataset too large (>{max_bionj_sequences} sequences)")

    if stats['pct_saturated'] > max_bionj_saturation_pct:
        bionj_suitable = False
        bionj_reasons.append(f"{stats['pct_saturated']:.1f}% saturated distances (max {max_bionj_saturation_pct}%)")

    if stats.get('mean_distance', 0) > max_bionj_mean_distance:
        bionj_suitable = False
        bionj_reasons.append(f"mean distance {stats['mean_distance']:.3f} too high (max {max_bionj_mean_distance})")

    if stats['pct_invalid'] > 10:
        bionj_suitable = False
        bionj_reasons.append(f"{stats['pct_invalid']:.1f}% invalid distances")

    if bionj_suitable:
        recommended.append("bionj")
    else:
        warnings.append(
            f"SKIP BioNJ: Not recommended for this dataset ({', '.join(bionj_reasons)}). "
            f"BioNJ may timeout or produce unreliable results."
        )

    # ML: Always recommended (most accurate, especially for divergent sequences)
    recommended.append("ml")

    # Create analysis object
    analysis = DatasetAnalysis(
        n_sequences=n_sequences,
        alignment_length=alignment_length,
        max_distance=stats.get('max_distance', np.nan),
        mean_distance=stats.get('mean_distance', np.nan),
        median_distance=stats.get('median_distance', np.nan),
        pct_saturated=stats.get('pct_saturated', 0.0),
        pct_invalid=stats.get('pct_invalid', 0.0),
        recommended_methods=recommended,
        warnings=warnings
    )

    # Print report
    if verbose:
        print("\n" + "=" * 80)
        print("DATASET ANALYSIS")
        print("=" * 80)
        print(f"\nSequences: {n_sequences}")
        print(f"Alignment length: {alignment_length} bp")
        print(f"\nDistance Matrix Statistics:")
        print(f"  Valid comparisons: {stats['n_valid']:,} / {stats['n_comparisons']:,} ({100 - stats['pct_invalid']:.1f}%)")

        if stats['n_valid'] > 0:
            print(f"  Distance range: {stats['min_distance']:.4f} - {stats['max_distance']:.4f}")
            print(f"  Mean distance: {stats['mean_distance']:.4f}")
            print(f"  Median distance: {stats['median_distance']:.4f}")
            print(f"  Saturated (p >= 0.75): {stats['n_saturated']:,} ({stats['pct_saturated']:.1f}%)")
            print(f"  High divergence (p >= 0.50): {stats['n_high_divergence']:,} ({stats['pct_high_divergence']:.1f}%)")

        print(f"\nRecommended Methods: {', '.join(m.upper() for m in recommended)}")

        if "bionj" not in recommended:
            print(f"Skipped Methods: BioNJ (not suitable for this dataset)")

        if warnings:
            print(f"\nWarnings & Recommendations:")
            for i, warning in enumerate(warnings, 1):
                print(f"  {i}. {warning}")

        print("=" * 80)

    return analysis

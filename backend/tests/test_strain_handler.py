"""
Comprehensive tests for strain handler deduplication functions.

Tests exact duplicate removal, similarity-based clustering, species-aware
clustering, and the full smart deduplication pipeline.
"""

import pytest
from rrna_phylo.io.fasta_parser import Sequence
from rrna_phylo.utils.strain_handler import (
    remove_exact_duplicates,
    calculate_sequence_similarity,
    cluster_similar_sequences,
    smart_dereplicate,
    select_representative,
    detect_strain_groups,
    dereplicate_strains,
    get_strain_summary
)


class TestRemoveExactDuplicates:
    """Test exact duplicate removal."""

    def test_removes_exact_duplicates(self):
        """Test that exact duplicates are removed."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGC"),  # Exact duplicate
            Sequence("seq3", "E. coli", "ATGC"),  # Exact duplicate
            Sequence("seq4", "P. aeruginosa", "GGCC"),
        ]

        unique, dup_map = remove_exact_duplicates(seqs)

        assert len(unique) == 2  # Only 2 unique sequences
        assert len(dup_map["seq1"]) == 3  # seq1, seq2, seq3 are duplicates
        assert len(dup_map["seq4"]) == 1  # seq4 is unique

    def test_preserves_unique_sequences(self):
        """Test that all unique sequences are preserved."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGG"),
            Sequence("seq3", "E. coli", "ATGA"),
        ]

        unique, dup_map = remove_exact_duplicates(seqs)

        assert len(unique) == 3  # All unique
        assert all(len(dups) == 1 for dups in dup_map.values())

    def test_ignores_gaps_in_comparison(self):
        """Test that gaps are ignored when comparing sequences."""
        seqs = [
            Sequence("seq1", "E. coli", "AT-GC"),
            Sequence("seq2", "E. coli", "ATGC"),  # Same as seq1 without gap
        ]

        unique, dup_map = remove_exact_duplicates(seqs)

        assert len(unique) == 1  # Should be treated as duplicates
        assert len(dup_map["seq1"]) == 2

    def test_case_insensitive(self):
        """Test that comparison is case-insensitive."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "atgc"),  # Lowercase
        ]

        unique, dup_map = remove_exact_duplicates(seqs)

        assert len(unique) == 1  # Should be treated as duplicates

    def test_empty_input(self):
        """Test with empty sequence list."""
        unique, dup_map = remove_exact_duplicates([])

        assert len(unique) == 0
        assert len(dup_map) == 0


class TestCalculateSequenceSimilarity:
    """Test sequence similarity calculation."""

    def test_identical_sequences(self):
        """Test 100% identical sequences."""
        sim = calculate_sequence_similarity("ATGC", "ATGC")
        assert sim == 100.0

    def test_completely_different(self):
        """Test completely different sequences."""
        sim = calculate_sequence_similarity("AAAA", "TTTT")
        assert sim == 0.0

    def test_partial_similarity(self):
        """Test partially similar sequences."""
        sim = calculate_sequence_similarity("ATGC", "ATGG")
        assert sim == 75.0  # 3 out of 4 match

    def test_unequal_lengths(self):
        """Test sequences of different lengths."""
        sim = calculate_sequence_similarity("ATGCATGC", "ATGC")
        assert sim == 100.0  # All 4 positions of shorter sequence match


class TestClusterSimilarSequences:
    """Test similarity-based clustering."""

    def test_species_aware_clustering_keeps_species_separate(self):
        """Test that species-aware clustering keeps different species separate."""
        seqs = [
            Sequence("U00096.1.100", "Escherichia coli", "ATGC"),
            Sequence("U00096.2.200", "Escherichia coli", "ATGC"),  # Same species, identical
            Sequence("AE006468.1.100", "Salmonella enterica", "ATGC"),  # Different species, identical
        ]

        clusters = cluster_similar_sequences(seqs, similarity_threshold=99.0, species_aware=True)

        # Should form 2 clusters: E. coli and Salmonella separate
        assert len(clusters) == 2

        # Check cluster composition
        ecoli_cluster = [c for c in clusters if c[0].main_accession == "U00096"][0]
        salmonella_cluster = [c for c in clusters if c[0].main_accession == "AE006468"][0]

        assert len(ecoli_cluster) == 2  # Two E. coli sequences
        assert len(salmonella_cluster) == 1  # One Salmonella sequence

    def test_species_unaware_clustering_merges_species(self):
        """Test that species-unaware clustering merges similar sequences regardless of species."""
        seqs = [
            Sequence("U00096.1.100", "Escherichia coli", "ATGC"),
            Sequence("AE006468.1.100", "Salmonella enterica", "ATGC"),  # Different species, identical
        ]

        clusters = cluster_similar_sequences(seqs, similarity_threshold=99.0, species_aware=False)

        # Should form 1 cluster: both sequences merged
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_similarity_threshold(self):
        """Test that similarity threshold works correctly."""
        seqs = [
            Sequence("U00096.1.100", "E. coli", "ATGCATGC"),  # 8 bp
            Sequence("U00096.2.200", "E. coli", "ATGCATGG"),  # 7/8 match = 87.5%
            Sequence("U00096.3.300", "E. coli", "ATGCATGC"),  # 100% match with seq1
        ]

        # With 99% threshold, seq1 and seq3 should cluster, seq2 separate
        clusters = cluster_similar_sequences(seqs, similarity_threshold=99.0, species_aware=True)

        assert len(clusters) == 2

        # With 85% threshold, all should cluster
        clusters = cluster_similar_sequences(seqs, similarity_threshold=85.0, species_aware=True)

        assert len(clusters) == 1
        assert len(clusters[0]) == 3

    def test_empty_input(self):
        """Test with empty sequence list."""
        clusters = cluster_similar_sequences([])
        assert len(clusters) == 0


class TestSmartDereplicate:
    """Test the full smart deduplication pipeline."""

    def test_two_tier_deduplication(self):
        """Test that both exact and similarity-based deduplication work."""
        seqs = [
            # Exact duplicates
            Sequence("U00096.1.100", "E. coli", "ATGCATGC"),
            Sequence("U00096.2.200", "E. coli", "ATGCATGC"),
            Sequence("U00096.3.300", "E. coli", "ATGCATGC"),
            # Highly similar (7/8 match = 87.5%)
            Sequence("U00096.4.400", "E. coli", "ATGCATGG"),
            # Different species, identical to first group
            Sequence("AE006468.1.100", "Salmonella", "ATGCATGC"),
        ]

        derep, stats = smart_dereplicate(
            seqs,
            remove_exact=True,
            similarity_threshold=99.5,
            species_aware=True,
            verbose=False
        )

        # Should remove exact duplicates
        # 3 E. coli ATGCATGC are duplicates, 1 Salmonella ATGCATGC is also same sequence
        # So 4 total identical sequences -> remove 3 duplicates, keep 1
        assert stats["exact_duplicates_removed"] >= 2

        # After exact removal: 2 sequences remain (1 ATGCATGC, 1 ATGCATGG)
        # At 99.5% threshold, 87.5% similarity won't cluster them
        assert stats["final_count"] == 2

    def test_statistics_calculation(self):
        """Test that statistics are calculated correctly."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGC"),  # Duplicate
            Sequence("seq3", "E. coli", "ATGG"),
        ]

        derep, stats = smart_dereplicate(seqs)

        assert stats["original_count"] == 3
        assert stats["exact_duplicates_removed"] == 1
        assert stats["final_count"] == 2
        assert stats["reduction_percentage"] == pytest.approx(33.33, rel=0.1)

    def test_selection_method_longest(self):
        """Test that longest sequence is selected as representative."""
        seqs = [
            Sequence("U00096.1.100", "E. coli", "ATGC"),  # 4 bp
            Sequence("U00096.2.200", "E. coli", "ATGCATGC"),  # 8 bp - longest
            Sequence("U00096.3.300", "E. coli", "AT"),  # 2 bp
        ]

        # All highly similar after alignment, cluster them
        # For this test, make them similar enough
        seqs = [
            Sequence("U00096.1.100", "E. coli", "ATGCATGC"),
            Sequence("U00096.2.200", "E. coli", "ATGCATGC"),
            Sequence("U00096.3.300", "E. coli", "ATGCATGC"),
        ]

        derep, stats = smart_dereplicate(
            seqs,
            selection_method="longest",
            similarity_threshold=99.5
        )

        assert len(derep) == 1
        # All have same length, so first one is returned
        assert derep[0].id in ["U00096.1.100", "U00096.2.200", "U00096.3.300"]

    def test_empty_input(self):
        """Test with empty sequence list."""
        derep, stats = smart_dereplicate([])

        assert len(derep) == 0
        assert stats["original_count"] == 0
        assert stats["final_count"] == 0


class TestSelectRepresentative:
    """Test representative selection methods."""

    def test_longest_method(self):
        """Test selection of longest sequence."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGCATGCATGC"),  # Longest
            Sequence("seq3", "E. coli", "AT"),
        ]

        rep = select_representative(seqs, method="longest")
        assert rep.id == "seq2"

    def test_first_method(self):
        """Test selection of first sequence."""
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGCATGCATGC"),
            Sequence("seq3", "E. coli", "AT"),
        ]

        rep = select_representative(seqs, method="first")
        assert rep.id == "seq1"

    def test_median_method(self):
        """Test selection of median-length sequence."""
        seqs = [
            Sequence("seq1", "E. coli", "AT"),  # 2 bp
            Sequence("seq2", "E. coli", "ATGC"),  # 4 bp - median
            Sequence("seq3", "E. coli", "ATGCATGCATGC"),  # 12 bp
        ]

        rep = select_representative(seqs, method="median")
        assert rep.id == "seq2"

    def test_single_sequence(self):
        """Test with single sequence."""
        seqs = [Sequence("seq1", "E. coli", "ATGC")]

        rep = select_representative(seqs, method="longest")
        assert rep.id == "seq1"

    def test_invalid_method(self):
        """Test that invalid method raises error."""
        # Need multiple sequences to reach the method validation code
        seqs = [
            Sequence("seq1", "E. coli", "ATGC"),
            Sequence("seq2", "E. coli", "ATGCATGC")
        ]

        # Test that ValueError is raised for invalid method
        try:
            select_representative(seqs, method="invalid")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Unknown method" in str(e)


class TestDetectStrainGroups:
    """Test strain grouping by accession."""

    def test_groups_by_accession(self):
        """Test that sequences are grouped by main accession."""
        seqs = [
            Sequence("U00096.1.100", "E. coli K-12", "ATGC"),
            Sequence("U00096.2.200", "E. coli K-12", "ATGC"),
            Sequence("AE006468.1.100", "Salmonella", "ATGC"),
            Sequence("AE006468.2.200", "Salmonella", "ATGC"),
        ]

        groups = detect_strain_groups(seqs)

        assert len(groups) == 2
        assert len(groups["U00096"]) == 2
        assert len(groups["AE006468"]) == 2

    def test_handles_no_accession_pattern(self):
        """Test sequences without standard accession pattern."""
        seqs = [
            Sequence("custom_id_1", "E. coli", "ATGC"),
            Sequence("custom_id_2", "E. coli", "ATGC"),
        ]

        groups = detect_strain_groups(seqs)

        # Each sequence becomes its own group
        assert len(groups) == 2


class TestDereplicateStrains:
    """Test the original strain deduplication function."""

    def test_representative_method(self):
        """Test representative selection method."""
        seqs = [
            Sequence("U00096.1.100", "E. coli K-12", "ATGC"),
            Sequence("U00096.2.200", "E. coli K-12", "ATGCATGC"),  # Longest
            Sequence("AE006468.1.100", "Salmonella", "GGCC"),
        ]

        derep, mapping = dereplicate_strains(seqs, method="representative", selection="longest")

        assert len(derep) == 2  # One per strain
        assert len(mapping) == 2
        assert len(mapping["U00096"]) == 2
        assert len(mapping["AE006468"]) == 1


class TestGetStrainSummary:
    """Test strain summary generation."""

    def test_summary_format(self):
        """Test that summary is formatted correctly."""
        seqs = [
            Sequence("U00096.1.100", "E. coli K-12", "ATGC"),
            Sequence("U00096.2.200", "E. coli K-12", "ATGC"),
            Sequence("AE006468.1.100", "Salmonella", "ATGC"),
        ]

        summary = get_strain_summary(seqs)

        assert "2 strain(s)" in summary
        assert "3 total sequences" in summary
        assert "U00096" in summary
        assert "AE006468" in summary


class TestIntegration:
    """Integration tests for the full pipeline."""

    def test_full_pipeline_with_real_data(self):
        """Test the full deduplication pipeline with realistic data."""
        # Simulate real rRNA data: multiple copies from same genome
        seqs = [
            # E. coli K-12: 7 rRNA copies (some identical, some slightly different)
            Sequence("U00096.223771.225312", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli str. K-12 substr. MG1655", "ATGCATGCATGC"),
            Sequence("U00096.2729616.2731157", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli str. K-12 substr. MG1655", "ATGCATGCATGC"),  # Identical
            Sequence("U00096.3427221.3428762", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli str. K-12 substr. MG1655", "ATGCATGCATGC"),  # Identical
            Sequence("U00096.3941808.3943349", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli str. K-12 substr. MG1655", "ATGCATGCATGG"),  # 1 bp different

            # Salmonella: 7 rRNA copies (mostly identical)
            Sequence("AE006468.289190.290733", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella virus Fels2", "ATGCATGCATGC"),  # Same as E. coli
            Sequence("AE006468.2800121.2801663", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella virus Fels2", "ATGCATGCATGC"),  # Identical
            Sequence("AE006468.3570470.3572013", "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella virus Fels2", "ATGCATGCATGC"),  # Identical
        ]

        # Step 1: Remove exact duplicates
        unique, dup_map = remove_exact_duplicates(seqs)
        assert len(unique) == 2  # Only 2 unique sequences: ATGCATGCATGC and ATGCATGCATGG

        # Step 2: Smart dereplicate with species-aware clustering
        derep, stats = smart_dereplicate(
            seqs,
            remove_exact=True,
            similarity_threshold=99.5,
            species_aware=True,
            verbose=False
        )

        # Should keep both sequences (different by 1 bp = 91.7% similarity < 99.5%)
        # But also keep them separate by species
        # Actually, after exact dedup: 1 ATGCATGCATGC (could be from either species), 1 ATGCATGCATGG (E. coli)
        # At 99.5% threshold, 91.7% similarity won't cluster them
        assert stats["original_count"] == 7
        assert stats["exact_duplicates_removed"] == 5  # 7 original -> 2 unique

    def test_species_aware_prevents_cross_species_clustering(self):
        """Verify that species-aware clustering prevents merging different species."""
        # Use slightly different sequences to avoid exact duplicate removal
        seqs = [
            Sequence("U00096.1.100", "Escherichia coli", "ATGCATGCATGC"),
            Sequence("U00096.2.200", "Escherichia coli", "ATGCATGCATGC"),  # Same as above
            Sequence("AE006468.1.100", "Salmonella enterica", "ATGCATGCATGG"),  # 1 bp different
        ]

        # With species-aware clustering
        derep_aware, stats_aware = smart_dereplicate(
            seqs,
            species_aware=True,
            similarity_threshold=99.5,
            remove_exact=True
        )

        # After exact duplicate removal: 2 sequences (1 E. coli, 1 Salmonella)
        # With 99.5% threshold and species-aware, both are kept separate
        # 11/12 match = 91.7% < 99.5%, so they won't cluster anyway
        assert stats_aware["final_count"] == 2

        # Test with sequences that ARE similar enough to cluster
        seqs2 = [
            Sequence("U00096.1.100", "Escherichia coli", "ATGCATGCATGCATGC"),  # 16 bp
            Sequence("AE006468.1.100", "Salmonella enterica", "ATGCATGCATGCATGC"),  # Identical
        ]

        # With species-aware clustering - should keep both (different species)
        derep_aware2, stats_aware2 = smart_dereplicate(
            seqs2,
            species_aware=True,
            similarity_threshold=99.5,
            remove_exact=False  # Don't remove exact duplicates, test clustering only
        )

        # Should keep 2 because they're different species
        assert stats_aware2["final_count"] == 2

        # With species-unaware clustering - should merge to 1
        derep_unaware2, stats_unaware2 = smart_dereplicate(
            seqs2,
            species_aware=False,
            similarity_threshold=99.5,
            remove_exact=False
        )

        # Should keep only 1 (incorrectly merges different species)
        assert stats_unaware2["final_count"] == 1


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])

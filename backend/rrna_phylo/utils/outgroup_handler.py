"""
Outgroup detection and tree rooting utilities.

An outgroup is a sequence or group of sequences that diverged earlier than
the ingroup being studied. It's used to root the tree and determine the
direction of evolution.
"""

from typing import List, Optional, Tuple
from rrna_phylo.io.fasta_parser import Sequence


# Common outgroup species for bacterial phylogenetics
BACTERIAL_OUTGROUPS = {
    # For Proteobacteria studies
    "proteobacteria": [
        "Bacillus",
        "Staphylococcus",
        "Clostridium",
        "Listeria"
    ],
    # For Firmicutes studies
    "firmicutes": [
        "Escherichia",
        "Pseudomonas",
        "Salmonella",
        "Vibrio"
    ],
    # For Enterobacteriaceae studies (narrow)
    "enterobacteriaceae": [
        "Pseudomonas",
        "Acinetobacter",
        "Vibrio",
        "Bacillus"
    ],
    # Universal bacterial outgroup (Archaea)
    "archaea": [
        "Methanocaldococcus",
        "Sulfolobus",
        "Halobacterium",
        "Thermococcus"
    ]
}


def detect_outgroup(sequences: List[Sequence], target_group: str = "auto") -> Optional[List[Sequence]]:
    """
    Detect potential outgroup sequences from the dataset.

    Args:
        sequences: List of all sequences
        target_group: Target group being studied ("auto" for automatic detection)

    Returns:
        List of outgroup sequences, or None if no outgroup detected

    Example:
        >>> seqs = [
        ...     Sequence("ecoli1", "Escherichia coli", "ATGC"),
        ...     Sequence("sal1", "Salmonella enterica", "ATGC"),
        ...     Sequence("bacillus1", "Bacillus subtilis", "ATGC")  # Outgroup
        ... ]
        >>> outgroup = detect_outgroup(seqs, "enterobacteriaceae")
        >>> len(outgroup)
        1  # Bacillus is the outgroup
    """
    if target_group == "auto":
        target_group = _infer_target_group(sequences)

    if target_group not in BACTERIAL_OUTGROUPS:
        return None

    outgroup_genera = BACTERIAL_OUTGROUPS[target_group]
    outgroup_seqs = []

    for seq in sequences:
        # Check if sequence description contains outgroup genus
        desc_lower = (seq.description or "").lower()
        for genus in outgroup_genera:
            if genus.lower() in desc_lower:
                outgroup_seqs.append(seq)
                break

    return outgroup_seqs if outgroup_seqs else None


def _infer_target_group(sequences: List[Sequence]) -> str:
    """
    Infer the target taxonomic group from sequence descriptions.

    Args:
        sequences: List of sequences

    Returns:
        Inferred target group name
    """
    descriptions = [seq.description.lower() if seq.description else "" for seq in sequences]

    # Count how many sequences are from each major group
    proteobacteria_count = sum(
        1 for d in descriptions
        if any(genus in d for genus in ["escherichia", "salmonella", "pseudomonas", "vibrio"])
    )

    firmicutes_count = sum(
        1 for d in descriptions
        if any(genus in d for genus in ["bacillus", "staphylococcus", "clostridium", "listeria"])
    )

    enterobacteriaceae_count = sum(
        1 for d in descriptions
        if any(genus in d for genus in ["escherichia", "salmonella", "shigella", "klebsiella"])
    )

    # Determine predominant group
    if enterobacteriaceae_count >= len(sequences) * 0.5:
        return "enterobacteriaceae"
    elif proteobacteria_count >= len(sequences) * 0.5:
        return "proteobacteria"
    elif firmicutes_count >= len(sequences) * 0.5:
        return "firmicutes"
    else:
        return "proteobacteria"  # Default


def validate_outgroup(sequences: List[Sequence], outgroup_ids: List[str]) -> Tuple[bool, str]:
    """
    Validate that specified outgroup sequences exist and are appropriate.

    Args:
        sequences: List of all sequences
        outgroup_ids: List of sequence IDs to use as outgroup

    Returns:
        Tuple of (is_valid, message)

    Example:
        >>> valid, msg = validate_outgroup(seqs, ["bacillus1"])
        >>> valid
        True
        >>> msg
        'Outgroup validated: 1 sequence(s)'
    """
    # Check all outgroup IDs exist
    seq_ids = {seq.id for seq in sequences}
    missing_ids = set(outgroup_ids) - seq_ids

    if missing_ids:
        return False, f"Outgroup sequence(s) not found: {', '.join(missing_ids)}"

    # Check we have at least one ingroup sequence
    if len(outgroup_ids) >= len(sequences):
        return False, "Outgroup cannot include all sequences"

    # Check we have at least 2 sequences remaining for ingroup
    if len(sequences) - len(outgroup_ids) < 2:
        return False, "Need at least 2 ingroup sequences"

    return True, f"Outgroup validated: {len(outgroup_ids)} sequence(s)"


def suggest_outgroup(sequences: List[Sequence]) -> str:
    """
    Suggest appropriate outgroup selection based on sequences.

    Args:
        sequences: List of all sequences

    Returns:
        Human-readable suggestion message

    Example:
        >>> suggestion = suggest_outgroup(sequences)
        >>> print(suggestion)
        Suggested outgroup for Enterobacteriaceae study:
          - Pseudomonas (4 sequences available)
          - Use: --outgroup "AE004091.*"
    """
    target_group = _infer_target_group(sequences)
    detected_outgroup = detect_outgroup(sequences, target_group)

    if not detected_outgroup:
        return (
            f"No outgroup detected for {target_group} study.\n"
            f"Consider adding sequences from: {', '.join(BACTERIAL_OUTGROUPS[target_group])}"
        )

    # Group outgroup sequences by genus
    genus_groups = {}
    for seq in detected_outgroup:
        desc = seq.description or seq.id
        # Extract genus (first word after species name patterns)
        genus = None
        for candidate in BACTERIAL_OUTGROUPS[target_group]:
            if candidate.lower() in desc.lower():
                genus = candidate
                break

        if genus:
            if genus not in genus_groups:
                genus_groups[genus] = []
            genus_groups[genus].append(seq)

    suggestions = [f"Suggested outgroup for {target_group} study:"]
    for genus, seqs in genus_groups.items():
        # Get common accession pattern
        example_id = seqs[0].id
        accession = example_id.split('.')[0] if '.' in example_id else example_id
        suggestions.append(f"  - {genus} ({len(seqs)} sequence(s))")
        suggestions.append(f"    Use: --outgroup \"{accession}.*\"")

    return '\n'.join(suggestions)


def get_outgroup_sequences(sequences: List[Sequence], outgroup_pattern: str) -> List[Sequence]:
    """
    Extract outgroup sequences matching a pattern.

    Args:
        sequences: List of all sequences
        outgroup_pattern: Pattern to match (e.g., "AJ276351.*" or "Bacillus")

    Returns:
        List of matching outgroup sequences

    Example:
        >>> outgroup = get_outgroup_sequences(seqs, "AJ276351.*")
        >>> # Returns all sequences with IDs starting with "AJ276351"
    """
    import re

    outgroup_seqs = []

    # Convert shell-style pattern to regex
    pattern = outgroup_pattern.replace("*", ".*").replace("?", ".")
    regex = re.compile(pattern, re.IGNORECASE)

    for seq in sequences:
        # Match against ID or description
        if regex.search(seq.id) or regex.search(seq.description or ""):
            outgroup_seqs.append(seq)

    return outgroup_seqs


# Rooting function placeholder - would integrate with tree structure
def root_tree_with_outgroup(tree, outgroup_ids: List[str]):
    """
    Root a phylogenetic tree using the outgroup.

    This is a placeholder - full implementation would require tree manipulation.

    Args:
        tree: Tree object to root
        outgroup_ids: List of sequence IDs in the outgroup

    Returns:
        Rooted tree

    Note:
        Most tree viewers (FigTree, iTOL) can root trees interactively.
        For programmatic rooting, use Bio.Phylo or dendropy.
    """
    # TODO: Implement tree rooting
    # For now, document in output that user should root at outgroup
    pass


if __name__ == "__main__":
    # Example usage
    from rrna_phylo.io.fasta_parser import Sequence

    test_seqs = [
        Sequence("U00096.1", "Escherichia coli K-12", "ATGC"),
        Sequence("AE006468.1", "Salmonella enterica", "ATGC"),
        Sequence("AE006468.2", "Salmonella enterica", "ATGC"),
        Sequence("AJ276351.1", "Bacillus subtilis", "ATGC"),  # Potential outgroup
    ]

    print(suggest_outgroup(test_seqs))
    print()

    outgroup = detect_outgroup(test_seqs, "enterobacteriaceae")
    if outgroup:
        print(f"Detected outgroup: {[s.id for s in outgroup]}")

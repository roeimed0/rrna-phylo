"""
Tree validation and repair utilities.

Ensures phylogenetic trees have valid branch lengths and structure.
"""

from rrna_phylo.core.tree import TreeNode


MIN_BRANCH_LENGTH = 1e-6  # Minimum branch length (prevents zero/negative)


def validate_and_fix_branch_lengths(tree: TreeNode, min_length: float = MIN_BRANCH_LENGTH) -> tuple[TreeNode, int]:
    """
    Validate and fix branch lengths in a phylogenetic tree.

    Issues fixed:
    - Zero-length branches (< min_length) set to min_length
    - Negative branch lengths set to min_length

    Args:
        tree: Root node of the tree
        min_length: Minimum allowed branch length (default: 1e-6)

    Returns:
        Tuple of (tree, num_fixes) where tree is the fixed tree and
        num_fixes is the number of branches that were corrected

    Example:
        >>> tree = build_tree(sequences)
        >>> fixed_tree, num_fixes = validate_and_fix_branch_lengths(tree)
        >>> if num_fixes > 0:
        >>>     print(f"Fixed {num_fixes} invalid branch lengths")
    """
    num_fixes = 0

    def fix_node(node: TreeNode) -> int:
        """Recursively fix branch lengths in tree."""
        fixes = 0

        # Fix this node's distance if invalid
        if node.distance < min_length:
            if node.distance != 0.0 or not node.is_leaf():  # Don't warn for leaf initialization
                fixes += 1
            node.distance = max(min_length, node.distance)

        # Recursively fix children
        if node.left:
            fixes += fix_node(node.left)
        if node.right:
            fixes += fix_node(node.right)

        return fixes

    num_fixes = fix_node(tree)

    return tree, num_fixes


def validate_tree_structure(tree: TreeNode) -> list[str]:
    """
    Validate tree structure and report any issues.

    Checks:
    - All internal nodes have exactly 2 children (binary tree)
    - No duplicate leaf names
    - All branch lengths are non-negative
    - Tree is properly rooted

    Args:
        tree: Root node of the tree

    Returns:
        List of validation error messages (empty if valid)

    Example:
        >>> errors = validate_tree_structure(tree)
        >>> if errors:
        >>>     for err in errors:
        >>>         print(f"ERROR: {err}")
    """
    errors = []

    # Check binary structure
    def check_binary(node: TreeNode, path: str = "root") -> None:
        """Ensure all internal nodes have exactly 2 children."""
        if not node.is_leaf():
            if node.left is None or node.right is None:
                errors.append(f"Internal node at {path} has <2 children")
            else:
                check_binary(node.left, f"{path}.left")
                check_binary(node.right, f"{path}.right")

    check_binary(tree)

    # Check for duplicate leaf names
    leaves = tree.get_leaves()
    leaf_names = [leaf.name for leaf in leaves]
    seen = set()
    duplicates = set()

    for name in leaf_names:
        if name in seen:
            duplicates.add(name)
        seen.add(name)

    if duplicates:
        errors.append(f"Duplicate leaf names found: {', '.join(sorted(duplicates))}")

    # Check all branch lengths are non-negative
    all_nodes = tree.get_all_nodes()
    negative_branches = []

    for node in all_nodes:
        if node.distance < 0:
            name = node.name if node.is_leaf() else "internal"
            negative_branches.append(f"{name}:{node.distance:.6f}")

    if negative_branches:
        errors.append(f"Negative branch lengths: {', '.join(negative_branches)}")

    # Check tree is rooted (root should have 2 children if >2 taxa)
    if tree.count_leaves() > 2:
        if tree.left is None or tree.right is None:
            errors.append("Tree is not properly rooted (root needs 2 children)")

    return errors


def get_branch_length_stats(tree: TreeNode) -> dict:
    """
    Get statistics about branch lengths in the tree.

    Returns:
        Dictionary with keys:
        - min: Minimum branch length
        - max: Maximum branch length
        - mean: Mean branch length
        - median: Median branch length
        - num_zero: Number of zero-length branches (< 1e-6)
        - num_negative: Number of negative branches

    Example:
        >>> stats = get_branch_length_stats(tree)
        >>> print(f"Branch length range: {stats['min']:.6f} - {stats['max']:.6f}")
        >>> print(f"Zero-length branches: {stats['num_zero']}")
    """
    all_nodes = tree.get_all_nodes()
    branch_lengths = [node.distance for node in all_nodes]

    if not branch_lengths:
        return {
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "median": 0.0,
            "num_zero": 0,
            "num_negative": 0
        }

    branch_lengths_sorted = sorted(branch_lengths)
    n = len(branch_lengths)

    return {
        "min": min(branch_lengths),
        "max": max(branch_lengths),
        "mean": sum(branch_lengths) / n,
        "median": branch_lengths_sorted[n // 2],
        "num_zero": sum(1 for bl in branch_lengths if bl < 1e-6),
        "num_negative": sum(1 for bl in branch_lengths if bl < 0)
    }

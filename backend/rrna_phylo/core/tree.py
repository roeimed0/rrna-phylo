"""
Core tree data structure for phylogenetic trees.

This module defines the TreeNode class used throughout the package.
"""


class TreeNode:
    """Represents a node in a phylogenetic tree."""

    def __init__(self, name: str = None, left=None, right=None, distance: float = 0.0, support: float = None):
        """Initialize tree node."""
        self.name = name
        self.left = left
        self.right = right
        self.distance = distance
        self.support = support  # NEW: Support value for consensus/bootstrap
        self.height = 0.0  # Used during tree construction

    def is_leaf(self) -> bool:
        """Check if this is a leaf node."""
        return self.left is None and self.right is None

    def to_newick(self, include_internal_names: bool = False) -> str:
        """Convert tree to Newick format string."""
        if self.is_leaf():
            # Quote name if it contains spaces or special characters
            quoted_name = f"'{self.name}'" if (' ' in self.name or '(' in self.name or ')' in self.name) else self.name
            return f"{quoted_name}:{self.distance:.6f}"
        else:
            left_str = self.left.to_newick(include_internal_names=include_internal_names)
            right_str = self.right.to_newick(include_internal_names=include_internal_names)
            # Only include internal node name if requested and present
            if include_internal_names and self.name:
                quoted_name = f"'{self.name}'" if (' ' in self.name or '(' in self.name or ')' in self.name) else self.name
                return f"({left_str},{right_str}){quoted_name}:{self.distance:.6f}"
            else:
                return f"({left_str},{right_str}):{self.distance:.6f}"

    def __repr__(self):
        """String representation for debugging."""
        if self.is_leaf():
            return f"Leaf({self.name})"
        else:
            left_name = self.left.name if self.left and self.left.is_leaf() else "Internal"
            right_name = self.right.name if self.right and self.right.is_leaf() else "Internal"
            return f"Internal(left={left_name}, right={right_name})"

    def count_leaves(self) -> int:
        """Count the number of leaf nodes (taxa) in the tree."""
        if self.is_leaf():
            return 1
        count = 0
        if self.left:
            count += self.left.count_leaves()
        if self.right:
            count += self.right.count_leaves()
        return count

    def get_leaves(self) -> list:
        """Get all leaf nodes in the tree."""
        if self.is_leaf():
            return [self]
        leaves = []
        if self.left:
            leaves.extend(self.left.get_leaves())
        if self.right:
            leaves.extend(self.right.get_leaves())
        return leaves

    def get_leaf_names(self) -> list:
        """Get names of all leaves in the tree."""
        return [leaf.name for leaf in self.get_leaves()]

    def copy(self):
        """Create a deep copy of the tree."""
        if self.is_leaf():
            return TreeNode(
                name=self.name,
                distance=self.distance,
                support=self.support
            )

        left_copy = self.left.copy() if self.left else None
        right_copy = self.right.copy() if self.right else None

        return TreeNode(
            name=self.name,
            left=left_copy,
            right=right_copy,
            distance=self.distance,
            support=self.support
        )

    def find_node(self, name: str):
        """Find a node by name in the tree."""
        if self.name == name:
            return self

        if self.left:
            result = self.left.find_node(name)
            if result:
                return result

        if self.right:
            result = self.right.find_node(name)
            if result:
                return result

        return None

    def get_internal_nodes(self) -> list:
        """Get all internal (non-leaf) nodes in the tree."""
        if self.is_leaf():
            return []

        nodes = [self]
        if self.left:
            nodes.extend(self.left.get_internal_nodes())
        if self.right:
            nodes.extend(self.right.get_internal_nodes())

        return nodes

    def get_all_nodes(self) -> list:
        """Get all nodes (internal and leaf) in the tree."""
        if self.is_leaf():
            return [self]

        nodes = [self]
        if self.left:
            nodes.extend(self.left.get_all_nodes())
        if self.right:
            nodes.extend(self.right.get_all_nodes())

        return nodes

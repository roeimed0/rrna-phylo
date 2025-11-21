"""
Core tree data structure for phylogenetic trees.

This module defines the TreeNode class used throughout the package.
"""


class TreeNode:
    """
    Represents a node in a phylogenetic tree.

    A phylogenetic tree is a binary tree where:
    - Leaf nodes represent sequences (taxa)
    - Internal nodes represent ancestral sequences
    - Branch lengths represent evolutionary distance

    Attributes:
        name: Name of the node (used for leaf nodes)
        left: Left child node
        right: Right child node
        distance: Branch length from parent to this node
        height: Height of node in tree (used during construction)
    """

    def __init__(self, name: str = None, left=None, right=None, distance: float = 0.0):
        """
        Initialize tree node.

        Args:
            name: Name of the node (for leaves, typically sequence ID)
            left: Left child TreeNode
            right: Right child TreeNode
            distance: Branch length from parent (evolutionary distance)
        """
        self.name = name
        self.left = left
        self.right = right
        self.distance = distance
        self.height = 0.0  # Used during tree construction

    def is_leaf(self) -> bool:
        """
        Check if this is a leaf node (terminal node representing a sequence).

        Returns:
            True if leaf node (no children), False if internal node
        """
        return self.left is None and self.right is None

    def to_newick(self) -> str:
        """
        Convert tree to Newick format string.

        Newick format is the standard for representing phylogenetic trees:
        - Leaves: name:distance
        - Internal: (left,right):distance
        - Example: ((A:0.1,B:0.2):0.3,C:0.4);

        Returns:
            Newick format string (without trailing semicolon)

        Example:
            >>> leaf = TreeNode("A", distance=0.1)
            >>> leaf.to_newick()
            'A:0.100000'
        """
        if self.is_leaf():
            return f"{self.name}:{self.distance:.6f}"
        else:
            left_str = self.left.to_newick()
            right_str = self.right.to_newick()
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
        """
        Count the number of leaf nodes (taxa) in the tree.

        Returns:
            Number of leaves in the subtree rooted at this node
        """
        if self.is_leaf():
            return 1
        count = 0
        if self.left:
            count += self.left.count_leaves()
        if self.right:
            count += self.right.count_leaves()
        return count

    def get_leaves(self) -> list:
        """
        Get all leaf nodes in the tree.

        Returns:
            List of leaf TreeNode objects
        """
        if self.is_leaf():
            return [self]
        leaves = []
        if self.left:
            leaves.extend(self.left.get_leaves())
        if self.right:
            leaves.extend(self.right.get_leaves())
        return leaves

    def get_leaf_names(self) -> list:
        """
        Get names of all leaves in the tree.

        Returns:
            List of leaf names (sequence IDs)
        """
        return [leaf.name for leaf in self.get_leaves()]

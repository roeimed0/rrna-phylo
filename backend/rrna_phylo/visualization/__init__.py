"""
Phylogenetic tree visualization using ETE3.

This module provides publication-quality tree visualization with support for:
- Circular and rectangular layouts
- Bootstrap value display
- Branch length annotations
- Custom node styling
- High-resolution output (PDF, SVG, PNG)
"""

try:
    from .ete3_viz import visualize_tree, ETE3_AVAILABLE
except ImportError:
    ETE3_AVAILABLE = False

    def visualize_tree(*args, **kwargs):
        raise ImportError(
            "ETE3 not installed. Install with: pip install ete3\n"
            "Note: May require Qt dependencies on some systems."
        )

__all__ = ['visualize_tree', 'ETE3_AVAILABLE']

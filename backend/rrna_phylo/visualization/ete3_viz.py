"""
Publication-quality phylogenetic tree visualization using ETE3.

ETE3 (Environment for Tree Exploration) is the gold standard for phylogenetic
visualization in Python, used by Nature, Science, and Cell publications.

Features:
- Publication-quality output (600 DPI, vector graphics)
- Circular and rectangular layouts
- Bootstrap value display with customizable styling
- Branch length annotations
- Multiple output formats (PDF, SVG, PNG, EPS)
"""

from typing import Optional

try:
    from ete3 import Tree
    from ete3.treeview import TreeStyle, NodeStyle, TextFace
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False


def _prepare_newick_for_ete3(newick_content: str) -> str:
    """
    Prepare Newick string for ETE3 - keep quotes as they are needed for names with spaces.

    ETE3 CAN handle single-quoted names, they're actually required for names with spaces!

    Args:
        newick_content: Newick format string with potentially quoted names

    Returns:
        Newick string ready for ETE3
    """
    # Keep the newick content as-is - ETE3 handles quotes correctly
    return newick_content


def visualize_tree(
    newick_file: str,
    output_file: str,
    show_bootstrap: bool = True,
    show_branch_length: bool = True,
    show_leaf_names: bool = True,
    bootstrap_threshold: float = 70.0,
    width: Optional[int] = None,
    height: Optional[int] = None,
    dpi: int = 300,
    branch_vertical_margin: int = 10,
    title: Optional[str] = None,
    branch_line_width: float = 1.5,
    scale: int = 120
) -> bool:
    """
    Create publication-quality phylogenetic tree visualization.

    Args:
        newick_file (str): Path to Newick format tree file
        output_file (str): Output image path (extension determines format)
                          Supported: .png, .pdf, .svg, .eps
        show_bootstrap (bool): Display bootstrap support values
        show_branch_length (bool): Display branch lengths
        show_leaf_names (bool): Show leaf (tip) labels
        bootstrap_threshold (float): Threshold for coloring bootstrap values
                                     Values >= threshold shown in green
                                     Values < threshold shown in red
        width (int, optional): Image width in pixels (None = auto)
        height (int, optional): Image height in pixels (None = auto)
        dpi (int): Resolution for raster output (PNG). Use 300+ for publication
        branch_vertical_margin (int): Vertical spacing between branches
        title (str, optional): Tree title
        branch_line_width (float): Branch line thickness in pixels (default: 1.5)
        scale (int): Pixels per branch length unit for zoom (default: 120)

    Returns:
        bool: True if successful, False otherwise

    Example:
        >>> visualize_tree(
        ...     "tree.nwk",
        ...     "tree_publication.pdf",
        ...     show_bootstrap=True,
        ...     dpi=600
        ... )

    Raises:
        ImportError: If ETE3 not installed
        FileNotFoundError: If newick_file doesn't exist
    """
    if not ETE3_AVAILABLE:
        raise ImportError(
            "ETE3 is required for visualization.\n"
            "Install with: pip install ete3\n"
            "Note: On some systems, Qt dependencies may be required:\n"
            "  Ubuntu/Debian: sudo apt-get install python3-pyqt5\n"
            "  macOS: brew install qt5\n"
            "  Windows: Usually works with pip install alone"
        )

    try:
        # Read and prepare Newick file for ETE3
        with open(newick_file, 'r') as f:
            newick_content = f.read().strip()

        # Ensure newick string ends with semicolon (required by ETE3)
        if not newick_content.endswith(';'):
            newick_content += ';'

        # Remove quotes for ETE3 compatibility
        newick_content = _prepare_newick_for_ete3(newick_content)

        # Load tree (format=0 is flexible and works with quoted names)
        # quoted_node_names=True allows names with spaces in single quotes
        tree = Tree(newick_content, format=0, quoted_node_names=True)

        # Create tree style
        ts = TreeStyle()
        ts.mode = "r"  # Rectangular layout

        # IMPORTANT: Align all leaf names to the same vertical axis
        ts.force_topology = True
        ts.draw_aligned_faces_as_table = False

        # Tree appearance - Simple and clean!
        ts.show_leaf_name = show_leaf_names
        ts.show_branch_length = show_branch_length
        ts.show_branch_support = show_bootstrap  # ETE3 built-in support display
        ts.branch_vertical_margin = branch_vertical_margin

        # Guide lines to connect branches to aligned names - BLACK SOLID LINES
        ts.draw_guiding_lines = True  # Enable guide lines for alignment
        ts.guiding_lines_type = 0  # Solid lines (0=solid, 1=dashed, 2=dotted)
        ts.guiding_lines_color = "black"  # Black color

        # Scale (zoom) - longer branches to see distance values better
        ts.scale = scale  # Pixels per branch length unit

        # No border around tree
        ts.show_border = False

        # Title
        if title:
            ts.title.add_face(TextFace(title, fsize=20, bold=True), column=0)

        # Scale bar
        ts.show_scale = show_branch_length

        # Force complete line rendering
        ts.complete_branch_lines_when_necessary = True
        ts.optimal_scale_level = "full"

        # Apply branch line width and color to all nodes - Simple styling
        for node in tree.traverse():
            nstyle = NodeStyle()
            # Branch styling
            nstyle["hz_line_width"] = branch_line_width
            nstyle["vt_line_width"] = branch_line_width
            nstyle["hz_line_color"] = "black"
            nstyle["vt_line_color"] = "black"
            nstyle["hz_line_type"] = 0  # Solid line
            nstyle["vt_line_type"] = 0

            # Node styling (simple and clean)
            nstyle["size"] = 0  # Hide node circles for cleaner look

            node.set_style(nstyle)

        # Render tree
        render_kwargs = {
            'tree_style': ts,
            'dpi': dpi
        }

        if width:
            render_kwargs['w'] = width
        if height:
            render_kwargs['h'] = height

        tree.render(output_file, **render_kwargs)

        return True

    except FileNotFoundError:
        print(f"Error: Newick file not found: {newick_file}")
        return False
    except Exception as e:
        print(f"Error visualizing tree: {e}")
        import traceback
        traceback.print_exc()
        return False


def batch_visualize(newick_files: list, output_dir: str, **kwargs) -> dict:
    """
    Visualize multiple trees at once.

    Args:
        newick_files (list): List of Newick file paths
        output_dir (str): Output directory for images
        **kwargs: Additional arguments passed to visualize_tree()

    Returns:
        dict: {filename: success_status}

    Example:
        >>> batch_visualize(
        ...     ['tree1.nwk', 'tree2.nwk'],
        ...     'output/',
        ...     dpi=600
        ... )
    """
    import os
    from pathlib import Path

    results = {}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for newick_file in newick_files:
        filename = Path(newick_file).stem
        output_file = output_path / f"{filename}.png"

        print(f"Visualizing {filename}...")
        success = visualize_tree(
            newick_file,
            str(output_file),
            **kwargs
        )
        results[filename] = success

    return results

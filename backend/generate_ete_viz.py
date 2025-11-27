"""
Generate ETE3 visualizations from test output trees.
"""
from ete3 import Tree, TreeStyle, NodeStyle
import sys

def visualize_tree(newick_file, output_file, title="Phylogenetic Tree"):
    """Generate ETE3 visualization from Newick file."""
    try:
        # Load tree with quoted_node_names flag for names with spaces
        tree = Tree(newick_file, format=1, quoted_node_names=True)

        # Configure tree style
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "c"  # Circular mode
        ts.title.add_face(title, column=0)
        ts.show_scale = True
        ts.scale = 100

        # Style nodes
        for node in tree.traverse():
            if node.is_leaf():
                nstyle = NodeStyle()
                nstyle["size"] = 3
                nstyle["fgcolor"] = "blue"
                node.set_style(nstyle)
            else:
                nstyle = NodeStyle()
                nstyle["size"] = 2
                nstyle["fgcolor"] = "black"
                node.set_style(nstyle)

        # Render to file
        tree.render(output_file, w=1200, h=1200, tree_style=ts, dpi=150)
        print(f"[OK] Saved visualization: {output_file}")
        return True

    except Exception as e:
        print(f"[ERROR] Failed to visualize {newick_file}: {e}")
        return False


if __name__ == "__main__":
    import os

    output_dir = "test_100species_smart_output"

    # UPGMA tree
    upgma_nwk = os.path.join(output_dir, "upgma_smart.nwk")
    upgma_pdf = os.path.join(output_dir, "upgma_smart_tree.pdf")
    upgma_png = os.path.join(output_dir, "upgma_smart_tree.png")

    if os.path.exists(upgma_nwk):
        print(f"Visualizing UPGMA tree...")
        visualize_tree(upgma_nwk, upgma_pdf, "UPGMA Tree (87 sequences)")
        visualize_tree(upgma_nwk, upgma_png, "UPGMA Tree (87 sequences)")

    # ML tree
    ml_nwk = os.path.join(output_dir, "ml_smart.nwk")
    ml_pdf = os.path.join(output_dir, "ml_smart_tree.pdf")
    ml_png = os.path.join(output_dir, "ml_smart_tree.png")

    if os.path.exists(ml_nwk):
        print(f"Visualizing ML tree...")
        visualize_tree(ml_nwk, ml_pdf, "ML Tree (GTR+G, 87 sequences)")
        visualize_tree(ml_nwk, ml_png, "ML Tree (GTR+G, 87 sequences)")

    print("\nAll visualizations complete!")

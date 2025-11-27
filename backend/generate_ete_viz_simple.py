"""
Generate ETE3 visualizations from test output trees using existing module.
"""
from rrna_phylo.visualization.ete3_viz import visualize_tree
import os

output_dir = "test_100species_smart_output"

# UPGMA tree
upgma_nwk = os.path.join(output_dir, "upgma_smart.nwk")
upgma_pdf = os.path.join(output_dir, "upgma_smart_tree.pdf")
upgma_png = os.path.join(output_dir, "upgma_smart_tree.png")

if os.path.exists(upgma_nwk):
    print(f"Visualizing UPGMA tree...")
    try:
        visualize_tree(
            upgma_nwk,
            upgma_pdf,
            show_bootstrap=False,  # No bootstrap for UPGMA
            show_branch_length=True,
            title="UPGMA Tree (87 bacterial species)",
            dpi=150,
            width=1400,
            height=1400
        )
        visualize_tree(
            upgma_nwk,
            upgma_png,
            show_bootstrap=False,
            show_branch_length=True,
            title="UPGMA Tree (87 bacterial species)",
            dpi=150,
            width=1400,
            height=1400
        )
        print(f"  [OK] UPGMA visualizations created")
    except Exception as e:
        print(f"  [ERROR] UPGMA visualization failed: {e}")

# ML tree
ml_nwk = os.path.join(output_dir, "ml_smart.nwk")
ml_pdf = os.path.join(output_dir, "ml_smart_tree.pdf")
ml_png = os.path.join(output_dir, "ml_smart_tree.png")

if os.path.exists(ml_nwk):
    print(f"Visualizing ML tree...")
    try:
        visualize_tree(
            ml_nwk,
            ml_pdf,
            show_bootstrap=False,  # No bootstrap in this test
            show_branch_length=True,
            title="ML Tree (GTR+G, 87 bacterial species)",
            dpi=150,
            width=1400,
            height=1400
        )
        visualize_tree(
            ml_nwk,
            ml_png,
            show_bootstrap=False,
            show_branch_length=True,
            title="ML Tree (GTR+G, 87 bacterial species)",
            dpi=150,
            width=1400,
            height=1400
        )
        print(f"  [OK] ML visualizations created")
    except Exception as e:
        print(f"  [ERROR] ML visualization failed: {e}")

print("\nAll visualizations complete!")
print(f"\nOutput files:")
print(f"  {upgma_pdf}")
print(f"  {upgma_png}")
print(f"  {ml_pdf}")
print(f"  {ml_png}")

---
name: tree-visualization
description: Comprehensive guide to phylogenetic tree visualization including matplotlib (Python), ggtree (R), and publication-quality figure generation. Covers rectangular, circular, radial layouts, bootstrap annotations, branch coloring, clade highlighting, and interactive visualizations. Based on best practices from Yulab's Tree Data Book.
---

# Phylogenetic Tree Visualization

## Purpose

Provide comprehensive guidance for creating publication-quality phylogenetic tree visualizations in both Python and R, with emphasis on modern visualization techniques using ggtree and matplotlib.

## When to Use

This skill activates when:
- Visualizing phylogenetic trees
- Creating tree comparison figures
- Adding bootstrap support annotations
- Coloring branches or clades
- Creating circular/radial layouts
- Generating publication-quality figures
- Exporting trees to image formats (PNG, JPG, SVG, PDF)
- Creating interactive tree visualizations
- Customizing tree aesthetics

## Reference Resources

**Primary Resource**: [Yulab Tree Data Book](https://yulab-smu.top/treedata-book/chapter4.html)
- Comprehensive guide to ggtree (R package)
- Modern tree visualization techniques
- Integration with phylogenetic data

## Visualization Approaches

### 1. Python with Matplotlib (Current Implementation)

**Advantages**:
- Pure Python (no R dependency)
- Full control over plotting
- Integrates with existing Python phylogenetic pipeline
- Customizable with matplotlib ecosystem

**Use Cases**:
- Quick visualizations during analysis
- Automated figure generation in pipelines
- Integration with Jupyter notebooks
- Batch processing multiple trees

**Implementation Pattern** (from `create_tree_comparison.py`):

```python
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def draw_tree_rectangular(ax, tree, title):
    """
    Draw phylogenetic tree in rectangular layout.

    Key Steps:
    1. Calculate node coordinates (x=depth, y=leaf position)
    2. Draw edges (horizontal + vertical lines)
    3. Add leaf labels
    4. Style with borders and formatting
    """
    # Get coordinates for all nodes
    coords = calculate_coordinates(tree)

    # Draw edges with proper topology
    draw_edges_recursive(tree, ax, coords)

    # Add leaf labels (italic for species names)
    for leaf_name, (x, y) in leaf_coords.items():
        ax.text(x + offset, y, f'  {leaf_name}',
                style='italic', fontsize=9)

    # Format axes
    ax.set_title(title, fontweight='bold')
    ax.axis('off')
```

**Circular Layout**:
```python
def draw_tree_circular(ax, tree, title):
    """
    Circular/radial tree layout - modern, space-efficient.

    Conversion:
    1. Assign angles to leaves (0-360°)
    2. Calculate radius based on branch lengths
    3. Convert polar (r, θ) to Cartesian (x, y)
    4. Draw arcs for angular connections
    """
    # Polar coordinates
    leaf_angles = assign_leaf_angles(tree)

    # Convert to Cartesian
    for node, (r, theta) in polar_coords.items():
        x = r * cos(theta)
        y = r * sin(theta)
        cartesian_coords[node] = (x, y)

    # Draw arcs at fixed radius
    draw_circular_edges(tree, ax)
```

### 2. R with ggtree (Recommended for Publications)

**Advantages**:
- Publication-quality default aesthetics
- Easy annotation and customization
- Integration with ggplot2 ecosystem
- Extensive tree manipulation functions
- Support for complex phylogenetic data

**Installation**:
```r
# Bioconductor installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("treeio")
```

**Basic ggtree Usage**:
```r
library(ggtree)
library(treeio)

# Read tree from Newick file
tree <- read.tree("tree.nwk")

# Basic rectangular tree
ggtree(tree) +
  geom_tiplab(fontface = "italic") +
  theme_tree2()

# Circular layout
ggtree(tree, layout = "circular") +
  geom_tiplab(aes(angle = angle), fontface = "italic")
```

**Advanced ggtree Features**:
```r
# Add bootstrap values
ggtree(tree) +
  geom_tiplab(fontface = "italic") +
  geom_nodelab(aes(label = support),
               hjust = -0.3, size = 3) +
  geom_treescale()

# Color branches by clade
ggtree(tree, aes(color = clade)) +
  geom_tiplab(fontface = "italic") +
  scale_color_manual(values = c("red", "blue", "green"))

# Highlight clades
ggtree(tree) +
  geom_hilight(node = 15, fill = "steelblue", alpha = 0.3) +
  geom_hilight(node = 23, fill = "darkgreen", alpha = 0.3) +
  geom_tiplab(fontface = "italic")
```

## Layout Types

### Rectangular/Phylogram
**When to use**: Standard phylogenetic trees, when branch lengths are important
```python
# Python
draw_tree_rectangular(ax, tree, title)

# R
ggtree(tree, layout = "rectangular")
```

### Circular/Radial
**When to use**: Large trees (>20 taxa), space efficiency, aesthetic appeal
```python
# Python
draw_tree_circular(ax, tree, title)

# R
ggtree(tree, layout = "circular")
```

### Fan/Unrooted
**When to use**: Unrooted trees, showing all relationships equally
```r
# R only (complex in matplotlib)
ggtree(tree, layout = "daylight")  # Fan layout
ggtree(tree, layout = "equal_angle")  # Unrooted
```

### Cladogram
**When to use**: Focus on topology, ignore branch lengths
```r
ggtree(tree, branch.length = "none") +
  geom_tiplab()
```

## Tree Annotations

### Bootstrap Support Values

**Python Approach**:
```python
def draw_tree_with_bootstrap(ax, tree):
    # Draw basic tree
    draw_tree_rectangular(ax, tree, "Tree with Bootstrap")

    # Add bootstrap values at nodes
    for node in tree.get_internal_nodes():
        if hasattr(node, 'support') and node.support is not None:
            x, y = coords[id(node)]
            ax.text(x, y + 0.02, f'{node.support:.0f}',
                   fontsize=8, ha='center',
                   bbox=dict(boxstyle='round', facecolor='wheat'))
```

**ggtree Approach**:
```r
ggtree(tree) +
  geom_tiplab(fontface = "italic") +
  geom_nodepoint(aes(subset = support > 75, size = support),
                 color = "red", alpha = 0.8) +
  geom_text2(aes(subset = !isTip, label = support),
            hjust = -0.3, size = 3)
```

### Branch Coloring by Group

```r
# Define groups
groupInfo <- split(tree$tip.label,
                   c(rep("Group1", 3), rep("Group2", 4)))

# Color branches
tree <- groupOTU(tree, groupInfo)

ggtree(tree, aes(color = group)) +
  geom_tiplab(fontface = "italic") +
  scale_color_manual(values = c("black", "red", "blue"))
```

### Time Scale (for dated trees)

```r
ggtree(tree, mrsd = "2020-01-01") +
  geom_tiplab(fontface = "italic") +
  theme_tree2() +
  scale_x_continuous(breaks = seq(2015, 2020, 1),
                    labels = seq(2015, 2020, 1))
```

## Tree Comparison Visualizations

### Side-by-Side Comparison (Python)

```python
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Tree Method Comparison', fontsize=16, fontweight='bold')

draw_tree_rectangular(axes[0], ml_tree, 'ML Level 4')
draw_tree_rectangular(axes[1], bionj_tree, 'BioNJ')
draw_tree_rectangular(axes[2], upgma_tree, 'UPGMA')

plt.savefig('comparison.jpg', dpi=300, bbox_inches='tight')
```

### Tanglegram (comparing two trees)

```r
library(dendextend)

# Convert to dendrogram objects
dend1 <- as.dendrogram(tree1)
dend2 <- as.dendrogram(tree2)

# Create tanglegram
tanglegram(dend1, dend2,
          highlight_distinct_edges = TRUE,
          common_subtrees_color_branches = TRUE)
```

### Consensus Tree Visualization

```r
# Multiple trees
trees <- read.tree("trees.nwk")  # Multiple trees in one file

# Calculate consensus
cons <- consensus(trees, p = 0.5)  # Majority-rule

# Visualize with support
ggtree(cons) +
  geom_tiplab(fontface = "italic") +
  geom_nodelab(aes(label = support),
               hjust = -0.3, color = "blue")
```

## Exporting Publication-Quality Figures

### Python Export Options

```python
# High-resolution raster
plt.savefig('tree.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('tree.jpg', dpi=300, bbox_inches='tight', facecolor='white')

# Vector format (scalable)
plt.savefig('tree.pdf', bbox_inches='tight')
plt.savefig('tree.svg', bbox_inches='tight')

# For presentations (lower DPI acceptable)
plt.savefig('tree_presentation.png', dpi=150)
```

### R/ggtree Export

```r
# Save as PDF (vector, publication-ready)
ggsave("tree.pdf", width = 8, height = 6)

# High-resolution PNG
ggsave("tree.png", width = 8, height = 6, dpi = 300)

# SVG for editing in Inkscape/Illustrator
ggsave("tree.svg", width = 8, height = 6)

# For specific journals (e.g., Nature: 89mm single column)
ggsave("tree_nature.pdf", width = 89, height = 120, units = "mm")
```

## Interactive Visualizations

### Python with Plotly

```python
import plotly.graph_objects as go

# Convert tree to plotly format
edges_x, edges_y = [], []
for edge in get_all_edges(tree):
    edges_x.extend([edge.start_x, edge.end_x, None])
    edges_y.extend([edge.start_y, edge.end_y, None])

# Create interactive plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=edges_x, y=edges_y,
                         mode='lines',
                         line=dict(color='black', width=2),
                         hoverinfo='none'))

# Add leaf labels
fig.add_trace(go.Scatter(x=leaf_x, y=leaf_y,
                         mode='text',
                         text=leaf_labels,
                         textposition='middle right'))

fig.update_layout(showlegend=False,
                 xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                 yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))

fig.write_html('tree_interactive.html')
```

### R with ggtree + plotly

```r
library(ggtree)
library(plotly)

p <- ggtree(tree) +
  geom_tiplab(fontface = "italic") +
  theme_tree2()

# Convert to interactive
ggplotly(p)
```

## Best Practices

### Aesthetics
1. **Italic species names**: Scientific convention
2. **Font sizes**: 8-12pt for labels, 10-14pt for titles
3. **Line widths**: 1-2pt for branches
4. **Colors**: Colorblind-friendly palettes (viridis, ColorBrewer)
5. **White background**: Standard for publications

### Layout Selection
- **Rectangular**: Most common, easy to read, <50 taxa
- **Circular**: Space-efficient, aesthetic, 20-100 taxa
- **Fan/Unrooted**: Unrooted networks, no clear outgroup
- **Cladogram**: Topology focus only

### Resolution
- **Publications**: 300 DPI minimum
- **Presentations**: 150 DPI acceptable
- **Web**: 96-150 DPI
- **Posters**: 300-600 DPI

### File Formats
- **PDF/SVG**: Vector, infinitely scalable, best for publications
- **PNG**: Raster, good compression, supports transparency
- **JPG**: Raster, smaller files, no transparency
- **TIFF**: Raster, high quality, some journals require

## Integration with rrna-phylo

### Python Pipeline (Current)

```python
# Build trees
from rrna_phylo.core.builder import PhylogeneticTreeBuilder

builder = PhylogeneticTreeBuilder()
ml_tree, logL = builder.build_ml_tree(sequences)
bionj_tree = builder.build_bionj_tree(sequences)

# Visualize comparison
from create_tree_comparison import create_comparison_image

create_comparison_image(
    sequences,
    seq_type_name='DNA',
    output_prefix='my_trees',
    layout='rectangular'
)
```

### Export to R/ggtree

```python
# Export tree to Newick format
with open('tree.nwk', 'w') as f:
    f.write(ml_tree.to_newick() + ';')
```

Then in R:
```r
library(ggtree)
tree <- read.tree("tree.nwk")

# Create publication figure
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(fontface = "italic", size = 4) +
  geom_treescale(x = 0, y = -1) +
  theme_tree2() +
  ggtitle("Maximum Likelihood Phylogeny") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("publication_tree.pdf", p, width = 6, height = 8)
```

## Common Patterns

### Pattern 1: Quick Comparison Figure
```python
# Generate all trees
trees = builder.build_all_trees(sequences)

# Create comparison
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for ax, (name, tree) in zip(axes, trees.items()):
    draw_tree_rectangular(ax, tree, name.upper())
plt.savefig('quick_comparison.png', dpi=150)
```

### Pattern 2: Publication Figure with Bootstrap
```python
# Build tree with bootstrap
from rrna_phylo.utils.bootstrap import bootstrap_tree

tree_with_support = bootstrap_tree(
    sequences,
    tree_builder_func,
    n_replicates=100
)

# Visualize with support values
draw_tree_with_bootstrap(ax, tree_with_support)
plt.savefig('tree_bootstrap.pdf', bbox_inches='tight')
```

### Pattern 3: Multi-Panel Figure for Paper
```python
fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

ax1 = fig.add_subplot(gs[0, :])  # Full width
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])

# Panel A: Main tree with bootstrap
draw_tree_rectangular(ax1, ml_tree_bootstrap,
                     'A) Maximum Likelihood Tree')

# Panel B: BioNJ comparison
draw_tree_rectangular(ax2, bionj_tree, 'B) BioNJ')

# Panel C: UPGMA comparison
draw_tree_rectangular(ax3, upgma_tree, 'C) UPGMA')

plt.savefig('figure1.pdf', bbox_inches='tight')
```

## Troubleshooting

### Overlapping Labels
```python
# Solution 1: Increase figure height
fig, ax = plt.subplots(figsize=(8, len(taxa) * 0.5))

# Solution 2: Rotate labels (circular layout)
ax.text(x, y, label, rotation=angle)

# Solution 3: Use circular layout for many taxa
draw_tree_circular(ax, tree, title)
```

### Branch Length Issues
```r
# Force equal branch lengths (cladogram)
ggtree(tree, branch.length = "none")

# Log-scale for very different lengths
ggtree(tree) + scale_x_log10()
```

### Color Blindness Accessibility
```python
# Use colorblind-friendly palettes
from matplotlib import cm
colors = cm.get_cmap('viridis')(np.linspace(0, 1, n_groups))

# Or ColorBrewer
from brewer2mpl import get_map
colors = get_map('Set2', 'qualitative', n_groups).mpl_colors
```

## Additional Resources

1. **ggtree Book**: https://yulab-smu.top/treedata-book/
2. **ggtree Paper**: Yu et al. (2017) Methods Ecol Evol
3. **Matplotlib Gallery**: https://matplotlib.org/stable/gallery/
4. **ETE Toolkit**: http://etetoolkit.org/ (Python alternative)
5. **FigTree**: http://tree.bio.ed.ac.uk/software/figtree/ (GUI tool)

## Example Workflow

```python
# Complete workflow: Analysis → Visualization → Export

# 1. Build tree with bootstrap
from rrna_phylo import Sequence
from rrna_phylo.core.builder import PhylogeneticTreeBuilder
from rrna_phylo.utils.bootstrap import bootstrap_tree

sequences = [...]  # Your sequences

builder = PhylogeneticTreeBuilder()
ml_tree, logL = builder.build_ml_tree(sequences)

def ml_builder(seqs):
    t, l, m = build_ml_tree_level4(seqs, model='auto')
    return t

tree_bootstrap = bootstrap_tree(sequences, ml_builder, n_replicates=100)

# 2. Create publication figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

draw_tree_rectangular(ax1, tree_bootstrap,
                     'A) Rectangular Layout')
draw_tree_circular(ax2, tree_bootstrap,
                  'B) Circular Layout')

plt.tight_layout()
plt.savefig('Figure1_phylogeny.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Figure1_phylogeny.png', dpi=300, bbox_inches='tight')

# 3. Export for ggtree (optional)
with open('tree_for_ggtree.nwk', 'w') as f:
    f.write(tree_bootstrap.to_newick() + ';')

print("Visualization complete!")
```

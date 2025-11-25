# Intra-Strain vs Multi-Strain Phylogeny

## Overview

Understanding the difference between analyzing multiple rRNA copies within a single genome (intra-strain) versus comparing sequences across different species/strains (multi-strain) is crucial for accurate phylogenetic interpretation.

## Intra-Strain Phylogeny

**Definition**: Phylogenetic analysis of multiple rRNA gene copies from the **same genome**.

### Key Characteristics:

1. **Multiple rRNA Operons**
   - Most bacteria have 3-15 copies of rRNA operons (16S-23S-5S)
   - E. coli K-12: 7 copies
   - Bacillus subtilis: 10 copies
   - Staphylococcus aureus: 5-6 copies

2. **Expected Branch Lengths**
   - Very short: 0.0000 - 0.0020 substitutions per site
   - Near-identical sequences (>99.5% identity)
   - Differences due to:
     - Recent gene conversion events
     - Slight mutational drift
     - Sequencing/assembly errors

3. **Biological Significance**
   - **Concerted evolution**: rRNA copies evolve together through gene conversion
   - **Homogenization**: Copies stay nearly identical within a genome
   - **Functional constraint**: High selection pressure maintains similarity

### Example from Real Data:

```
E. coli K-12 MG1655 copies (7 total):
├── U00096.223771.225312 (dist: 0.0052)
└── Internal (dist: 0.0056)
    ├── Internal (dist: 0.0032)
    │   ├── Internal (dist: 0.0007)
    │   │   ├── U00096.4166659.4168200 (dist: 0.0000)
    │   │   └── U00096.4035531.4037072 (dist: 0.0007)
    │   └── U00096.4208147.4209688 (dist: 0.0000)
    └── Internal (dist: 0.0010)
        ├── Internal (dist: 0.0008)
        │   ├── U00096.3941808.3943349 (dist: 0.0016)
        │   └── U00096.3427221.3428762 (dist: 0.0017)
        └── U00096.2729616.2731157 (dist: 0.0014)
```

**Interpretation**: All 7 E. coli copies cluster tightly (max distance: 0.0052), confirming they're from the same genome.

## Multi-Strain Phylogeny

**Definition**: Phylogenetic analysis comparing sequences from **different species or strains**.

### Key Characteristics:

1. **Evolutionary Divergence**
   - Branch lengths: 0.01 - 0.30+ substitutions per site
   - Reflects millions of years of evolution
   - Used for taxonomic classification

2. **Expected Branch Lengths**
   - Same genus: 0.01 - 0.05 (E. coli vs Salmonella)
   - Same family: 0.05 - 0.15 (Enterobacteriaceae)
   - Same phylum: 0.15 - 0.30+ (Proteobacteria)
   - Different phyla: 0.30+ (Firmicutes vs Proteobacteria)

3. **Biological Significance**
   - **Species boundaries**: 16S rRNA <97% identity = different species
   - **Taxonomic relationships**: Higher-level classification
   - **Evolutionary history**: Speciation and diversification events

### Example from Real Data:

```
Multi-strain tree:
├── Gram-positive (Firmicutes/Bacillota) - dist: 0.0896
│   ├── Staphylococcus aureus (5 copies cluster: 0.0000-0.0012)
│   └── Bacillus subtilis (0.0421 from Staph)
│
└── Gram-negative (Proteobacteria)
    ├── Enterobacteriaceae - dist: 0.0650
    │   ├── Salmonella (7 copies: 0.0000-0.0018)
    │   └── E. coli (7 copies: 0.0000-0.0052)
    │       └── Distance between species: 0.0078
    │
    └── Pseudomonadaceae - dist: 0.0795
        └── P. aeruginosa (4 copies: 0.0000-0.0007)
            └── Distance to Enterobacteriaceae: 0.0686
```

**Interpretation**:
- **Within-species clustering** (intra-strain): 0.0000-0.0052
- **Between closely related species** (E. coli vs Salmonella): 0.0078
- **Between families** (Enterobacteriaceae vs Pseudomonadaceae): 0.0686
- **Between phyla** (Firmicutes vs Proteobacteria): 0.0896

## Practical Guidelines

### When You Have Multiple Copies from Same Genome:

1. **Expect tight clustering**
   - All copies should have very short branches
   - If distances >0.01, check for:
     - Contamination
     - Misannotation
     - Sequencing errors

2. **Use consensus sequence**
   - For downstream analysis, use one representative copy
   - Or create consensus from all copies
   - Reduces computational burden

3. **Validate genome assembly**
   - Excessive variation suggests assembly issues
   - Check for chimeric contigs

### When Comparing Different Strains/Species:

1. **Use one copy per strain**
   - Select the longest/highest quality copy
   - Or use consensus from all copies
   - Avoids overrepresentation bias

2. **Interpret branch lengths**
   - Short branches (0.00-0.02): Same/closely related species
   - Medium branches (0.02-0.10): Different species, same genus/family
   - Long branches (0.10+): Different families or higher taxonomic levels

3. **Bootstrap analysis**
   - Essential for multi-strain comparisons
   - Indicates confidence in branching order
   - Use 100-1000 replicates for publication

## Common Pitfalls

### Pitfall 1: Treating Intra-Strain Variation as Real Diversity
**Problem**: Including all 7 E. coli copies creates false impression of diversity
**Solution**: Use one representative copy per strain

### Pitfall 2: Expecting Identical Copies
**Problem**: Assuming all copies will be 100% identical
**Reality**: Small differences (0.1-0.5%) are normal due to gene conversion

### Pitfall 3: Ignoring Paralog vs Ortholog Distinction
**Problem**: Comparing wrong copies between strains
**Solution**: Ensure you're comparing orthologous copies (same genomic position)

## Implementation in rRNA-Phylo

### Analyzing Both Types:

```bash
# Full dataset with all copies (shows both intra- and multi-strain)
rrna-phylo all_copies.fasta -o full_results/

# One representative per strain (clean multi-strain phylogeny)
rrna-phylo representatives.fasta -o clean_results/ --bootstrap 100
```

### Interpreting Output:

1. **Look for tight clusters** = Same strain/genome
2. **Check branch lengths**:
   - <0.01 = Likely same genome or very recent divergence
   - 0.01-0.10 = Different species, related genus
   - >0.10 = Different families or higher taxonomy

3. **Use ASCII trees for quick inspection**:
   ```
   # Tight clustering = intra-strain
   |-- Copy1 (dist: 0.0001)
   |-- Copy2 (dist: 0.0002)
   `-- Copy3 (dist: 0.0003)

   # Separated branches = multi-strain
   |-- SpeciesA (dist: 0.0456)
   `-- SpeciesB (dist: 0.0523)
   ```

## References

1. **Gene Conversion in rRNA**: Liao (1999) "Gene conversion drives within genic sequences"
2. **Concerted Evolution**: Dover (1982) "Molecular drive: a cohesive mode of species evolution"
3. **16S rRNA Taxonomy**: Yarza et al. (2014) "Uniting the classification of cultured and uncultured bacteria"

## Key Takeaways

✅ **Intra-strain**: Multiple copies from same genome cluster tightly (dist <0.01)
✅ **Multi-strain**: Different species separate clearly (dist >0.01)
✅ **Best practice**: Use one representative per strain for clean phylogenies
✅ **Validation**: Check that genome copies cluster together before comparing strains

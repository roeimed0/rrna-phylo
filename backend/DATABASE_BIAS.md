# Database Bias in Phylogenetic Analysis

## The Problem

Public sequence databases (GenBank, RefSeq, SILVA) have **massive sampling bias**:

### Examples of Bias:

1. **Species-Level Bias**
   - **Humans**: ~10 billion sequences
   - **E. coli**: ~1 million complete genomes
   - **Mouse**: ~500,000 sequences
   - **Rare tropical species**: 1-10 sequences

2. **Geographic Bias**
   - **USA/Europe**: 85% of sequences
   - **Africa/South America**: <5% of sequences
   - **Deep ocean/cave bacteria**: Extremely underrepresented

3. **Temporal Bias**
   - **2020-2024**: 70% of all sequences (COVID-era sequencing boom)
   - **Pre-2010**: <10% of current database

4. **Methodological Bias**
   - Illumina short reads dominate
   - PacBio/Nanopore long reads underrepresented
   - Culturable bacteria overrepresented vs. unculturable

## Why This Matters for Phylogenetics

### Problem 1: Tree Topology Distortion

**Scenario**: You want to study bacterial evolution

**Bad sampling**:
```
Dataset:
- 1000 E. coli sequences
- 500 Salmonella sequences
- 3 rare species

Result: Tree dominated by E. coli/Salmonella
→ Rare species relationships obscured
→ Branch lengths biased toward model organisms
```

**Good sampling**:
```
Dataset (stratified):
- 10 E. coli representatives
- 10 Salmonella representatives
- 3 rare species

Result: Balanced tree
→ All species weighted equally
→ True evolutionary relationships visible
```

### Problem 2: Bootstrap Support Inflation

**With bias**:
- 1000 E. coli sequences → Bootstrap support artificially high
- Tree "confident" but wrong due to overrepresentation

**Without bias**:
- 10 E. coli representatives → Bootstrap reflects true uncertainty

### Problem 3: Evolutionary Rate Estimation

Overrepresented species can:
- Distort molecular clock estimates
- Bias substitution rate calculations
- Inflate or deflate branch lengths

## Solutions

### 1. Outgroup Selection (For Tree Rooting)

**Why outgroups are essential**:
- Root the tree (determine direction of evolution)
- Detect long-branch attraction artifacts
- Validate topology

**How to choose an outgroup**:

```bash
# For Enterobacteriaceae study (E. coli, Salmonella)
# Use: Pseudomonas or Vibrio (same phylum, different family)
rrna-phylo sequences.fasta --outgroup "Pseudomonas*"

# For Firmicutes study (Bacillus, Staphylococcus)
# Use: E. coli or Pseudomonas (different phylum)
rrna-phylo sequences.fasta --outgroup "U00096.*"

# For all bacteria
# Use: Archaea (different domain)
rrna-phylo sequences.fasta --outgroup "Archaea*"
```

**Rules**:
1. **Outgroup must be monophyletic** (form a single clade)
2. **Outgroup should diverge BEFORE the ingroup** (not too distant, not too close)
3. **Use multiple outgroup taxa** if possible (2-3 species)
4. **Ideal distance**: 10-30% sequence divergence from ingroup

### 2. Stratified Sampling (Balance Species)

**Strategy**: Cap overrepresented species, preserve rare species

```bash
# Cap each species at 10 sequences
rrna-phylo sequences.fasta --stratify --max-per-species 10

# Keep at least 1 sequence per species
rrna-phylo sequences.fasta --stratify --max-per-species 10 --min-per-species 1

# Combined with dereplication
rrna-phylo sequences.fasta --dereplicate --stratify --max-per-species 10
```

**When to use**:
- ✅ Comparing multiple species with unequal representation
- ✅ Model organisms (E. coli, human, mouse) in your dataset
- ✅ Database downloads without pre-filtering
- ❌ Single-species population genetics
- ❌ Within-genome diversity studies

### 3. Taxonomic Coverage Sampling

**Strategy**: Maximize phylogenetic diversity

```bash
# Select 100 most diverse species
rrna-phylo sequences.fasta --max-diversity 100 --prioritize-rare

# Balanced selection
rrna-phylo sequences.fasta --max-diversity 100 --balance
```

**Use case**: Large meta-analysis with thousands of species

### 4. Random Subsampling (For Computational Limits)

**Strategy**: Random sampling to reduce dataset size

```bash
# Randomly sample 1000 sequences
rrna-phylo sequences.fasta --random-sample 1000 --seed 42
```

**Warning**: ⚠️ Random sampling does NOT correct bias!
- Use only if dataset is already balanced
- Or combine with stratification

## Recommended Workflows

### Workflow 1: Multi-Species Comparison (RECOMMENDED)

```bash
# Step 1: Dereplicate (remove gene copy duplicates)
# Step 2: Stratify (balance species representation)
# Step 3: Add outgroup (root the tree)
# Step 4: Bootstrap (assess confidence)

rrna-phylo sequences.fasta \
  --dereplicate \
  --stratify --max-per-species 10 \
  --outgroup "Pseudomonas*" \
  --bootstrap 100 \
  -o results/
```

### Workflow 2: Large Database Analysis

```bash
# For SILVA/RDP downloads with 1000s of species
rrna-phylo sequences.fasta \
  --dereplicate \
  --max-diversity 200 \
  --prioritize-rare \
  --outgroup "Archaea*" \
  --bootstrap 50 \
  -o results/
```

### Workflow 3: Publication-Quality Tree

```bash
# Maximum rigor
rrna-phylo sequences.fasta \
  --dereplicate \
  --stratify --max-per-species 5 \
  --outgroup "outgroup_species*" \
  --bootstrap 1000 \
  -o publication/ \
  --prefix final_tree
```

## Detecting Bias in Your Dataset

### Automated Detection

```bash
# Check for bias (no tree building)
rrna-phylo sequences.fasta --check-bias --quiet
```

Output:
```
=== SAMPLING RECOMMENDATIONS ===

Dataset: 10503 sequences from 23 species

⚠️  OVERREPRESENTED SPECIES DETECTED:
  - Escherichia coli (71.2% of dataset, 7482 seqs)
  - Homo sapiens (18.5% of dataset, 1944 seqs)
  - Mus musculus (8.1% of dataset, 851 seqs)

RECOMMENDATION: Use stratified sampling
  Command: --stratify --max-per-species 10

📊 17 rare species detected (each <1% of dataset)

=== RECOMMENDED WORKFLOW ===
1. Dereplicate: --dereplicate
2. Stratify: --stratify --max-per-species 10
3. Add outgroup: --outgroup "pattern"
4. Bootstrap: --bootstrap 100
```

### Manual Inspection

Check sequence headers:
```bash
# Count sequences per species
grep ">" sequences.fasta | cut -d';' -f8 | sort | uniq -c | sort -rn | head -20
```

## Best Practices

### ✅ DO:

1. **Always dereplicate** if you have multiple rRNA copies per genome
2. **Use outgroups** for all multi-species comparisons
3. **Stratify** when model organisms present
4. **Report sampling strategy** in methods section
5. **Check bias** before building trees

### ❌ DON'T:

1. **Don't use all sequences** from biased databases blindly
2. **Don't ignore rare species** - they may be phylogenetically important
3. **Don't use sister species as outgroup** - won't root correctly
4. **Don't random sample** without stratification
5. **Don't mix different marker genes** (16S, 23S, etc.)

## Literature References

### Database Bias:

1. **Schloss (2016)**: "Application of a database-independent approach to assess the quality of operational taxonomic unit picking methods"
   - Shows how database composition affects OTU clustering

2. **Parks et al. (2020)**: "A complete domain-to-species taxonomy for Bacteria and Archaea"
   - GTDB: addresses RefSeq/GenBank bias

3. **Hinchliff et al. (2015)**: "Synthesis of phylogeny and taxonomy into a comprehensive tree of life"
   - Open Tree of Life: taxonomic sampling issues

### Outgroup Selection:

4. **Huelsenbeck et al. (2002)**: "Potential applications and pitfalls of Bayesian inference"
   - Importance of outgroup choice in rooting

5. **Wheeler (1990)**: "Nucleic acid sequence phylogeny and random outgroups"
   - Theory of outgroup selection

### Sampling Strategies:

6. **Townsend (2007)**: "Profiling phylogenetic informativeness"
   - Optimizing taxon sampling for phylogenetics

7. **Heath et al. (2008)**: "Taxon sampling and the accuracy of phylogenetic analyses"
   - How sampling affects tree topology

## Future Improvements

This module provides basic bias detection and stratified sampling. Future versions could include:

1. **Geographic balancing** (requires location metadata)
2. **Phylogenetic diversity maximization** (PD optimization)
3. **Temporal balancing** (account for sequencing date)
4. **Coverage-based sampling** (weight by sequencing depth)
5. **Integration with SILVA/RDP taxonomies** (automatic classification)

## Quick Reference

| Problem | Solution | Command |
|---------|----------|---------|
| Multiple rRNA copies per genome | Dereplication | `--dereplicate` |
| E. coli/human overrepresented | Stratified sampling | `--stratify --max-per-species 10` |
| Need to root tree | Add outgroup | `--outgroup "pattern"` |
| 1000s of species | Diversity sampling | `--max-diversity 200` |
| Unknown bias | Check first | `--check-bias` |
| Everything | Full pipeline | `--dereplicate --stratify --outgroup X --bootstrap 100` |

---

**Remember**: Good phylogenetic analysis starts with good sampling. Garbage in = garbage out!

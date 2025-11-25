# rRNA-Phylo Usage Guide

Complete guide to building phylogenetic trees from rRNA sequences with proper handling of real-world data challenges.

---

## Quick Start

```bash
# Basic usage (automatic alignment)
rrna-phylo sequences.fasta -o results/

# Recommended workflow (publication-quality)
rrna-phylo sequences.fasta \
  --dereplicate \
  --outgroup "Pseudomonas*" \
  --bootstrap 100 \
  -o results/
```

---

## The Three Key Problems with Real Data

### Problem 1: Multiple rRNA Copies per Genome
**Issue**: Most bacteria have 3-15 rRNA operons in their genome
- E. coli K-12: 7 copies
- Salmonella: 7 copies
- Creates **intra-strain clustering** that obscures true evolution

**Solution**: `--dereplicate`
```bash
rrna-phylo sequences.fasta --dereplicate -o results/
```

**What it does**:
- Groups sequences by genome accession (e.g., all `U00096.*` = E. coli K-12)
- Keeps 1 representative per genome (longest sequence)
- 24 sequences → 5 clean representatives

**Before vs After**:
```
Before (24 seqs):          After (5 seqs):
├─ E. coli copy 1          ├─ E. coli (1 rep)
├─ E. coli copy 2          ├─ Salmonella (1 rep)
├─ E. coli copy 3          ├─ Bacillus (1 rep)
├─ E. coli copy 4          ├─ Staph (1 rep)
├─ E. coli copy 5          └─ Pseudomonas (1 rep)
├─ E. coli copy 6
├─ E. coli copy 7
├─ Salmonella copy 1
... (cluttered!)
```

### Problem 2: Database Bias
**Issue**: Real databases have massive overrepresentation
- 10 billion human sequences
- 1 million E. coli genomes
- 3 sequences of rare tropical species

**Solution**: `--stratify`
```bash
rrna-phylo sequences.fasta --stratify --max-per-species 10 -o results/
```

**What it does**:
- Caps each species at N sequences (default: 10)
- Preserves ALL rare species (≥1 sequence)
- Balances tree without losing diversity

**Example**:
```
Dataset: 10,000 sequences
- 7,000 E. coli (70%)
- 2,000 Human (20%)
- 950 Mouse (9.5%)
- 50 other species (0.5%)

After stratification (max=10):
- 10 E. coli representatives
- 10 Human representatives
- 10 Mouse representatives
- 50 rare species (all kept!)

Result: 80 sequences, balanced representation
```

### Problem 3: Tree Rooting (Outgroups)
**Issue**: Unrooted trees can't show direction of evolution

**Solution**: `--outgroup`
```bash
rrna-phylo sequences.fasta --outgroup "Pseudomonas*" -o results/
```

**What it does**:
- Extracts outgroup sequences by pattern matching
- Outgroup = distant species that diverged BEFORE your ingroup
- Allows proper tree rooting

**Example outgroup choices**:
```
Studying Enterobacteriaceae (E. coli, Salmonella)?
  → Use: --outgroup "Pseudomonas*" (different family)

Studying Gram-positives (Bacillus, Staph)?
  → Use: --outgroup "U00096*" (E. coli = Gram-negative)

Studying all Bacteria?
  → Use: --outgroup "Archaea*" (different domain)
```

---

## Command Reference

### Analysis & Recommendations

#### Check for Database Bias
```bash
rrna-phylo sequences.fasta --check-bias
```

**Output**:
```
=== SAMPLING RECOMMENDATIONS ===
Dataset: 24 sequences from 5 species

[!] OVERREPRESENTED SPECIES DETECTED:
  - Escherichia coli (29.2%, 7 seqs)
  - Salmonella (29.2%, 7 seqs)

RECOMMENDATION: Use stratified sampling
  Command: --stratify --max-per-species 10

=== RECOMMENDED WORKFLOW ===
1. Dereplicate: --dereplicate
2. Stratify: --stratify --max-per-species 10
3. Add outgroup: --outgroup "pattern"
4. Bootstrap: --bootstrap 100
```

#### Suggest Appropriate Outgroups
```bash
rrna-phylo sequences.fasta --suggest-outgroup
```

**Output**:
```
Suggested outgroup for enterobacteriaceae study:
  - Pseudomonas (4 sequences)
    Use: --outgroup "AE004091.*"
  - Bacillus (1 sequence)
    Use: --outgroup "AJ276351.*"
```

### Data Preprocessing

#### Dereplication (Remove Multiple Copies)
```bash
# Keep longest sequence per genome
rrna-phylo seqs.fasta --dereplicate -o results/

# Create consensus from all copies (slower)
rrna-phylo seqs.fasta --dereplicate --derep-method consensus -o results/
```

#### Stratified Sampling (Balance Species)
```bash
# Cap at 10 per species, keep ≥1
rrna-phylo seqs.fasta --stratify --max-per-species 10 -o results/

# Cap at 5, keep at least 2 per species
rrna-phylo seqs.fasta --stratify --max-per-species 5 --min-per-species 2 -o results/
```

#### Outgroup Specification
```bash
# Use wildcard pattern
rrna-phylo seqs.fasta --outgroup "Pseudomonas*" -o results/

# Use specific accession
rrna-phylo seqs.fasta --outgroup "AE004091*" -o results/

# Multiple patterns (separate runs, OR logic in patterns)
rrna-phylo seqs.fasta --outgroup "Pseudomonas*" -o results/
```

### Tree Building

#### Method Selection
```bash
# Maximum Likelihood (recommended, slow)
rrna-phylo seqs.fasta --method ml -o results/

# BioNJ (fast, distance-based)
rrna-phylo seqs.fasta --method bionj -o results/

# UPGMA (fast, assumes molecular clock)
rrna-phylo seqs.fasta --method upgma -o results/

# All methods
rrna-phylo seqs.fasta --method all -o results/
```

#### Bootstrap Analysis
```bash
# 100 bootstrap replicates (recommended for publication)
rrna-phylo seqs.fasta --bootstrap 100 -o results/

# 1000 replicates (gold standard, very slow)
rrna-phylo seqs.fasta --bootstrap 1000 -o results/

# Quick test (10 replicates)
rrna-phylo seqs.fasta --bootstrap 10 -o results/
```

### Output Formats

```bash
# ASCII tree (terminal-friendly)
rrna-phylo seqs.fasta --output-format ascii -o results/

# Newick format (for FigTree, iTOL)
rrna-phylo seqs.fasta --output-format newick -o results/

# Both formats (default)
rrna-phylo seqs.fasta --output-format both -o results/
```

### Alignment Control

```bash
# Skip alignment (sequences already aligned)
rrna-phylo aligned.fasta --skip-align -o results/

# Save aligned sequences to custom location
rrna-phylo seqs.fasta --aligned-output my_aligned.fasta -o results/
```

---

## Complete Workflows

### Workflow 1: Quick Analysis (Single Species)
```bash
# No dereplication needed if all sequences are different species
rrna-phylo sequences.fasta --method ml --bootstrap 10 -o quick_results/
```

**Use case**: Testing, exploration, single-species population genetics

### Workflow 2: Multi-Species Comparison (Standard)
```bash
# Dereplicate + outgroup + bootstrap
rrna-phylo sequences.fasta \
  --dereplicate \
  --outgroup "Pseudomonas*" \
  --method ml \
  --bootstrap 100 \
  --output-format both \
  -o standard_results/
```

**Use case**: Most common scenario, comparing bacterial species

### Workflow 3: Large Database Analysis
```bash
# Check bias first
rrna-phylo silva_download.fasta --check-bias

# Then run with stratification + dereplication
rrna-phylo silva_download.fasta \
  --dereplicate \
  --stratify --max-per-species 10 \
  --outgroup "Archaea*" \
  --method ml \
  --bootstrap 50 \
  -o database_results/
```

**Use case**: Analyzing thousands of sequences from SILVA/RDP/GenBank

### Workflow 4: Publication-Quality Tree
```bash
# Maximum rigor: stratify + dereplicate + high bootstrap
rrna-phylo sequences.fasta \
  --dereplicate \
  --stratify --max-per-species 5 \
  --outgroup "known_outgroup*" \
  --method ml \
  --bootstrap 1000 \
  --output-format both \
  --prefix final_tree \
  -o publication/
```

**Use case**: Final analysis for paper submission

### Workflow 5: Comparing Preprocessing Methods
```bash
# Run 1: No preprocessing (see the problem)
rrna-phylo seqs.fasta --method ml -o raw/

# Run 2: Dereplication only
rrna-phylo seqs.fasta --dereplicate --method ml -o dereplicated/

# Run 3: Full pipeline (see the solution)
rrna-phylo seqs.fasta \
  --dereplicate \
  --stratify --max-per-species 10 \
  --outgroup "Pseudomonas*" \
  --method ml \
  --bootstrap 100 \
  -o full_pipeline/
```

**Use case**: Demonstrating importance of proper data handling

---

## Understanding the Output

### ASCII Tree Format
```
ML Tree:
----------------------------------------
`-- Internal (dist: 0.0000)
    |-- Internal (dist: 0.0026)
    |   |-- Internal (dist: 0.0471)
    |   |   |-- Internal (dist: 0.1379)
    |   |   |   |-- AJ276351.1.1517 (dist: 0.0434)      # Bacillus
    |   |   |   `-- CP000253.2242800.2244356 (dist: 0.0441)  # Staph
    |   |   `-- AE004091.722096.723631 (dist: 0.0743)   # Pseudomonas (OUTGROUP)
    |   `-- AE006468.4394688.4396232 (dist: 0.0128)     # Salmonella
    `-- U00096.223771.225312 (dist: 0.0026)             # E. coli
```

**Interpreting distances**:
- `0.0000-0.0020`: Same genome (intra-strain variation)
- `0.0020-0.0100`: Very closely related species (same genus)
- `0.0100-0.1000`: Different genera, same family
- `0.1000+`: Different families or higher taxonomy

### Bootstrap Support Values
```
(Species_A:0.01,(Species_B:0.02,Species_C:0.03)95:0.01)100:0.0;
                                              ^^          ^^^
                                              Support values (out of 100)
```

**Interpretation**:
- `90-100`: Strong support (trustworthy)
- `70-90`: Moderate support (likely correct)
- `<70`: Weak support (uncertain branching)

### Newick Format
```
((E_coli:0.0026,(Salmonella:0.0128,((Bacillus:0.0434,Staph:0.0441):0.1379,Pseudomonas:0.0743):0.0471):0.0026):0.0000);
```

**Use with**:
- **FigTree**: Desktop app, publication-quality figures
- **iTOL**: Web-based, interactive trees
- **Dendroscope**: Network/tree comparison
- **ggtree (R)**: Programmatic visualization

---

## Troubleshooting

### "No sequences matched outgroup pattern"
**Problem**: Your outgroup pattern didn't match any sequences

**Solutions**:
1. Use `--suggest-outgroup` first to see available patterns
2. Check sequence IDs in your FASTA file
3. Try broader pattern: `"Pseudomonas*"` instead of `"Pseudomonas_aeruginosa"`

### "Need at least 3 sequences"
**Problem**: After dereplication/stratification, too few sequences remain

**Solutions**:
1. Check if your input has enough diversity
2. Reduce `--max-per-species` threshold
3. Don't use `--stratify` if you only have a few sequences

### Alignment is slow
**Problem**: MUSCLE alignment takes forever with 1000+ sequences

**Solutions**:
1. Use stratified sampling to reduce dataset first
2. Dereplicate before alignment
3. Pre-align with faster tools (Clustal Omega) and use `--skip-align`

### Tree has very short branches (everything near 0.000)
**Problem**: All sequences are from same genome (intra-strain)

**Solutions**:
1. Use `--dereplicate` to remove duplicates
2. Check that you're comparing different species, not just one

### Tree topology doesn't match expected phylogeny
**Problem**: Tree shows weird relationships

**Solutions**:
1. Check for overrepresentation bias with `--check-bias`
2. Use stratified sampling to balance species
3. Add outgroup with `--outgroup`
4. Increase bootstrap replicates to check confidence

---

## Best Practices

### Always:
✅ Run `--check-bias` first to understand your data
✅ Use `--dereplicate` if you have genomes with multiple rRNA copies
✅ Specify `--outgroup` for multi-species comparisons
✅ Use `--bootstrap 100` minimum for publication

### Often:
🔄 Use `--stratify` when model organisms (E. coli, human) are present
🔄 Save aligned sequences with `--aligned-output` for record-keeping
🔄 Use `--method all` to compare distance vs ML trees

### Never:
❌ Don't ignore bias warnings
❌ Don't use all 7 E. coli rRNA copies without dereplication
❌ Don't skip bootstrap for publication figures
❌ Don't use random sampling when you need stratified sampling
❌ Don't mix different rRNA genes (16S, 23S, etc.)

---

## Literature References

**Dereplication**:
- Parks et al. (2020) *Nature Biotechnology* - GTDB uses one representative per species

**Stratified Sampling**:
- Schloss (2016) *mSphere* - Database composition affects OTU methods
- Used by QIIME2, mothur, SILVA

**Outgroup Selection**:
- Huelsenbeck et al. (2002) *Systematic Biology* - Importance of outgroup choice

**Bootstrap Analysis**:
- Felsenstein (1985) *Evolution* - Bootstrap confidence limits

---

## Example Test Run

Using the included test dataset (`test_real_rrana.fasta` - 24 sequences, 5 bacterial species):

```bash
# Step 1: Check the data
cd backend
python -m rrna_phylo.cli test_real_rrana.fasta --check-bias

# Step 2: Get outgroup recommendation
python -m rrna_phylo.cli test_real_rrana.fasta --suggest-outgroup

# Step 3: Run full pipeline
python -m rrna_phylo.cli test_real_rrana.fasta \
  --dereplicate \
  --outgroup "AE004091*" \
  --method ml \
  --bootstrap 10 \
  --output-format both \
  -o test_full_workflow/
```

**Expected result**:
- 24 sequences → 5 dereplicated (one per species)
- Clean tree with proper phylogenetic relationships
- Pseudomonas as outgroup for rooting

---

## Getting Help

```bash
# Show all options
rrna-phylo --help

# Test the tool
rrna-phylo test_real_rrana.fasta --check-bias

# Report issues
https://github.com/yourusername/rrna-phylo/issues
```

---

**Remember**: Good phylogenetic analysis starts with good data handling. Always check bias, dereplicate when needed, and use outgroups! 🌳

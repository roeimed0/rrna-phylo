# rRNA-Phylo

Production-ready phylogenetic tree building from rRNA sequences with smart preprocessing and publication-quality visualization.

## Project Goals

- **Phylogenetic Analysis**: Build high-quality phylogenetic trees from rRNA sequences using multiple methods (distance-based, Maximum Likelihood)
- **Generative AI for Trees**: Use machine learning to synthesize trees from sequences and ensemble multiple tree methods (Graph Neural Networks, Transformers)
- **Smart Preprocessing**: Handle real-world data challenges (duplicates, database bias, multiple rRNA operons per genome)
- **Publication-Quality Output**: Generate professional visualizations suitable for scientific publications
- **Production-Ready**: Fast, reliable, and well-tested implementation in pure Python

## Current Features [OK]

### Production-Ready (December 2025)

✅ **Phylogenetic Tree Building**
- UPGMA (distance-based, molecular clock)
- BioNJ (variance-weighted neighbor-joining)
- **ML (Maximum Likelihood with GTR+Gamma model) - 14x FASTER!** ⚡
  - Vectorized Numba-accelerated likelihood calculation (171x faster)
  - 93% biologically reasonable results
  - ~2 minutes for 24 sequences (previously 28 minutes)
  - Production-ready with comprehensive testing
- Automatic model selection (AIC/BIC)
- NNI tree search with optimized branch re-optimization

✅ **Bootstrap Support Analysis** (NEW - December 2025)
- Publication-quality statistical support (Felsenstein 1985)
- 10-1000 bootstrap replicates
- Automatic clade frequency counting
- Support values annotated on tree nodes
- Standard interpretation (Strong ≥95%, Moderate 70-94%, Weak 50-69%)
- Integrated with production CLI

✅ **ETE3 Visualization**
- Publication-quality output (PDF, PNG, SVG, EPS)
- 300-600 DPI resolution
- Aligned leaf names for professional appearance
- Bootstrap support values with color coding
- Used by Nature, Science, and Cell publications

✅ **Smart Deduplication**
- Two-tier strategy (exact + similarity at 99.5%)
- Species-aware clustering prevents phylogenetic errors
- 40-60% dataset reduction without information loss
- 28 comprehensive tests

✅ **Display Names**
- Human-readable species names in trees
- Format: "Species name (Accession) #N"
- Example: "Escherichia coli str. K-12 substr. MG1655 (U00096) #1"

✅ **Complete CLI**
- 20+ command-line options
- Full control over preprocessing, tree building, and visualization
- Comprehensive user guide with real-world workflows

### Quick Start

```bash
# Build a tree with bootstrap and visualization
cd backend
python -m rrna_phylo.cli test_real_rrana.fasta \
    --method ml \
    --bootstrap 100 \
    --dereplicate \
    --visualize \
    -o results/
```

**See [USAGE_GUIDE.md](backend/docs/user-guide/usage-guide.md) for complete examples and workflows.**

## Documentation

### Architecture & Design

📚 **[ARCHITECTURE.md](ARCHITECTURE.md)** - Comprehensive architecture documentation
- Module hierarchy and responsibilities
- Data flow diagrams
- ML likelihood implementation evolution (1x → 171x speedup journey)
- Substitution model hierarchy (JC69 → GTR)
- Performance optimization techniques
- Key interfaces and APIs

📊 **[ARCHITECTURE_DIAGRAMS.md](ARCHITECTURE_DIAGRAMS.md)** - Visual diagrams
- Module dependency graphs
- Likelihood calculation flow
- Tree flattening (vectorization key)
- Site pattern compression
- Performance timeline

### Refactoring Reports

✅ **[REFACTORING_COMPLETE.md](REFACTORING_COMPLETE.md)** - December 2024 refactoring summary
- Eliminated 288 lines of duplicate code
- Unified model interfaces (all 5 substitution models)
- 100% backward compatible, zero performance impact

📋 **[CODE_OPTIMIZATION_REVIEW.md](CODE_OPTIMIZATION_REVIEW.md)** - Detailed code analysis
- Identifies optimization opportunities
- Priority rankings (Critical → Optional)
- Effort estimates and risk assessment

## Development Skills

This project uses Claude Code skills for development guidance:

### Core Skills

1. **[rRNA Prediction Patterns](.claude/skills/rrna-prediction-patterns/SKILL.md)**
   - rRNA detection methods (HMM, BLAST, patterns)
   - Quality scoring and validation
   - Length and completeness assessment

2. **[Phylogenetic Methods](.claude/skills/phylogenetic-methods/SKILL.md)**
   - Distance-based methods (UPGMA, Neighbor-Joining)
   - Maximum likelihood (RAxML, IQ-TREE)
   - Bayesian inference (MrBayes)
   - Tree formats and consensus

3. **[Project Architecture Patterns](.claude/skills/project-architecture-patterns/SKILL.md)**
   - FastAPI backend design
   - Service layer organization
   - Testing strategies
   - API conventions

4. **[ML Integration Patterns](.claude/skills/ml-integration-patterns/SKILL.md)**
   - rRNA classification (supervised learning)
   - Multi-tree consensus (ensemble methods)
   - Generative tree synthesis (GNNs/transformers)

## Getting Started

Ready to start implementing? The skills contain all the patterns and best practices you need.

Simply mention your goal, and the relevant skill will activate to guide you!

## License

MIT License - See LICENSE file for details
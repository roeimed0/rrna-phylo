# rRNA-Phylo

rRNA prediction and phylogenetic tree builder for prokaryotes and eukaryotes.

## Project Goals

- **rRNA Detection**: Identify 16S, 18S, 23S, 28S, 5S, and 5.8S rRNA in genomic sequences
- **Phylogenetic Analysis**: Build trees using multiple methods (distance, ML, Bayesian)
- **Multi-Tree Consensus**: Combine trees for forensics-grade reliability
- **ML Integration**: Enhance detection and tree building with machine learning

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
# rRNA Prediction Patterns Skill

## Overview

This skill provides comprehensive guidelines for implementing rRNA (ribosomal RNA) detection and prediction in both prokaryotic and eukaryotic organisms.

## Activation Triggers

The skill automatically activates when:

### Keywords
- rRNA, ribosomal, ribosomal RNA
- Specific rRNA types: 16S, 18S, 23S, 28S, 5S, 5.8S
- Detection methods: HMM model, BLAST search, conserved regions
- Classification: prokaryotic rRNA, eukaryotic rRNA, SSU, LSU

### Intent Patterns
- Phrases like "detect rRNA", "identify ribosomal", "predict 16S"
- Mentions of rRNA analysis, detection, or prediction
- Questions about sequence identification

### File Patterns
- Files in `backend/rrna/**/*.py`
- Files matching `**/rrna*.py`
- Test files: `tests/**/test_rrna*.py`

### Content Patterns (when editing files)
- Classes or functions related to rRNA
- References to HMMER, HMM profiles
- Database references: SILVA, RDP, Rfam

## What This Skill Covers

- **rRNA Types**: All major prokaryotic (16S, 23S, 5S) and eukaryotic (18S, 28S, 5.8S, 5S) rRNA
- **Detection Methods**: HMM profiles, BLAST searches, pattern-based detection
- **Validation**: Length validation, quality scoring, conserved region checking
- **Error Handling**: Ambiguous predictions, failed detections, edge cases
- **Best Practices**: Multi-method validation, quality assessment, proper error handling

## Usage

When working on rRNA-related code, the skill will automatically be suggested. Use it to ensure:
- Proper detection methodology
- Comprehensive quality scoring
- Robust error handling
- Validation of sequence length and completeness
- Appropriate use of specialized databases (SILVA, RDP, Rfam)

## Configuration

Configuration is managed in [skill-rules.json](../skill-rules.json).

**Type**: Domain skill
**Enforcement**: Suggest (advisory, not blocking)
**Priority**: High

## Testing

You can test the skill activation using:

```bash
cd .claude/hooks
echo '{"session_id":"test","prompt":"I need to implement 16S rRNA detection","cwd":"PROJECT_DIR","transcript_path":"","permission_mode":"auto"}' | npx tsx skill-activation-prompt.ts
```

## Related Skills

- [phylogenetic-methods](../phylogenetic-methods/SKILL.md) - For phylogenetic analysis using rRNA sequences

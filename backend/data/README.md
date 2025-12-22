# Data Folder

Place your FASTA files here for easy access from the interactive menu.

## Supported File Formats

- `.fasta` - Standard FASTA format
- `.fa` - Alternative FASTA extension
- `.fna` - FASTA nucleic acid format

## Example

```
data/
├── my_sequences.fasta
├── bacterial_16s.fasta
└── eukaryotic_18s.fasta
```

When you use the interactive menu (Options 1, 2, or 3), files from this folder will be displayed for quick selection.

## File Requirements

- **Format**: FASTA format with sequence headers starting with `>`
- **Type**: DNA or RNA sequences (rRNA recommended)
- **Count**: At least 3 sequences for tree building
- **Length**: Sequences should be aligned or will be aligned automatically

## Obtaining Example Data

You can download rRNA sequences from:

- **SILVA** (https://www.arb-silva.de/) - Comprehensive rRNA database
- **PR2** (https://pr2-database.org/) - Eukaryotic 18S rRNA
- **NCBI** (https://www.ncbi.nlm.nih.gov/) - General sequence database
- **Greengenes** (https://greengenes.secondgenome.com/) - Bacterial 16S rRNA

## Quick Test

Create a test dataset with:

```bash
python rrna_phylo_test.py
```

This creates `test_sequences.fasta` with 50 sequences (1000 bp each) for testing.

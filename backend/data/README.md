# Data Folder

Place your FASTA files here for easy access from the interactive menu.

## Standard Header Format

**For best results, use this standard format:**

```
>ID|Species_name
```

**Example:**
```
>AB571241|Homo_sapiens
>NC_012920|Callithrix_jacchus
>K03432|Mus_musculus
```

**Benefits:**
- **Traceability**: Keep original accession/ID for reference
- **Readability**: Phylogenetic trees show species names instead of cryptic IDs
- **Consistency**: All trees use the same labeling scheme

**To prepare your FASTA files (deduplicate + clean headers):**
```bash
# Interactive menu (recommended)
python app.py
# Select Option 1: Prepare FASTA File

# Or direct command
python prepare_fasta.py your_data.fasta cleaned_data.fasta
```

## Test Datasets

The `test/` folder contains small FASTA files for quick testing:

- **mammals_test.fasta** - 5 mammalian species
- **birds_test.fasta** - 5 bird species
- **fish_test.fasta** - 5 fish species
- **bacteria_test.fasta** - 5 bacterial species (16S rRNA)
- **fungi_test.fasta** - 5 fungal species (18S rRNA)

These files:
- Ship with the project for immediate testing
- Use standard `ID|Species_name` format
- Are protected from cleanup operations
- Allow you to test tree building without downloading data

## Supported File Formats

- `.fasta` - Standard FASTA format
- `.fa` - Alternative FASTA extension
- `.fna` - FASTA nucleic acid format

## Example Directory

```
data/
├── test/
│   ├── mammals_test.fasta
│   ├── birds_test.fasta
│   ├── fish_test.fasta
│   ├── bacteria_test.fasta
│   └── fungi_test.fasta
├── pr2_Mammalia_2025-12-22.fasta
├── bacterial_16s.fasta
└── eukaryotic_18s.fasta
```

When you use the interactive menu (Options 1, 2, or 3), files from this folder will be displayed for quick selection.

## File Requirements

- **Format**: FASTA format with headers starting with `>`
- **Header**: `>ID|Species_name` (recommended) or any format (will be auto-parsed)
- **Type**: DNA or RNA sequences (rRNA recommended)
- **Count**: At least 3 sequences for tree building
- **Length**: Sequences will be aligned automatically with MUSCLE if needed

## Obtaining Example Data

You can download rRNA sequences from:

- **SILVA** (https://www.arb-silva.de/) - Comprehensive rRNA database
- **PR2** (https://pr2-database.org/) - Eukaryotic 18S rRNA
- **NCBI** (https://www.ncbi.nlm.nih.gov/) - General sequence database


## Quick Test

Create a test dataset with:

```bash
python rrna_phylo_test.py
```

This creates `test_sequences.fasta` with 50 sequences (1000 bp each) for testing.

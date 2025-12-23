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

The `test/` folder contains real biological datasets for testing:

### Unaligned (for testing MUSCLE alignment):
- **mammals_test.fasta** - 34 mammalian species (60 KB)
- **Arcosauria_test.fasta** - 111 reptile/bird species (148 KB)
- **cartilaginous_fish_test.fasta** - 28 shark/ray species (50 KB)

### Pre-aligned (for fast tree building):
- **mammals_test_aligned.fasta** - 34 species, 2132 bp aligned (73 KB)
- **Arcosauria_test_aligned.fasta** - 111 species, 1916 bp aligned (237 KB)
- **cartilaginous_fish_test_aligned.fasta** - 28 species, 1827 bp aligned (57 KB)

**Performance comparison:**
- Unaligned: Runs MUSCLE alignment (slow but automatic)
- Aligned: Use `--pre-aligned` flag (5-10x faster!)

**Usage:**
```bash
# Slow: Runs MUSCLE alignment
python rrna_phylo_cli.py data/test/mammals_test.fasta

# Fast: Skips MUSCLE (use pre-aligned file)
python rrna_phylo_cli.py data/test/mammals_test_aligned.fasta --pre-aligned
```

These files:
- Ship with the project for immediate testing
- Use standard `ID|Species_name` format
- Are protected from cleanup operations
- Contain real 18S rRNA sequences from PR2 database
- Perfect for testing all tree building methods

## Supported File Formats

- `.fasta` - Standard FASTA format
- `.fa` - Alternative FASTA extension
- `.fna` - FASTA nucleic acid format

## Example Directory

```
data/
├── test/
│   ├── mammals_test.fasta (34 species)
│   ├── Arcosauria_test.fasta (111 species)
│   └── cartilaginous_fish_test.fasta (28 species)
├── pr2_Mammalia_2025-12-22.fasta
└── your_custom_data.fasta
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

# IEDB MHC Class II — Negatives Processing Pipeline

An R pipeline for extracting, annotating, and cleaning **negative** MHC Class II peptide binders from the IEDB dataset. Retrieves source protein sequences from UniProt, builds an annotated FASTA, and produces a cleaned, column-selected output CSV.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Dependencies](#dependencies)
4. [Input](#input)
5. [Output Files](#output-files)
6. [Pipeline Sections](#pipeline-sections)
7. [Utility Functions](#utility-functions)
8. [Configuration Notes](#configuration-notes)

---

## Overview

This pipeline focuses exclusively on **negative** MHC Class II ligand records from IEDB (human, Class II). It filters, cleans, and deduplicates peptides; retrieves their source protein sequences from the UniProt REST API in bulk; builds an annotated FASTA file; and saves a peptide-to-protein mapping table and a column-selected CSV for downstream use.

---

## How It Works

```
Raw IEDB CSV
     │
     ▼
 Filter: human · Class II · Negative · valid start position
 Deduplicate rows
     │
     ▼
 Extract & validate UniProt accessions from source IRI
     │
     ▼
 Clean peptide sequences
 Remove PTMs, non-standard amino acids, and empty entries
     │
     ▼
 Fetch protein sequences from UniProt REST API (chunked, 100 IDs/call)
     │
     ▼
 Build annotated FASTA
 (protein header + sequence + peptide position comments)
     │
     ▼
 Save outputs:
   protein_with_peptides.fasta   — annotated protein FASTA
   peptide_mapping.csv           — peptide ↔ protein position table
   sel_cols.csv                  — column-selected filtered table
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading, filtering, and writing |
| `httr` | UniProt REST API calls |
| `Biostrings` | Parsing FASTA sequence files |
| `readr` | Reading input CSV for column selection step |

Install with:

```r
install.packages(c("data.table", "httr", "readr"))
BiocManager::install("Biostrings")
```

---

## Input

| File | Description |
|---|---|
| `mhc_ligand_full.csv` | Full IEDB MHC ligand export (first row skipped — metadata header) |

Expected path:
```
/gpfs/home4/skleynhans/examples_captan/mhc_ligand_full.csv
```

---

## Output Files

| File | Description |
|---|---|
| `protein_with_peptides.fasta` | Annotated FASTA — one entry per protein, with peptide position comments |
| `peptide_mapping.csv` | Mapping table: UniProt ID, peptide sequence, start, end positions |
| `sel_cols.csv` | Column-selected filtered table (`Name.6`, `uniprot_id`, `peptide`) |

---

## Pipeline Sections

### Section 1 — Load and filter
Reads the IEDB CSV (skipping the metadata row) and filters to human, Class II, **Negative** records with a valid peptide start position. Deduplicates immediately after loading.

### Section 2 — Extract UniProt IDs
Parses accessions from the `Source Molecule IRI` column by stripping the base URL and any version suffixes. Validates accessions against the canonical UniProt format regex and drops any that don't match.

### Section 3 — Clean peptide info
Renames the peptide column, drops empty peptides, computes start/end positions, and deduplicates by `(uniprot_id, peptide, start, end)`. Also removes peptides with PTM annotations (`+`) and non-standard amino acid characters.

### Section 4 — Fetch protein sequences
Splits UniProt IDs into chunks of 100 and queries the `rest.uniprot.org/uniprotkb/accessions` endpoint for each chunk. Responses are concatenated into a single FASTA string, written to a temp file, and parsed with `Biostrings::readAAStringSet`. UniProt IDs are extracted from the `|ID|` portion of FASTA headers.

### Section 5 — Build protein lookup
Converts the parsed sequences into a named list (`id → sequence`) for fast lookup. Includes a sanity check printing how many `df` UniProt IDs match the retrieved proteins.

### Section 6 — Build annotated FASTA
Iterates over each retrieved protein and writes:
- A `>UNIPROT_ID` FASTA header
- The full protein sequence
- One comment line per peptide (`# peptide=… start=… end=…`)

### Section 7 — Save mapping table
Writes a simple four-column CSV (`uniprot_id`, `peptide`, `start`, `end`) for all filtered peptides.

### Section 8 — Column selection
Selects only the columns needed downstream (`Name.6`, `uniprot_id`, `peptide`) and writes a clean CSV. Uses a reusable helper function that validates column existence before writing.

---

## Utility Functions

### `merge_ranges(starts, ends)`
Merges overlapping integer ranges (e.g. peptide binding regions) into non-overlapping intervals. Useful for masking or removing peptide-covered regions from protein sequences.

```r
merge_ranges(c(10, 15, 30), c(20, 25, 40))
# Returns: data.table with merged intervals
```

### `create_csv_with_selected_columns(input_file, output_file, selected_columns)`
Reads a CSV, validates that all requested columns exist, subsets to those columns, and writes a new CSV. Wraps everything in `tryCatch` and prints a clear error if any column is missing.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| UniProt chunk size | `chunk_size <- 100` | 100 | Reduce if hitting API rate limits |
| API politeness delay | `Sys.sleep(0.2)` | 0.2s | Increase on shared HPC nodes |
| Input path | `fread(...)` | See script | Update to match your environment |
| Output path | `out_file <- ...` | Working directory | Set `setwd()` or use absolute paths |
| Columns to keep | `columns_to_keep` | `Name.6, uniprot_id, peptide` | Edit to match your downstream needs |

# Shuffle-Based Negative Peptide Generation Pipeline

An R pipeline that generates synthetic **negative** peptides by amino acid shuffling, producing a 5:1 negative-to-positive dataset as both a plain TXT file and a FASTA file.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Dependencies](#dependencies)
4. [Input](#input)
5. [Output Files](#output-files)
6. [Pipeline Sections](#pipeline-sections)
7. [Helper Functions](#helper-functions)
8. [Configuration Notes](#configuration-notes)

---

## Overview

This pipeline takes a file of positive MHC Class II peptides (with or without a flanking context column) and generates shuffled synthetic negatives at a 5:1 ratio. Each negative is produced by randomly permuting the amino acids of a source positive, with guarantees that the result differs from the original and is not already present in either the positive set or the generated negatives.

---

## How It Works

```
Positive peptide TXT (with or without context column)
     │
     ▼
 Strip context column if present
 Deduplicate · drop empty peptides
     │
     ▼
 For each positive peptide:
   └─ Shuffle amino acids up to 100 attempts
   └─ Accept only if: differs from original
                      not in positive set
                      not already generated
   └─ Repeat until 5 unique negatives produced
     │
     ▼
 Collect all negatives → data.table
     │
     ├──► negatives_shuffled.txt   (peptide column, no header)
     └──► negatives_shuffled.faa   (FASTA, one entry per peptide)
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast file reading and writing |

Install with:

```r
install.packages("data.table")
```

---

## Input

| File | Description |
|---|---|
| `DQ_positives.txt` | Tab-separated positive peptides, with or without a context column |

The script auto-detects whether a second (context) column is present:
- **2 columns** — treated as `peptide` + `context12`; context is stripped before processing
- **1 column** — treated as peptide-only; used as-is

Expected input path:
```
/gpfs/home4/skleynhans/report/iedb_data/official_files/DQ_positives.txt
```

---

## Output Files

All outputs are written to `/gpfs/home4/skleynhans/report/iedb_data/shuffled/`.

| File | Description |
|---|---|
| `negatives_shuffled.txt` | Plain peptide list, no header, one peptide per line |
| `negatives_shuffled.faa` | FASTA file — headers follow `DQ_0_{i}` format |

---

## Pipeline Sections

### Section 1 — Read positives and strip context
Reads the input TXT file and detects whether it has one or two columns. If two columns are present, the context column is dropped. Deduplicates and removes empty peptides, then reports the final positive count.

### Section 2 — Generate 5:1 negatives
Calls `generate_negatives()` with `n_per_peptide = 5`. For each positive, up to `n × 100` shuffling attempts are made. A shuffled peptide is accepted only if it is unique across both the positive set and all already-generated negatives. A warning is printed for any peptide that could not yield the full 5 unique negatives within the attempt limit.

### Section 3 — Write outputs
Writes the generated negatives as a plain TXT (peptide column only, no header) and as a FASTA file with sequential headers in `DQ_0_{i}` format.

---

## Helper Functions

### `shuffle_peptide(peptide)`
Randomly permutes the characters of a single peptide string. Retries up to 100 times until the result differs from the input. Returns the original if no unique shuffle is found within the limit (rare for peptides > 2 AA).

### `generate_negatives(peptides, n_per_peptide)`
Iterates over all positive peptides and collects `n_per_peptide` unique shuffled negatives per peptide. Uniqueness is checked globally — a shuffled sequence is rejected if it appears in either the positive set or any previously generated negative. Prints a warning for peptides that could not reach the target count.

### `write_peptide_txt(dt, path)`
Writes a single-column peptide data.table to a headerless TXT file using `fwrite`.

### `write_peptide_fasta(dt, prefix, path)`
Writes peptides to a FASTA file. Headers are formatted as `>{prefix}_{i}` where `i` is a sequential index. Uses `rbind` interleaving of headers and sequences for efficient line construction.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Input file | `INPUT_FILE` | See script | Update locus prefix (DQ/DR/DP) as needed |
| Output TXT | `OUTPUT_TXT` | See script | Update to match locus |
| Output FASTA | `OUTPUT_FAA` | See script | Update to match locus |
| Negatives ratio | `n_per_peptide = 5L` | 5 | Change to adjust ratio |
| Shuffle attempts per peptide | `attempts < 100L` | 100 | Increase for very short peptides |
| Max generation trials | `n_per_peptide * 100L` | 500 | Scales with ratio setting |
| FASTA header prefix | `"DQ_0"` | `"DQ_0"` | Change locus/label to match input |

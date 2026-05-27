# Protein-Derived Negative Peptide Generation Pipeline

An R pipeline that generates **biologically grounded negative peptides** by sampling random windows from known DQ source proteins, masking all experimentally confirmed binding regions before sampling.

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

Unlike shuffle-based negatives (which permute amino acids from known positives), this pipeline generates negatives that are **real protein sub-sequences** — they come from the same source proteins as DQ positives but are sampled from regions that have not been observed to bind. This makes them harder, more realistic negatives for training or evaluating binding predictors.

The strategy:
1. Load DQ positive records from the annotated IEDB table
2. For each source protein, record all known binding intervals
3. Sample random 12–25 AA windows from **unmasked** (non-binding) regions only
4. Reject any candidate that matches a known DQ positive

---

## How It Works

```
Annotated IEDB full CSV
     │
     ▼
 Subset to DQ positives (locus == "DQ", label == 1)
     │
     ├──► Known DQ positive peptide set  (for rejection filter)
     ├──► Per-protein binding region map (start/end intervals)
     └──► Protein sequence lookup        (uniprot_id → sequence)
     │
     ▼
 For each source protein:
   └─ Sample random window length: 12–25 AA
   └─ Sample random start position
   └─ Reject if window overlaps any known binding interval
   └─ Reject if candidate matches any known DQ positive
   └─ Reject if candidate already generated
   └─ Accept → collect; repeat until N_PER_PROTEIN reached
     │
     ▼
 Collect all negatives → data.table
     │
     ├──► DQ_derived_negatives.txt   (one peptide per line, no header)
     └──► DQ_derived_negatives.faa   (FASTA with DQ_0_{i} headers)
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading, filtering, and writing |

Install with:

```r
install.packages("data.table")
```

---

## Input

| File | Description |
|---|---|
| `context_iedb_full.csv` | Annotated IEDB table produced by the main IEDB pipeline — must contain `locus`, `label`, `uniprot_id`, `sequence`, `peptide`, `start`, `end` columns |

Expected path:
```
/gpfs/home4/skleynhans/report/iedb_data/official_files/context_iedb_full.csv
```

---

## Output Files

All outputs are written to `/gpfs/home4/skleynhans/report/iedb_data/rm_binding_shuffle/`.

| File | Description |
|---|---|
| `DQ_derived_negatives.txt` | Plain peptide list, one per line, no header |
| `DQ_derived_negatives.faa` | FASTA file — headers follow `DQ_0_{i}` format |

---

## Pipeline Sections

### Section 1 — Load annotated IEDB data
Reads the full annotated IEDB CSV produced by the upstream pipeline. This file must include sequence, position, locus, and label columns.

### Section 2 — Subset to DQ positives
Filters to `locus == "DQ"` and `label == 1`. Reports the number of DQ positive records and the number of unique source proteins available for sampling.

### Section 3 — Build known positive peptide set
Collects all unique DQ positive peptide sequences into a set used as a rejection filter during sampling. Any sampled window matching a known positive is discarded.

### Section 4 — Build per-protein binding region map
Groups DQ positive records by `uniprot_id` and stores a data frame of `(start, end)` intervals per protein. These intervals define the masked zones that candidate windows must not overlap.

### Section 5 — Build protein sequence lookup
Deduplicates the protein sequence table to one row per `uniprot_id` and drops any rows with missing or empty sequences. Reports how many proteins are available.

### Section 6 — Overlap detection
The `overlaps_binding()` helper checks whether a candidate window `[s, e]` overlaps any interval in a protein's binding map. Two intervals overlap if `s <= interval_end` AND `e >= interval_start`.

### Section 7 — Sample negatives from unmasked regions
Iterates over each source protein. For each, repeatedly:
- Samples a random peptide length (12–25 AA)
- Samples a random start position
- Rejects if the window overlaps a binding interval
- Rejects if the candidate matches a known positive or is already generated
- Accepts if all filters pass

Continues until `N_PER_PROTEIN` negatives are collected or `MAX_ATTEMPTS` is reached. Prints a warning for any protein that could not reach the target count.

### Section 8 — Write outputs
Writes the collected negatives as a plain TXT file (one peptide per line) and a FASTA file with sequential `DQ_0_{i}` headers.

---

## Helper Functions

### `overlaps_binding(s, e, intervals)`
Returns `TRUE` if the window `[s, e]` overlaps at least one interval in the binding map data frame. Used inline during the sampling loop to mask known binding regions.

### `write_peptide_fasta(dt, prefix, path)`
Writes a data.table of peptides to a FASTA file. Headers are formatted as `>{prefix}_{i}` using sequential indices. Uses `rbind` interleaving of header and sequence vectors for efficient line construction.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Input CSV | `FULL_CSV` | See script | Must be the annotated IEDB full table |
| Output TXT | `OUTPUT_TXT` | See script | Update locus prefix if adapting for DR/DP |
| Output FASTA | `OUTPUT_FAA` | See script | Update locus prefix if adapting for DR/DP |
| Negatives per protein | `N_PER_PROTEIN <- 10L` | 10 | Increase for larger negative sets |
| Max sampling attempts | `MAX_ATTEMPTS <- 500L` | 500 | Increase for short proteins or dense binding maps |
| Random seed | `set.seed(42)` | 42 | Ensures reproducibility |
| FASTA header prefix | `"DQ_0"` | `"DQ_0"` | Change to `"DR_0"` or `"DP_0"` to match locus |
| Peptide length range | `sample(12L:25L, 1L)` | 12–25 AA | Adjust to match your predictor's requirements |

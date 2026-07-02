# Protein-Derived Negative Peptide Generation Pipeline

An R pipeline that generates **biologically grounded negative peptides** by masking known DQ binding regions (plus flanking anchor context) from source proteins, shuffling the remaining residues, and sampling random windows from the shuffled pool.

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

This pipeline generates negatives that are **chemically realistic but positionally scrambled**. For each DQ source protein, all experimentally confirmed binding regions are masked — extended by `FLANK_BUFFER` residues on each side to also exclude the immediate anchor-position context flanking each binding core. The residues outside those extended masked regions are concatenated into a pool and shuffled — preserving the per-protein amino acid composition while destroying any positional binding signal. Random 12–25 AA windows are then sampled from this shuffled pool.

This approach minimises false positive risk on two levels: binding regions and their flanking anchor context are excluded entirely, and the shuffled order makes it extremely unlikely that a sampled window will coincidentally form a real binder.

---

## How It Works

```
Annotated IEDB full CSV
     │
     ▼
 Subset to DQ positives (locus == "DQ", label == 1)
     │
     ├──► Known DQ positive peptide set   (rejection filter)
     ├──► Per-protein binding intervals   (start/end per uniprot_id)
     └──► Protein sequence lookup         (uniprot_id → sequence)
     │
     ▼
 For each source protein:
   └─ Merge overlapping binding intervals
   └─ Extend each interval by FLANK_BUFFER residues on each side
   └─ Extract all residues OUTSIDE extended masked intervals → unmasked pool
   └─ Shuffle the unmasked residue pool (preserves AA composition)
     │
     ▼
 For each shuffled pool:
   └─ Sample random window length: 12–25 AA
   └─ Sample random start position in pool
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
Reads the full annotated IEDB CSV from the upstream pipeline. Must include sequence, position, locus, and label columns.

### Section 2 — Subset to DQ positives
Filters to `locus == "DQ"` and `label == 1`. Reports record count and unique source protein count.

### Section 3 — Build known positive peptide set
Collects all unique DQ positive peptide sequences into a character vector used as a rejection filter during sampling. Any sampled window matching a known positive is discarded.

### Section 4 — Build protein sequence lookup
Deduplicates to one row per `uniprot_id` and drops any rows with missing or empty sequences. Reports how many proteins are available for sampling.

### Section 5 — Build per-protein binding region map
Groups DQ positive records by `uniprot_id` and stores the lists of `start` and `end` positions per protein, used when masking each sequence. If a protein has multiple binding peptides, all of their intervals are collected here and passed together to the masking step.

### Section 6 — Build shuffled unmasked residue pools
For each protein, merges all overlapping binding intervals, then extends each merged interval by `FLANK_BUFFER` residues on both the N- and C-terminal sides (clamped to sequence boundaries). This excludes not only the binding core but also the immediate flanking anchor context that is known to influence HLA-DQ peptide register and groove occupancy. All residues outside the extended masked intervals are concatenated into a single string and shuffled, producing a per-protein pool with realistic amino acid composition but no positional binding signal. Proteins where the extended masking leaves fewer than 12 residues in the pool are dropped.

### Section 7 — Sample negatives from shuffled pools
For each protein's shuffled pool, repeatedly samples random 12–25 AA windows until `N_PER_PROTEIN` unique negatives are collected or `MAX_ATTEMPTS` is reached. Candidates are rejected if they match any known DQ positive or are already generated. A warning is printed for proteins that cannot reach the target count. Note that proteins with dense or closely spaced binding regions may produce smaller pools after flanked masking and are more likely to trigger this warning.

### Section 8 — Write outputs
Writes the collected negatives as a plain TXT file (one peptide per line) and a FASTA file with sequential `DQ_0_{i}` headers.

---

## Helper Functions

### `merge_intervals(starts, ends)`
Sorts intervals by start position and merges any that overlap or are adjacent. Returns a `data.table` with columns `start` and `end`. Handles all binding peptides for a given protein simultaneously, so multiple binders on the same protein are correctly consolidated before masking.

### `get_unmasked_residues(sequence, starts, ends, flank = FLANK_BUFFER)`
Calls `merge_intervals`, then extends each merged interval by `flank` residues on each side (clamped to sequence boundaries) before building the position mask. Concatenates all residues at unmasked positions into a single string. Returns an empty string if the entire protein is masked. The `flank` argument defaults to the global `FLANK_BUFFER` parameter but can be overridden per call.

### `write_peptide_fasta(dt, prefix, path)`
Writes a data.table of peptides to a FASTA file. Headers are formatted as `>{prefix}_{i}` using sequential indices.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Input CSV | `FULL_CSV` | See script | Must be the annotated IEDB full table |
| Output TXT | `OUTPUT_TXT` | See script | Update locus prefix if adapting for DR/DP |
| Output FASTA | `OUTPUT_FAA` | See script | Update locus prefix if adapting for DR/DP |
| Negatives per protein | `N_PER_PROTEIN <- 10L` | 10 | Increase for larger negative sets |
| Max sampling attempts | `MAX_ATTEMPTS <- 500L` | 500 | Increase for short proteins or dense binding maps |
| Flank buffer | `FLANK_BUFFER <- 2L` | 2 | Residues masked on each side of every binding interval; increase for more conservative exclusion of anchor context |
| Random seed | `set.seed(42)` | 42 | Ensures reproducibility |
| FASTA header prefix | `"DQ_0"` | `"DQ_0"` | Change to `"DR_0"` or `"DP_0"` to match locus |
| Peptide length range | `sample(12L:25L, 1L)` | 12–25 AA | Adjust to match your predictor's requirements |

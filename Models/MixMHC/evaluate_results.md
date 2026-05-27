# MixMHC2pred Performance Assessment — DQ / DR / DP

An R script that evaluates MixMHC2pred binding predictions per locus, producing ROC-AUC, PR-AUC, and distribution plots for each of DQ, DR, and DP.

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

This script processes one MixMHC2pred prediction file per locus. Labels (positive/negative) are assigned by specifying explicit row ranges in the `LOCI` config block — no separate label file is needed. Blacklisted peptides (negatives confirmed to overlap with known positives) are removed before any metrics are computed. Per-locus plots and a summary table are saved to `OUT_DIR`.

---

## How It Works

```
MixMHC2pred prediction TSV (per locus)
     │
     ▼
 Skip comment lines (#)
 Assign labels from explicit neg_rows / pos_rows ranges
 Validate row numbers are within bounds
     │
     ▼
 Remove blacklisted peptides
 Clean %Rank_best (coerce to numeric, drop NAs)
     │
     ├──► %Rank_best histogram (pos vs neg overlay)
     ├──► ROC curve + AUC
     └──► PR curve + AP
     │
     ▼
 allele_metrics_summary.csv
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading, manipulation, and writing |
| `ggplot2` | All diagnostic plots |
| `pROC` | ROC-AUC computation |
| `PRROC` | Precision-recall AUC computation |

Install with:

```r
install.packages(c("data.table", "ggplot2", "pROC", "PRROC"))
```

---

## Input

One prediction TSV per locus. Each file is the direct output of `MixMHC2pred_unix`. Comment lines starting with `#` are skipped automatically.

| Locus | Variable | Description |
|---|---|---|
| DQ | `pred_file` in `LOCI$DQ` | MixMHC2pred output TSV |
| DR | `pred_file` in `LOCI$DR` | MixMHC2pred output TSV *(currently commented out)* |
| DP | `pred_file` in `LOCI$DP` | MixMHC2pred output TSV *(currently commented out)* |

Labels are **not read from a file** — they are derived from the `neg_rows` and `pos_rows` ranges specified per locus in the `LOCI` config block.

---

## Output Files

All files are written to `OUT_DIR`.

| File | Description |
|---|---|
| `{locus}_labelled.csv` | Prediction rows with assigned Label column |
| `{locus}_Rank_Dist.png` | `%Rank_best` histogram — positive vs negative overlay |
| `{locus}_ROC.png` | ROC curve with AUC annotation |
| `{locus}_PR.png` | Precision-recall curve with AP annotation |
| `allele_metrics_summary.csv` | AUC, AP, N_pos, N_neg per locus |

---

## Pipeline Sections

### CONFIG — Paths, locus definitions, and blacklist
Three configuration blocks are set at the top of the script:

**`OUT_DIR`** — output directory, created automatically if it does not exist.

**`LOCI`** — one entry per locus with:
- `pred_file` — path to the MixMHC2pred output TSV
- `neg_rows` — 1-based row indices for negative peptides (ranges, vectors, or mixed)
- `pos_rows` — 1-based row indices for positive peptides

DR and DP are currently commented out. Uncomment and fill in paths and row ranges to activate them.

**`BLACKLIST`** — per-locus character vectors of peptide sequences to exclude before metrics are computed. These are negatives confirmed to overlap with known positives in the UniProt-derived dataset. Add sequences as new overlaps are identified.

### Per-locus loop

#### Step 1 — Load predictions
Calls `load_pred()` to read the TSV, skipping `#` comment lines. Checks that the file exists before proceeding; skips the locus with a warning if not found.

#### Step 2 — Assign labels
Assigns `Label = 0` to `neg_rows` and `Label = 1` to `pos_rows`. Validates that all supplied row numbers are within the file's bounds and stops with a clear error message if not. Rows not covered by either range are flagged and dropped during cleaning.

#### Step 3 — Remove blacklisted peptides
Removes any peptide in `BLACKLIST[[locus]]` from the prediction table. Reports how many were removed, which sequences were found, and which blacklisted sequences were not present in the file.

#### Step 4 — Clean `%Rank_best`
Coerces `%Rank_best` to numeric (silencing non-numeric warnings) and drops rows where it is `NA` or where `Label` is `NA`.

#### Step 5 — `%Rank_best` histogram
Overlaid histogram of `%Rank_best` for positives and negatives. Lower `%Rank_best` = stronger predicted binder. Saved as `{locus}_Rank_Dist.png`.

#### Step 6 — ROC curve
Computes ROC-AUC using `-(%Rank_best)` as the predictor score (negating so that lower rank = higher score). Saves the curve as `{locus}_ROC.png`.

#### Step 7 — PR curve
Computes PR-AUC using `PRROC::pr.curve()` with negated `%Rank_best`. Saves the curve as `{locus}_PR.png`.

#### Step 8 — Summary table
Accumulates AUC, AP, N_pos, and N_neg per locus. After all loci are processed, prints and saves the combined table to `allele_metrics_summary.csv`.

---

## Helper Functions

### `load_pred(path)`
Reads a MixMHC2pred output file, filters out lines starting with `#`, and parses the remainder as a tab-separated table with a header row.

### `save_plot(p, filename, width, height)`
Saves a ggplot object to `OUT_DIR` at 150 dpi. Prints the saved path to console.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Output directory | `OUT_DIR` | See script | Created automatically |
| Locus pred files | `LOCI$*/pred_file` | See script | One path per locus |
| Negative row range | `LOCI$*/neg_rows` | See script | 1-based; accepts ranges, vectors, or mixed |
| Positive row range | `LOCI$*/pos_rows` | See script | 1-based; accepts ranges, vectors, or mixed |
| Blacklisted peptides | `BLACKLIST$*` | See script | Add sequences as new overlaps are found |
| Active loci | `LOCI` block | DQ only | Uncomment DR and DP entries to activate |
| Plot DPI | `dpi = 150` | 150 | Increase to 300 for publication figures |
| `%Rank_best` direction | `-(%Rank_best)` | Negated | Lower rank = stronger binder; negation makes it a score |

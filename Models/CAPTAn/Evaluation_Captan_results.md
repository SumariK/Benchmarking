# CAPTAn DQ Evaluation Pipeline

An R script that evaluates CAPTAn MHC Class II binding predictions against IEDB ground truth labels for the DQ locus, producing ROC-AUC, PR-AUC, score summary statistics, and five diagnostic plots.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Dependencies](#dependencies)
4. [Input](#input)
5. [Output Files](#output-files)
6. [Pipeline Sections](#pipeline-sections)
7. [Configuration Notes](#configuration-notes)

---

## Overview

This script takes the downsampled IEDB ground truth table and a CAPTAn core prediction CSV, aligns ground truth peptides to CAPTAn scoring windows by positional overlap, and evaluates binding prediction performance. It recovers the original FASTA headers by stripping the trailing index from CAPTAn's unique headers — no separate header mapping file is needed.

---

## How It Works

```
summary_CAPTAn_core.csv              context_iedb_downsampled.csv
  (CAPTAn predictions)                 (ground truth labels)
         │                                      │
         ▼                                      ▼
  Filter: DQ proteins + DQ alleles       Filter to DQ locus
  Recover fasta_header                   Extract peptide, start,
  by stripping trailing _{i}             end, label, fasta_header
         │                                      │
         ▼                                      │
  Collapse to best allele score                 │
  per protein window (Start, End)               │
         │                                      │
         └──────────────┬───────────────────────┘
                        ▼
         Match ground truth peptides to CAPTAn windows
         by fasta_header + positional overlap [Start,End]
         Take highest-scoring overlapping window per peptide
         Unmatched peptides → Score = 0
                        │
                        ▼
         ROC-AUC + 95% CI
         PR-AUC
         Score summary table
                        │
                        ▼
         5 plots (saved as PNG + printed):
           Score distribution
           ScoreCore distribution
           ROC curve
           PR curve
           Score vs ScoreCore scatter
                        │
                        ▼
         captan_dq_evaluation_results.csv
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading, manipulation, and writing |
| `ggplot2` | All diagnostic plots |
| `pROC` | ROC-AUC computation and confidence intervals |
| `PRROC` | Precision-recall AUC computation |

Install with:

```r
install.packages(c("data.table", "ggplot2", "pROC", "PRROC"))
```

---

## Input

| Variable | File | Description |
|---|---|---|
| `CAPTAN_CSV` | `DQ/summary_CAPTAn_core.csv` | CAPTAn core prediction summary |
| `GT_CSV` | `context_iedb_downsampled.csv` | Downsampled IEDB annotated table — ground truth labels |

Both files are resolved relative to `BASE_DIR`.

---

## Output Files

All files are written to `BASE_DIR`.

| File | Description |
|---|---|
| `captan_dq_evaluation_results.csv` | Full matched evaluation table |
| `captan_dq_score_dist.png` | Score density by label |
| `captan_dq_scorecore_dist.png` | ScoreCore density by label |
| `captan_dq_roc.png` | ROC curve with AUC and 95% CI |
| `captan_dq_pr_curve.png` | Precision-recall curve with AUC |
| `captan_dq_score_vs_scorecore.png` | Score vs ScoreCore scatter by label |

---

## Pipeline Sections

### Section 1 — Load CAPTAn output
Loads the CAPTAn CSV and classifies each row by allele type (DQ, DP, DR). Filters to rows where the protein name starts with `DQ_` and the allele contains `DQA` or `DQB`, ensuring DQ proteins are only evaluated against DQ alleles.

### Section 2 — Recover original fasta_header
CAPTAn appends a trailing `_{i}` index to guarantee unique FASTA headers. This section strips that index to recover the original `fasta_header`:
```
DQ_1_03_7  →  DQ_1_03
```
No separate header mapping file is required.

### Section 3 — Load ground truth
Loads the downsampled IEDB CSV, filters to `locus == "DQ"`, and retains `peptide`, `start`, `end`, `label`, and `fasta_header`. Prints the recovered CAPTAn headers alongside ground truth headers for a quick visual sanity check.

### Section 4 — Collapse to best allele score per window
For each unique `(fasta_header, Start, End)` combination, takes the maximum `Score`, `ScoreCore`, and `ScoreContext` across all DQ alleles and records which allele achieved the best score. This collapses the many-alleles-per-window structure down to one row per protein window.

### Section 5 — Match ground truth to CAPTAn windows
For each ground truth peptide, finds all CAPTAn windows on the same protein (`fasta_header` match) that positionally overlap the peptide's `[start, end]` interval. Selects the highest-scoring overlapping window. Three match outcomes are possible:

| `match_type` | Meaning |
|---|---|
| `matched` | At least one overlapping CAPTAn window found |
| `no_overlapping_window` | Protein found but no window overlaps the peptide |
| `protein_not_found` | No CAPTAn rows for this protein at all |

Unmatched peptides (`NA` score) are set to `Score = 0`.

### Section 6 — Evaluation
Computes ROC-AUC with 95% CI using `pROC::roc()` and PR-AUC using `PRROC::pr.curve()`. Prints both metrics to console.

### Section 7 — Score summary table
Prints per-label descriptive statistics: N, Min, Median, Mean, Max, and SD of `Score`.

### Section 8 — Score distribution plot
Density plot of CAPTAn `Score` split by label. Saved as `captan_dq_score_dist.png`.

### Section 9 — ScoreCore distribution plot
Density plot of `ScoreCore` split by label. Saved as `captan_dq_scorecore_dist.png`.

### Section 10 — ROC curve
ROC curve with AUC and 95% CI in the title, plotted against a diagonal random classifier baseline. Saved as `captan_dq_roc.png`.

### Section 11 — PR curve
Precision-recall curve with PR-AUC in the title and a dashed horizontal baseline at the positive class prevalence rate. Saved as `captan_dq_pr_curve.png`.

### Section 12 — Score vs ScoreCore scatter
Scatter plot of `Score` vs `ScoreCore` coloured by label, useful for visualising how core and full scores relate and whether they separate labels similarly. Saved as `captan_dq_score_vs_scorecore.png`.

### Section 13 — Save results
Writes the full matched evaluation table to `captan_dq_evaluation_results.csv`.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Base directory | `BASE_DIR` | See script | Update to your working directory |
| CAPTAn predictions | `CAPTAN_CSV` | `DQ/summary_CAPTAn_core.csv` | Relative to `BASE_DIR` |
| Ground truth CSV | `GT_CSV` | `context_iedb_downsampled.csv` | Relative to `BASE_DIR` |
| Locus filter (CAPTAn) | `grepl("^DQ_", Protein)` | `"DQ"` | Update for DR/DP evaluation |
| Allele filter (CAPTAn) | `grepl("DQA\|DQB", Allele)` | DQ alleles | Update to match locus |
| Locus filter (ground truth) | `locus == "DQ"` | `"DQ"` | Update to match locus |
| Plot dimensions | `ggsave(..., width, height)` | 6–7 × 5–6 in | Adjust as needed |
| Plot resolution | `dpi = 300` | 300 | Reduce for faster previews |

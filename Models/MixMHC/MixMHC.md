# MixMHC2pred DQ — SLURM Job Script

A Bash/SLURM script that runs **MixMHC2pred** binding predictions for all available DQ alleles against a shuffled negative peptide input file, without flanking context.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Dependencies](#dependencies)
4. [Input](#input)
5. [Output Files](#output-files)
6. [Pipeline Sections](#pipeline-sections)
7. [SLURM Configuration](#slurm-configuration)
8. [Configuration Notes](#configuration-notes)

---

## Overview

This script automatically discovers all DQ allele PWM definition files in `~/PWMdef/`, builds the allele list dynamically, and passes it to MixMHC2pred. The `--no_context` flag is used because the input file contains peptide sequences without flanking context columns.

---

## How It Works

```
~/PWMdef/DQA1*.txt files
     │
     ▼
 Extract allele names by stripping path and .txt extension
 (e.g. DQA1_01_01.txt → DQA1_01_01)
     │
     ▼
 MixMHC2pred_unix
   --input    shuffle_mix.txt      (shuffled peptides, no context)
   --output   dq_noctx_shlfd_out.txt
   --alleles  all discovered DQ alleles
   --no_context
     │
     ▼
 dq_noctx_shlfd_out.txt
```

---

## Dependencies

| Requirement | Details |
|---|---|
| SLURM | Workload manager — submit with `sbatch` |
| `MixMHC2pred_unix` | Must be present at `$HOME/MixMHC2pred_unix` |
| `~/PWMdef/` | Directory containing allele PWM definition `.txt` files |

To submit the job:

```bash
sbatch mixmhc_dq.sh
```

---

## Input

| File | Description |
|---|---|
| `shuffle_mix.txt` | Shuffled peptide sequences — plain TXT, no context column |
| `~/PWMdef/DQA1*.txt` | PWM definition files — one per DQ allele, auto-discovered |

---

## Output Files

| File | Description |
|---|---|
| `dq_noctx_shlfd_out.txt` | MixMHC2pred predictions for all DQ alleles |
| `mixmhc_dq.out` | SLURM standard output log |
| `mixmhc_dq.err` | SLURM standard error log |

---

## Pipeline Sections

### Section 1 — Allele discovery
Uses `ls` to find all `DQA1*.txt` files in `~/PWMdef/`, strips the path with `basename`, and strips the `.txt` extension with `sed`. The resulting allele names are stored in `$DQ_ALLELES` as a space-separated list passed directly to the `-a` flag.

### Section 2 — MixMHC2pred
Runs `MixMHC2pred_unix` from `$HOME` with:
- `-i` — input peptide file (shuffled, no context)
- `-o` — output predictions file
- `--allelesFolder` — path to PWM definition files
- `-a` — space-separated list of all discovered DQ alleles
- `--no_context` — disables flanking context scoring; required when input has no context columns

---

## SLURM Configuration

| Directive | Value | Notes |
|---|---|---|
| `--job-name` | `mixmhc_dq` | Name visible in `squeue` |
| `--output` | `mixmhc_dq.out` | Written to submission directory |
| `--error` | `mixmhc_dq.err` | Written to submission directory |
| `--time` | `06:00:00` | 6-hour wall time limit |
| `--mem` | `30G` | Memory allocation per node |

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Input peptides | `-i` flag | See script | Update path to your peptide TXT file |
| Output file | `-o` flag | See script | Update path and filename as needed |
| Allele pattern | `DQA1*.txt` | `DQA1*` | Change to `DPA1*` or `DRB1*` for other loci |
| Alleles folder | `--allelesFolder` | `~/PWMdef` | Update if PWM files are stored elsewhere |
| MixMHC2pred binary | `./MixMHC2pred_unix` | `$HOME` | Must be executable; check with `ls -l ~/MixMHC2pred_unix` |
| Context mode | `--no_context` | off | Remove flag if input file includes flanking context columns |
| Wall time | `--time` | 6 hours | Increase for larger inputs or more alleles |

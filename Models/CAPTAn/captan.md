# CAPTAn SLURM Job Script — DQ Locus

A Bash/SLURM script that runs the full **CAPTAn** MHC Class II binding prediction pipeline on a compute cluster, executing all four CAPTAn stages in sequence.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Dependencies](#dependencies)
4. [Input](#input)
5. [Output](#output)
6. [Pipeline Stages](#pipeline-stages)
7. [SLURM Configuration](#slurm-configuration)
8. [Configuration Notes](#configuration-notes)

---

## Overview

This script submits a CAPTAn binding prediction job to a SLURM cluster. It activates the `captan37` conda environment and runs all four CAPTAn stages — core scoring, context scoring, merging, and final scoring — against a FASTA file of input peptides/proteins.

---

## How It Works

```
SLURM job scheduler
     │
     ▼
 Activate conda environment (captan37)
 Set model and I/O paths
     │
     ▼
 Stage 1: captan-core     — core binding predictions
     │
     ▼
 Stage 2: captan-context  — context-aware predictions (threshold 0.2)
     │
     ▼
 Stage 3: captan-merge    — merge core + context results
     │
     ▼
 Stage 4: captan-score    — final scoring and output
     │
     ▼
 Results written to $OUT directory
```

---

## Dependencies

| Requirement | Details |
|---|---|
| SLURM | Workload manager — submit with `sbatch` |
| Miniconda3 | Installed at `~/miniconda3` |
| `captan37` conda env | Must contain `captan-core`, `captan-context`, `captan-merge`, `captan-score` |

To submit the job:

```bash
sbatch captan_dq.sh
```

---

## Input

| Variable | Path | Description |
|---|---|---|
| `IN_FA` | `uniprot_iedb_results/no_ctx_mixed.faa` | FASTA file of input peptides or proteins |

---

## Output

| Variable | Path | Description |
|---|---|---|
| `OUT` | `uniprot_iedb_results/captan_uniprot` | Output directory — all CAPTAn result files written here |

SLURM log files are written to the submission directory:

| File | Description |
|---|---|
| `captan_{job_id}.out` | Standard output log |
| `captan_{job_id}.err` | Standard error log |

---

## Pipeline Stages

### Stage 1 — `captan-core`
Runs the core binding prediction models against the input FASTA. Uses the models in `$MODEL_CORE_DIR`. The `--details` flag writes per-peptide detail files alongside summary output.

### Stage 2 — `captan-context`
Runs context-aware prediction models. Uses the models in `$MODEL_CTX_DIR`. The `--threshold 0.2` flag sets the binding score cutoff for context predictions. The `--details` flag writes per-peptide detail files.

### Stage 3 — `captan-merge`
Merges the core and context prediction results into a unified output. Uses the context model directory for merge parameters.

### Stage 4 — `captan-score`
Produces the final binding scores from the merged predictions. Uses the scoring model in `$MODEL_DIR`.

---

## SLURM Configuration

| Directive | Value | Notes |
|---|---|---|
| `--job-name` | `captan-dq` | Name visible in `squeue` |
| `--output` | `captan_%j.out` | `%j` is replaced by the job ID |
| `--error` | `captan_%j.err` | `%j` is replaced by the job ID |
| `--time` | `72:00:00` | 72-hour wall time limit |
| `--mem` | `30G` | Memory allocation per node |

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Input FASTA | `IN_FA` | See script | Update to point to your input file |
| Output directory | `OUT` | See script | Update to your desired output path |
| Core models | `MODEL_CORE_DIR` | `~/captan/models/models_core` | Update if models are elsewhere |
| Context models | `MODEL_CTX_DIR` | `~/captan/models/models_context` | Update if models are elsewhere |
| Scoring model | `MODEL_DIR` | `~/captan/bin/captan-score` | Update if binary is elsewhere |
| Context threshold | `--threshold 0.2` | 0.2 | Adjust binding score cutoff as needed |
| Conda env | `captan37` | `captan37` | Must match your installed environment name |
| Wall time | `--time` | 72 hours | Increase for larger inputs |
| Memory | `--mem` | 30G | Increase if jobs are killed with OOM errors |

### CUDA / GPU note
The script contains commented-out `module load` commands for CUDA and cuDNN. If your cluster requires GPU modules, uncomment and update those lines with the exact module versions available on your system (use `module spider CUDA` to list available versions).

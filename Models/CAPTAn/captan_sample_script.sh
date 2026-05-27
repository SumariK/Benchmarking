#!/bin/bash
# ============================================================
# CAPTAn SLURM job — DQ locus
# ============================================================
#
# PURPOSE
# -------
# Runs the full CAPTAn MHC Class II binding prediction pipeline
# on a SLURM compute cluster. Executes all four CAPTAn stages:
#   1. captan-core     — core binding predictions
#   2. captan-context  — context-aware predictions
#   3. captan-merge    — merge core + context results
#   4. captan-score    — final scoring and output
#
# USAGE
#   sbatch captan_dq.sh
# ============================================================

#SBATCH --job-name=captan-dq
#SBATCH --output=captan_%j.out   # %j = SLURM job ID
#SBATCH --error=captan_%j.err
#SBATCH --time=72:00:00          # wall time limit
#SBATCH --mem=30G                # memory per node

## ============================================================
## SECTION 1: Shell safety options
## ============================================================
## set -x  prints every command before executing it — useful
##         for debugging; remove in production if too verbose.
## set -e  exits immediately if any command returns non-zero —
##         prevents silent failures mid-pipeline.
## ============================================================

set -x
set -e

## ============================================================
## SECTION 2: Environment setup
## ============================================================
## Activates the captan37 conda environment.
##
## CUDA / GPU note:
## If your cluster requires GPU modules, uncomment the lines
## below and replace with the exact versions available on your
## system. Use 'module spider CUDA' to list options.
## ============================================================

# module --ignore_cache load "2022"
# module spider CUDA
# module load CUDA/11.0.2-GCC-10.2.0      # replace with version from module spider
# module load cuDNN/8.0.4.30-CUDA-11.0.2  # replace with matching cuDNN version

source ~/miniconda3/etc/profile.d/conda.sh
conda activate captan37

## ============================================================
## SECTION 3: Model paths and I/O variables
## ============================================================
## Update IN_FA and OUT to point to your input FASTA and
## desired output directory before submitting.
## ============================================================

## Model directories
MODEL_CORE_DIR="$HOME/captan/models/models_core"
MODEL_CTX_DIR="$HOME/captan/models/models_context"
MODEL_DIR="$HOME/captan/bin/captan-score"

## Input FASTA and output directory
IN_FA="/gpfs/home4/skleynhans/report/iedb_data/uniprot_iedb_results/mixed.faa"
OUT="/gpfs/home4/skleynhans/report/iedb_data/uniprot_iedb_results/captan_uniprot"

## ============================================================
## STAGE 1: Core binding predictions
## ============================================================
## Runs the core CAPTAn binding prediction models against the
## input FASTA. --details writes per-peptide detail files
## alongside the summary output.
## ============================================================

captan-core \
  --models  $MODEL_CORE_DIR \
  --input   $IN_FA \
  --output  $OUT \
  --details
## ============================================================
## STAGE 2: Context-aware predictions
## ============================================================
## Runs context-aware prediction models using the context model
## directory. --threshold 0.2 sets the binding score cutoff;
## adjust this value to tune sensitivity vs specificity.
## ============================================================

captan-context \
  --models    $MODEL_CTX_DIR \
  --input     $IN_FA \
  --output    $OUT \
  --threshold 0.2 \
  --details

 ## ============================================================
## STAGE 3: Merge core + context results
## ============================================================
## Combines the core and context prediction outputs into a
## single unified result set using the context model directory
## for merge parameters.
## ============================================================

captan-merge \
  --models  $MODEL_CTX_DIR \
  --input   $IN_FA \
  --output  $OUT \
  --details

 ## ============================================================
## STAGE 4: Final scoring
## ============================================================
## Produces the final binding scores from the merged predictions.
## Note the flag order differs from earlier stages:
##   -o  output directory
##   -m  scoring model path
##   -i  input FASTA
## ============================================================

captan-score \
  -o $OUT \
  -m $MODEL_DIR \
  -i $IN_FA 
  


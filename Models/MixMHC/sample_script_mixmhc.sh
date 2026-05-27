#!/bin/bash
# ============================================================
# MixMHC2pred DQ — SLURM job script
# ============================================================
#
# PURPOSE
# -------
# Runs MixMHC2pred binding predictions for all DQ alleles
# against a shuffled negative peptide file. Alleles are
# discovered automatically from ~/PWMdef/DQA1*.txt files.
# Uses --no_context since the input has no flanking columns.
#
# USAGE
#   sbatch mixmhc_dq.sh
#
# INPUT
#   shuffle_mix.txt        — shuffled peptides, no context
#   ~/PWMdef/DQA1*.txt     — DQ allele PWM definitions
#
# OUTPUT
#   dq_noctx_shlfd_out.txt — MixMHC2pred predictions
# ============================================================

#SBATCH --job-name=mixmhc_dq
#SBATCH --output=mixmhc_dq.out   # written to submission directory
#SBATCH --error=mixmhc_dq.err
#SBATCH --time=06:00:00           # 6-hour wall time limit
#SBATCH --mem=30G                 # memory per node

cd $HOME

# ============================================================
# SECTION 1: Allele discovery
# ============================================================
# Finds all DQA1*.txt files in ~/PWMdef/, strips the path
# with basename, and strips the .txt extension with sed.
# The result is a space-separated list of allele names
# passed directly to the -a flag of MixMHC2pred.
#
# Example: DQA1_01_01.txt → DQA1_01_01
#
# To run for a different locus, change the glob pattern:
#   DPA1* for DP alleles
#   DRB1* for DR alleles
# ============================================================

DQ_ALLELES=$(ls ~/PWMdef/DQA1*.txt | xargs -n1 basename | sed 's/\.txt//')

# ============================================================
# SECTION 2: Run MixMHC2pred
# ============================================================
# Runs MixMHC2pred_unix from $HOME with all discovered alleles.
#
# Flag reference:
#   -i               input peptide file (no context)
#   -o               output predictions file
#   --allelesFolder  path to PWM definition files
#   -a               space-separated allele list from Section 1
#   --no_context     disables flanking context scoring;
#                    required when input has no context columns
#                    (remove if the input includes context)
# ============================================================

./MixMHC2pred_unix \
  -i /gpfs/home4/skleynhans/report/iedb_data/shuffled/shuffle_mix.txt \
  -o /gpfs/home4/skleynhans/report/iedb_data/shuffled/mixmhc/results_out.txt \
  --allelesFolder ~/PWMdef \
  -a $DQ_ALLELES \
  --no_context

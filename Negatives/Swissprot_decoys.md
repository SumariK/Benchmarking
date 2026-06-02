# SwissProt-Derived Negative Peptide Generation

A Python pipeline that generates decoy negative peptides by randomly sampling windows from SwissProt protein sequences, matched to the length distribution of known positive peptides.

---

## Table of Contents

1. [Overview](#overview)
2. [Important Note on Study Design](#important-note-on-study-design)
3. [How It Works](#how-it-works)
4. [Dependencies](#dependencies)
5. [Setup](#setup)
6. [Input](#input)
7. [Output Files](#output-files)
8. [Pipeline Sections](#pipeline-sections)
9. [Configuration Notes](#configuration-notes)

---

## Overview

This pipeline samples random peptide-length windows from the SwissProt reviewed protein database to serve as decoy negatives for MHC-II binding prediction benchmarking. Negatives are matched to the length distribution of the positive set, filtered to canonical amino acids only, and screened against both the local positive set and all human Class II confirmed binders from the IEDB MHC ligand report before being accepted. The final output is a shuffled, labelled benchmark dataset ready for model evaluation.

---

## Important Note on Study Design

This pipeline operates under a **positive-unlabelled (PU) learning** assumption. Because SwissProt contains human proteins and MHC-II presents endogenous human peptides, some randomly sampled decoys will inevitably be genuine binders that are simply unlabelled. This has two consequences:

- Model performance metrics such as AUC will be **conservatively underestimated** — the model is penalised for correctly predicting unlabelled true binders as positive
- The problem is made **artificially harder** than it would be with experimentally confirmed non-binders

The IEDB screen integrated into the sampling loop reduces this contamination by rejecting any candidate matching a confirmed human Class II binder in the IEDB MHC ligand report. However, this filter is best-effort — binders that are unpublished, understudied, or restricted to alleles not well represented in IEDB will not be caught. The PU learning caveat therefore still applies and should be noted in any reported results.

---

## How It Works

```
uniprot_sprot.fasta     positives.txt      mhc_ligand_full.csv
        │                    │                      │
        ▼                    ▼                      ▼
 Load protein          Load positive         Load IEDB human
 sequences             peptides              Class II positives
 (plain strings)       (set for fast         Clean + filter to
                       lookup)               canonical AA only
        │                    │                      │
        │                    └──────────────────────┤
        │                                           ▼
        │                              Combined rejection set
        │                          (local positives + IEDB binders)
        │                                           │
        └───────────────────────────────────────────┤
                                                    ▼
                                   Build length distribution
                                   (count positives per length)
                                                    │
                                                    ▼
                                For each length → sample target count:
                                  └─ Pick random protein
                                  └─ Pick random start position
                                  └─ Extract window of target length
                                  └─ Reject if non-canonical AA present
                                  └─ Reject if in combined rejection set
                                  └─ Reject if already in negative set
                                  └─ Accept → add to negatives
                                                    │
                                 ┌──────────────────┴──────────────┐
                                 ▼                                 ▼
                          negatives.txt                 benchmark_dataset.txt
                       (raw negative list)          (shuffled positives + negatives
                                                        with tab-separated labels)
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `biopython` | Parsing SwissProt FASTA |
| `pandas` | Loading and filtering IEDB MHC ligand report |

Install with:

```bash
pip install biopython pandas
```

---

## Setup

### 1. Download SwissProt FASTA

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```

### 2. Download IEDB MHC ligand report

Download the full MHC ligand export from [https://www.iedb.org](https://www.iedb.org) and save as `mhc_ligand_full.csv`.

### 3. Update paths in the CONFIG section

Open `swissprot_negatives.py` and set `SWISSPROT_FASTA`, `POSITIVES_FILE`, and `IEDB_CSV` to your local paths.

### 4. Run

```bash
python swissprot_negatives.py
```

---

## Input

| File | Description |
|---|---|
| `uniprot_sprot.fasta` | SwissProt reviewed protein sequences, downloaded from UniProt FTP |
| `positives.txt` | One positive peptide per line, no header |
| `mhc_ligand_full.csv` | Full IEDB MHC ligand export, used to build the IEDB rejection set |

---

## Output Files

| File | Description |
|---|---|
| `negatives.txt` | Raw decoy negative peptides, one per line, no header |
| `negatives.faa` | FASTA file with deduplicated indexed headers (`neg_1`, `neg_2`, ...) |
| `benchmark_dataset.txt` | Shuffled combined dataset with tab-separated labels — `<peptide>\t<label>` where 1 = positive, 0 = negative |

---

## Pipeline Sections

### Section 1 — Load positive peptides
Reads the positive peptide file, one peptide per line. Stored as a set for O(1) rejection lookup during sampling. Blank lines are skipped automatically.

### Section 2 — Build length distribution
Counts how many positives exist at each peptide length. This distribution drives the sampling target — for each length, `count × NEGATIVES_PER_POSITIVE` decoys are generated, so the negative set mirrors the length profile of the positive set.

### Section 3 — Build IEDB confirmed binder rejection set
Loads the full IEDB MHC ligand report and extracts all human Class II positive peptides as an additional rejection filter. The report is filtered to `Class == II`, `Source Organism == Homo sapiens`, and `Qualitative Measurement == Positive`. Peptides containing PTM modifications (indicated by `+`) and non-canonical amino acids are removed, matching the cleaning logic used in the upstream IEDB processing pipeline. The resulting set is merged with the local positives to form a combined rejection set. This filter is best-effort — binders absent from IEDB will not be caught — so the PU learning caveat still applies.

### Section 4 — Load SwissProt sequences
Parses the full SwissProt FASTA using BioPython and stores sequences as plain strings. Only sequences are retained; headers are discarded. The entire database is held in memory — expect ~1–2 GB RAM usage.

### Section 5 — Sample negatives by length
For each peptide length, randomly selects proteins and start positions until the target count is reached or the safety ceiling (`target × 1000` attempts) is hit. Three rejection criteria are applied in order:

- **Non-canonical amino acids** — peptides containing X, U, O, B, Z, or J are discarded. These residues cause undefined behaviour in most MHC-II predictors: unknown token mapping, silent encoding errors, or crashes. Only the 20 standard amino acids (ACDEFGHIKLMNPQRSTVWY) are accepted.
- **Combined rejection set** — any candidate matching a peptide in the merged local positives + IEDB binder set is discarded.
- **Duplicate negatives** — candidates already accepted in the current run are discarded.

A warning is printed for any length where the target cannot be reached.

### Section 6 — Write raw negatives list
Writes all accepted negatives to a plain text file, one peptide per line, no header. This file can be used independently for downstream processing.

### Section 7 — Build shuffled benchmark dataset
Combines positives (label 1) and negatives (label 0), shuffles the combined list with the fixed random seed, and writes a tab-separated file. Shuffling removes any positional ordering bias that could affect model training or evaluation metrics.

---

## Configuration Notes

| Setting | Default | Notes |
|---|---|---|
| `SWISSPROT_FASTA` | See script | Path to gunzipped SwissProt FASTA |
| `POSITIVES_FILE` | See script | One peptide per line, no header |
| `IEDB_CSV` | See script | Full IEDB MHC ligand export, skip=1 header row |
| `NEGATIVES_OUT` | See script | Output path for raw negatives |
| `NEGATIVES_FAA` | See script | Output path for FASTA file |
| `BENCHMARK_OUT` | See script | Output path for labelled benchmark dataset |
| `NEGATIVES_PER_POSITIVE` | 5 | Decoys generated per positive at each length; increase for larger negative sets |
| `RANDOM_SEED` | 42 | Fixed seed for reproducibility across runs |
| `CANONICAL_AA` | 20 standard AA | Modify only if your predictor handles extended alphabets |

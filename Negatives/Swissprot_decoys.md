# SwissProt-Derived Negative Peptide Generation

An R pipeline that generates decoy negative peptides by randomly sampling windows from SwissProt protein sequences, matched to the length distribution of known positive peptides.

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
        │                          stored as hashed environment
        │                          for O(1) lookup
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
                                  └─ Accept → collect
                                                    │
                          ┌─────────────────────────┴──────────────────┐
                          ▼                         ▼                  ▼
                   negatives.txt              negatives.faa    benchmark_dataset.txt
                (raw negative list)       (FASTA, indexed    (shuffled + labelled,
                                          dedup headers)      tab-separated)
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading, filtering, and writing |
| `Biostrings` | Parsing SwissProt FASTA |

Install with:

```r
install.packages("data.table")
BiocManager::install("Biostrings")
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

Open `swissprot_negatives.R` and set `SWISSPROT_FASTA`, `POSITIVES_FILE`, and `IEDB_CSV` to your local paths.

### 4. Run

```bash
Rscript swissprot_negatives.R
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
Reads the positive peptide file, one peptide per line. Whitespace is stripped and blank lines are skipped. Stored as a character vector and deduplicated for use as a rejection lookup.

### Section 2 — Build length distribution
Counts how many positives exist at each peptide length using `table()`. This distribution drives the sampling target — for each length, `count × NEGATIVES_PER_POSITIVE` decoys are generated, so the negative set mirrors the length profile of the positive set.

### Section 3 — Build IEDB confirmed binder rejection set
Loads the full IEDB MHC ligand report using `data.table::fread` (skip = 1) and filters to human Class II positives. PTM-modified peptides (containing `+`) and non-canonical sequences are removed. The cleaned set is merged with the local positives into a combined rejection set, which is then stored as a hashed R environment for O(1) lookup during sampling. This filter is best-effort — binders absent from IEDB will not be caught.

### Section 4 — Load SwissProt sequences
Parses the full SwissProt FASTA using `Biostrings::readAAStringSet` and converts sequences to plain character strings. The `AAStringSet` object is removed from memory immediately after conversion to free RAM. Expect ~1–2 GB memory usage while the FASTA is being parsed.

### Section 5 — Sample negatives by length
For each peptide length, randomly selects proteins and start positions using `sample.int` until the target count is reached or `target × MAX_ATTEMPTS_FACTOR` attempts are exhausted. A second hashed environment tracks already-accepted negatives for deduplication. Three rejection criteria are applied in order:

- **Non-canonical amino acids** — sequences not matching `^[ACDEFGHIKLMNPQRSTVWY]+$` are discarded.
- **Combined rejection set** — any candidate found in the local positives + IEDB binder environment is discarded.
- **Duplicate negatives** — candidates already accepted in the current run are discarded via the negatives environment.

A warning is printed for any length where the target cannot be reached.

### Section 6 — Write raw negatives list
Writes all accepted negatives to a plain text file using `writeLines`, one peptide per line, no header.

### Section 7 — Write FASTA with deduplicated indexed headers
Writes the accepted negatives to a FASTA file. Headers use a sequential base name (`neg_1`, `neg_2`, ...) which is unique by design. The deduplication logic also handles the general case where protein-name-based headers are used — it tracks occurrences of each base name and appends an incrementing suffix (e.g. `Albumin_1`, `Albumin_2`) when a collision is detected. A summary reports how many headers required disambiguation.

### Section 8 — Build and write shuffled benchmark dataset
Combines positives (label 1) and negatives (label 0) into a `data.table`, shuffles rows using `sample(.N)`, and writes a tab-separated file with no column header using `fwrite`. Shuffling removes any positional ordering bias that could affect model training or evaluation metrics.

---

## Configuration Notes

| Setting | Default | Notes |
|---|---|---|
| `SWISSPROT_FASTA` | See script | Path to gunzipped SwissProt FASTA |
| `POSITIVES_FILE` | See script | One peptide per line, no header |
| `IEDB_CSV` | See script | Full IEDB MHC ligand export, skip = 1 header row |
| `NEGATIVES_OUT` | See script | Output path for raw negatives TXT |
| `NEGATIVES_FAA` | See script | Output path for FASTA file |
| `BENCHMARK_OUT` | See script | Output path for labelled benchmark dataset |
| `NEGATIVES_PER_POSITIVE` | 5 | Decoys generated per positive at each length |
| `MAX_ATTEMPTS_FACTOR` | 1000 | Safety ceiling multiplier per length |
| `RANDOM_SEED` | 42 | Fixed seed passed to `set.seed()` for reproducibility |
| `CANONICAL_AA` | `^[ACDEFGHIKLMNPQRSTVWY]+$` | Regex for canonical AA filter; modify only if your predictor handles extended alphabets |

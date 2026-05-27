# IEDB MHC Class II Ligand Pipeline

An R pipeline for processing the IEDB MHC ligand dataset into peptide input files and protein FASTA files, ready for use with **MixMHC2pred** and **CAPTAn**.

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

This pipeline takes a raw IEDB MHC ligand export and transforms it into clean, analysis-ready files. It handles the full journey from raw CSV to tool-ready FASTA: filtering, protein sequence retrieval from UniProt and NCBI, flanking context computation, deduplication, downsampling, and finally generating CAPTAn-compatible FASTA files with guaranteed unique headers.

Two output sets are produced:

- **Full dataset** — all deduplicated peptides and proteins
- **Downsampled dataset** — biologically realistic 1:5 positive-to-negative ratio

---

## How It Works

### The biological problem

MHC Class II molecules (HLA-DR, DQ, DP) present peptides to CD4+ T cells. The IEDB database records which peptides have been experimentally shown to bind (positives) or not bind (negatives). To train or evaluate binding predictors like MixMHC2pred and CAPTAn, you need:

1. Clean peptide sequences with flanking amino acid context
2. The source protein sequences they came from
3. A realistic class balance (far more negatives than positives exist in biology)

This pipeline automates all of that.

### The data flow

```
Raw IEDB CSV
     │
     ▼
 Filter (human, Class II, HLA-DR/DQ/DP, valid peptide lengths 12–25 AA)
     │
     ▼
 Extract protein accessions (UniProt or NCBI)
     │
     ▼
 Fetch protein sequences via UniProt REST API + NCBI Entrez
 (including recovery of misrouted UniProt IDs filed under NCBI IRIs)
     │
     ▼
 Compute 12-AA flanking context per peptide
 (3 upstream + 3 N-terminal + 3 C-terminal + 3 downstream)
     │
     ▼
 Annotate: locus (DR/DQ/DP), label (1=positive, 0=negative), haplotype
     │
     ▼
 Deduplicate within each locus+label group
     │
     ├──► Full dataset outputs (peptide files, FASTA files, per locus)
     │
     ▼
 Downsample positives to 1:5 ratio (floor(n_negatives / 5), seed=42)
     │
     ├──► Downsampled dataset outputs
     │
     ▼
 Generate CAPTAn-ready FASTA with unique headers + mapping CSVs
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `data.table` | Fast data loading and manipulation |
| `httr` | UniProt REST API calls |
| `Biostrings` | Parsing FASTA sequence files |
| `rentrez` | NCBI Entrez protein sequence retrieval |

Install with:

```r
install.packages(c("data.table", "httr", "rentrez"))
BiocManager::install("Biostrings")
```

---

## Input

| File | Description |
|---|---|
| `mhc_ligand_full.csv` | Full IEDB MHC ligand export (first row skipped — it is a metadata header) |

Place the input file at:
```
/gpfs/home4/skleynhans/examples_captan/mhc_ligand_full.csv
```

---

## Output Files

All outputs are written to `/gpfs/home4/skleynhans/report/iedb_data/`.

### Full dataset

| File | Description |
|---|---|
| `context_iedb_full.csv` | Full annotated table with all columns |
| `mixmhc_context_all.txt` | All peptides + 12-AA context (tab-separated, no header) |
| `iedb_proteins_all.faa` | All source proteins (FASTA) |
| `DR/DQ/DP_positives.txt` | Per-locus positive peptides + context |
| `DR/DQ/DP_negatives.txt` | Per-locus negative peptides + context |
| `DR/DQ/DP_positives.faa` | Per-locus positive source proteins (FASTA) |
| `DR/DQ/DP_negatives.faa` | Per-locus negative source proteins (FASTA) |

### Downsampled dataset (1 positive per 5 negatives)

| File | Description |
|---|---|
| `context_iedb_downsampled.csv` | Downsampled annotated table |
| `mixmhc_context_all_downsampled.txt` | Downsampled peptides + context |
| `iedb_proteins_all_downsampled.faa` | Downsampled source proteins (FASTA) |
| `DR/DQ/DP_positives_downsampled.txt` | Per-locus downsampled positives + context |
| `DR/DQ/DP_negatives_downsampled.txt` | Per-locus negatives + context (unchanged) |
| `DR/DQ/DP_positives_downsampled.faa` | Per-locus downsampled positive proteins |
| `DR/DQ/DP_negatives_downsampled.faa` | Per-locus negative proteins (unchanged) |

### CAPTAn-ready files

| File | Description |
|---|---|
| `DR/DQ/DP_1_captan.faa` | Unique-header positive proteins per locus |
| `DR/DQ/DP_0_captan.faa` | Unique-header negative proteins per locus |
| `DR/DQ/DP_header_mapping.csv` | Old vs. new header mapping per locus |
| `all_loci_header_mapping.csv` | Combined mapping across all loci |

CAPTAn FASTA headers follow the format:
```
{locus}_{label}_{haplotype}_{i}
# e.g. DQ_1_03_1, DQ_0_unknown_12
```
The trailing index `i` guarantees uniqueness within each locus+label group.

---

## Pipeline Sections

### Section 1 — Load and filter
Reads the IEDB CSV and keeps only human, MHC Class II, HLA-DR/DQ/DP entries. Peptides with ambiguous characters, modifications (`+`), or lengths outside 12–25 amino acids are removed.

### Section 2 — Protein ID extraction
Extracts UniProt or NCBI accessions from the source molecule IRI field. Validates accession format and retains negatives even without a valid ID (since negative peptides don't always have a well-characterized source protein).

### Section 3 — Sequence retrieval
Fetches protein sequences in batches:
- **UniProt** — via the `rest.uniprot.org` REST API (100 IDs per chunk)
- **NCBI** — via `rentrez::entrez_fetch` (200 IDs per batch)

A known quirk in the IEDB data is handled here: some UniProt accessions are mistakenly stored under NCBI IRIs. These are detected and rerouted to the UniProt fetcher automatically.

### Section 4 — Flanking context computation
For each peptide, extracts 12 characters of context from the source protein:
- 3 AA upstream of the peptide
- First 3 AA of the peptide (N-terminal)
- Last 3 AA of the peptide (C-terminal)
- 3 AA downstream of the peptide

Positions that fall outside the protein are padded with `"-"` (configurable to `"X"`).

### Section 5 — Annotation
Assigns each row a locus (`DR`, `DQ`, `DP`), a binary label (`1` = positive, `0` = negative), a haplotype extracted from the allele name, and a FASTA header string.

### Section 6 — Deduplication
Removes duplicate peptides within each locus+label group, keeping one representative entry per unique peptide sequence.

### Section 7 — Downsampling
For each locus, positives are randomly downsampled to `floor(n_negatives / 5)` using `set.seed(42)` for reproducibility. Negatives are never downsampled.

### Section 8 — Sequence lookup
Builds a combined named vector of all retrieved sequences (UniProt + NCBI) for use when writing FASTA files.

### Section 9 — Save outputs
Writes all full and downsampled peptide text files and protein FASTA files, organized by locus and label.

### Section 10 — CAPTAn-ready FASTA
Deduplicates by unique protein sequence (not accession), assigns guaranteed-unique headers, writes per-locus FASTA files, and saves header mapping CSVs so predictions can be traced back to their source peptides after CAPTAn evaluation.

---

## Configuration Notes

| Setting | Location | Default | Notes |
|---|---|---|---|
| Flanking pad character | `PAD <- "-"` | `"-"` | Change to `"X"` for some MixMHC2pred versions |
| Random seed | `set.seed(42)` | `42` | Controls downsampling reproducibility |
| Downsampling ratio | `floor(n_neg / 5)` | 1:5 | Hardcoded per locus |
| UniProt chunk size | `chunk_size <- 100` | 100 | Reduce if hitting API rate limits |
| NCBI batch size | `batch_size = 200` | 200 | Reduce if getting Entrez errors |
| Input path | `setwd(...)` / `fread(...)` | See script | Update to match your environment |
| Output path | `out_dir <- ...` | See script | Update to match your environment |

## ============================================================
## Shuffle-based negative peptide generation
## ============================================================
##
## PURPOSE
## -------
## Generates synthetic negative peptides by amino acid shuffling,
## producing a 5:1 negative-to-positive dataset as both TXT and
## FASTA output files.
##
## INPUT
##   DQ_positives.txt — Tab-separated positives (with or without
##                      a flanking context column)
##
## OUTPUT
##   negatives_shuffled.txt — Peptide list, no header
##   negatives_shuffled.faa — FASTA with DQ_0_{i} headers
##
## DEPENDENCIES
##   data.table
## ============================================================

library(data.table)

## ── Paths ──────────────────────────────────────────────────
INPUT_FILE <- "/gpfs/home4/skleynhans/report/iedb_data/official_files/DQ_positives.txt"
OUTPUT_TXT <- "/gpfs/home4/skleynhans/report/iedb_data/shuffled/negatives_shuffled.txt"
OUTPUT_FAA <- "/gpfs/home4/skleynhans/report/iedb_data/shuffled/negatives_shuffled.faa"

## ============================================================
## HELPERS
## ============================================================

## ── Write peptide-only TXT (no header) ─────────────────────
write_peptide_txt <- function(dt, path) {
  fwrite(dt[, .(peptide)], path, col.names = FALSE)
  cat("TXT written:", nrow(dt), "peptides ->", path, "\n")
}

## ── Write peptide-only FASTA ────────────────────────────────
## Headers follow the format: >{prefix}_{i}
write_peptide_fasta <- function(dt, prefix, path) {
  headers     <- paste0(">", prefix, "_", seq_len(nrow(dt)))
  fasta_lines <- as.vector(rbind(headers, dt$peptide))
  writeLines(fasta_lines, path)
  cat("FASTA written:", nrow(dt), "entries ->", path, "\n")
}

## ── Shuffle a single peptide until it differs from original ─
## Retries up to 100 times; returns original if limit reached.
shuffle_peptide <- function(peptide) {
  peptide_chars <- strsplit(peptide, "")[[1]]
  shuffled      <- peptide
  attempts      <- 0L

  while (shuffled == peptide && attempts < 100L) {
    shuffled <- paste0(sample(peptide_chars), collapse = "")
    attempts <- attempts + 1L
  }
  return(shuffled)
}

## ── Generate N unique shuffled negatives per positive ───────
## A shuffled sequence is accepted only if it:
##   (1) differs from the original positive
##   (2) is not already in the positive set
##   (3) is not already in the generated negatives
## Warns if the attempt limit is reached before n_per_peptide
## unique negatives are found.
generate_negatives <- function(peptides, n_per_peptide) {
  positives_set <- as.character(peptides)
  negatives     <- character(0)

  for (pep in peptides) {
    generated  <- 0L
    trials     <- 0L
    max_trials <- n_per_peptide * 100L

    while (generated < n_per_peptide && trials < max_trials) {
      shuffled <- shuffle_peptide(pep)

      if (!(shuffled %in% positives_set) && !(shuffled %in% negatives)) {
        negatives <- c(negatives, shuffled)
        generated <- generated + 1L
      }
      trials <- trials + 1L
    }

    if (generated < n_per_peptide) {
      cat(sprintf(
        "Warning: Could only generate %d/%d unique negatives for peptide %s\n",
        generated, n_per_peptide, pep
      ))
    }
  }
  return(negatives)
}

## ============================================================
## SECTION 1: Read positives and strip context column
## ============================================================
## Auto-detects whether the input file has one column (peptide
## only) or two columns (peptide + context12). Strips the
## context column if present, deduplicates, and drops empties.
## ============================================================

pos_raw <- fread(INPUT_FILE, header = FALSE, sep = "\t")

if (ncol(pos_raw) >= 2) {
  setnames(pos_raw, c("peptide", "context12"))
  pos_dt <- unique(pos_raw[, .(peptide)])
  cat("Context column detected and removed\n")
} else {
  setnames(pos_raw, "peptide")
  pos_dt <- unique(pos_raw)
  cat("No context column detected, proceeding as-is\n")
}

pos_dt <- pos_dt[nchar(peptide) > 0]
cat("Positives loaded:", nrow(pos_dt), "\n")

## ============================================================
## SECTION 2: Generate 5:1 negatives
## ============================================================
## Calls generate_negatives() with n_per_peptide = 5.
## For each positive, up to n * 100 shuffle attempts are made.
## Collected negatives are stored in a data.table.
## ============================================================

cat("Generating 5:1 negatives...\n")

negatives_5to1 <- generate_negatives(
  pos_dt$peptide,
  n_per_peptide = 5L
)

cat("5:1 negatives generated:", length(negatives_5to1), "\n")

neg_dt <- data.table(peptide = negatives_5to1)

## ============================================================
## SECTION 3: Write outputs
## ============================================================
## Writes the generated negatives in two formats:
##   TXT  — plain peptide list, no header
##   FASTA — sequential headers: DQ_0_1, DQ_0_2, ...
##
## Update the FASTA prefix ("DQ_0") to match the locus and
## label of the input file (e.g. "DR_0", "DP_0").
## ============================================================

write_peptide_txt(neg_dt, OUTPUT_TXT)
write_peptide_fasta(neg_dt, "DQ_0", OUTPUT_FAA)

cat("\nDONE\n")

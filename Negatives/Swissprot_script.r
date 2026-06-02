## ============================================================
## DQ-Focused Benchmark Dataset Generation
## ============================================================
##
## PURPOSE
## -------
## Builds a balanced benchmark dataset of 1000 DQ positives
## and 5000 SwissProt-derived decoy negatives.
##
## Positives are drawn directly from the annotated IEDB full
## CSV produced by the upstream processing pipeline (locus = DQ,
## label = 1). No separate positives file is required.
##
## Negatives are sampled as random windows from SwissProt
## protein sequences, length-matched to the DQ positive set.
## Candidates are screened against confirmed positives across
## ALL Class II loci (DQ, DP, DR) to prevent any known binder
## from entering the negative set regardless of locus.
##
## NOTE ON STUDY DESIGN
## --------------------
## Despite the cross-locus rejection screen, some sampled
## decoys may be genuine binders absent from IEDB (unpublished,
## understudied alleles, novel epitopes). This is a
## positive-unlabelled (PU) learning setup — model performance
## metrics such as AUC will be conservative underestimates
## of true performance.
##
## INPUT
##   context_iedb_full.csv  — Annotated IEDB table from upstream
##                            pipeline. Must contain: locus,
##                            label, peptide columns.
##   uniprot_sprot.fasta    — SwissProt FASTA (UniProt FTP)
##
## OUTPUT
##   DQ_benchmark_positives.txt  — 1000 DQ positive peptides
##   DQ_benchmark_negatives.txt  — 5000 SwissProt decoy negatives
##   DQ_benchmark_negatives.faa  — FASTA with indexed headers
##   DQ_benchmark_dataset.txt    — Shuffled combined dataset
##                                 Format: <peptide>\t<label>
##                                 label 1 = positive, 0 = negative
##
## DEPENDENCIES
##   data.table, Biostrings
##   Install: install.packages("data.table")
##            BiocManager::install("Biostrings")
##
## USAGE
##   Rscript swissprot_negatives_DQ.R
##   Adjust paths in the CONFIG section below.
## ============================================================

library(data.table)
library(Biostrings)

## ── Config ──────────────────────────────────────────────────
IEDB_FULL_CSV    <- "/gpfs/home4/skleynhans/report/iedb_data/context_iedb_full.csv"
SWISSPROT_FASTA  <- "/gpfs/home4/skleynhans/report/iedb_data/swissprot/uniprot_sprot.fasta"

OUT_DIR          <- "/gpfs/home4/skleynhans/report/iedb_data/swissprot/"

N_POSITIVES      <- 1000L   ## DQ positives to sample
N_NEGATIVES      <- 5000L   ## SwissProt decoy negatives to generate

MAX_ATTEMPTS_FACTOR <- 1000L  ## safety ceiling = N_NEGATIVES * this
RANDOM_SEED         <- 42L

## Canonical amino acids only
CANONICAL_AA <- "^[ACDEFGHIKLMNPQRSTVWY]+$"

set.seed(RANDOM_SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## SECTION 1: Load annotated IEDB data
## ============================================================
## Reads the full annotated IEDB CSV produced by the upstream
## processing pipeline. Must contain locus, label, and peptide
## columns. Reports counts per locus and label as a sanity check.
## ============================================================

cat("Loading annotated IEDB data...\n")
iedb <- fread(IEDB_FULL_CSV)

cat("Records per locus and label:\n")
print(iedb[, .N, by = .(locus, label)][order(locus, label)])

## ============================================================
## SECTION 2: Sample DQ positives
## ============================================================
## Filters to locus == "DQ" and label == 1. Deduplicates by
## peptide sequence. Samples N_POSITIVES randomly. If fewer
## than N_POSITIVES unique DQ positives are available, all
## are used and a warning is printed.
## ============================================================

cat("\nSampling DQ positives...\n")

dq_pos_all <- unique(iedb[locus == "DQ" & label == 1L], by = "peptide")
cat(sprintf("  Unique DQ positives available: %d\n", nrow(dq_pos_all)))

if (nrow(dq_pos_all) < N_POSITIVES) {
  warning(sprintf(
    "Only %d unique DQ positives available — using all.",
    nrow(dq_pos_all)
  ))
  dq_pos_sample <- dq_pos_all
} else {
  dq_pos_sample <- dq_pos_all[sample(.N, N_POSITIVES)]
}

cat(sprintf("  DQ positives sampled: %d\n", nrow(dq_pos_sample)))

## ============================================================
## SECTION 3: Build cross-locus rejection set
## ============================================================
## Collects all confirmed positives across DQ, DP, and DR to
## build a comprehensive rejection filter. This prevents any
## peptide that is a confirmed binder at ANY Class II locus
## from entering the negative set, regardless of which locus
## it was originally annotated against.
##
## The rejection set is stored as a hashed R environment for
## O(1) lookup during sampling.
## ============================================================

cat("\nBuilding cross-locus rejection set (DQ + DP + DR positives)...\n")

all_class2_pos <- unique(
  iedb[locus %in% c("DQ", "DP", "DR") & label == 1L, peptide]
)
cat(sprintf("  DQ positives:            %d\n",
            uniqueN(iedb[locus == "DQ" & label == 1L, peptide])))
cat(sprintf("  DP positives:            %d\n",
            uniqueN(iedb[locus == "DP" & label == 1L, peptide])))
cat(sprintf("  DR positives:            %d\n",
            uniqueN(iedb[locus == "DR" & label == 1L, peptide])))
cat(sprintf("  Combined rejection set:  %d unique peptides\n",
            length(all_class2_pos)))

rejection_env <- new.env(hash = TRUE, size = length(all_class2_pos))
for (pep in all_class2_pos) assign(pep, TRUE, envir = rejection_env)

## ============================================================
## SECTION 4: Build length distribution from DQ positives
## ============================================================
## Counts how many DQ positives exist at each peptide length.
## SwissProt windows are sampled to match this distribution,
## generating (count / total) * N_NEGATIVES decoys per length.
## Fractional targets are rounded to integers that sum exactly
## to N_NEGATIVES.
## ============================================================

cat("\nBuilding length distribution from DQ positives...\n")

pep_lengths   <- nchar(dq_pos_sample$peptide)
length_counts <- table(pep_lengths)

## Compute proportional targets summing exactly to N_NEGATIVES
props         <- as.numeric(length_counts) / sum(length_counts)
raw_targets   <- props * N_NEGATIVES
floor_targets <- floor(raw_targets)
remainder     <- N_NEGATIVES - sum(floor_targets)

## Distribute remainder to lengths with largest fractional parts
frac_parts    <- raw_targets - floor_targets
top_idx       <- order(frac_parts, decreasing = TRUE)[seq_len(remainder)]
floor_targets[top_idx] <- floor_targets[top_idx] + 1L
length_targets <- as.integer(floor_targets)
names(length_targets) <- names(length_counts)

cat("  Length targets for sampling:\n")
for (i in seq_along(length_targets)) {
  cat(sprintf("    Length %3s: %d positives -> %d negatives\n",
              names(length_targets)[i],
              as.integer(length_counts[i]),
              length_targets[i]))
}
cat(sprintf("  Total negatives to generate: %d\n", sum(length_targets)))

## ============================================================
## SECTION 5: Load SwissProt protein sequences
## ============================================================
## Parses the full SwissProt FASTA using Biostrings and stores
## sequences as plain character strings. The AAStringSet object
## is removed immediately after conversion to free memory.
## Expect ~1-2 GB RAM usage during parsing.
## ============================================================

cat("\nLoading SwissProt sequences...\n")
fa       <- readAAStringSet(SWISSPROT_FASTA)
proteins <- as.character(fa)
rm(fa)
cat(sprintf("  Proteins loaded: %s\n", format(length(proteins), big.mark = ",")))

## ============================================================
## SECTION 6: Sample negatives from SwissProt
## ============================================================
## For each peptide length, randomly samples windows from
## random SwissProt proteins until the per-length target is
## reached or the safety ceiling is hit.
##
## Each candidate is rejected if it:
##   (a) contains non-canonical amino acids
##   (b) matches any confirmed positive in the cross-locus
##       rejection set (DQ + DP + DR)
##   (c) is already in the current negative set (deduplication)
##
## A second hashed environment tracks accepted negatives for
## O(1) deduplication lookup.
## ============================================================

cat("\nSampling negatives from SwissProt...\n")

negatives     <- character(0)
negatives_env <- new.env(hash = TRUE)
n_proteins    <- length(proteins)

for (i in seq_along(length_targets)) {
  
  length       <- as.integer(names(length_targets)[i])
  target       <- length_targets[i]
  max_attempts <- target * MAX_ATTEMPTS_FACTOR
  generated    <- 0L
  attempts     <- 0L
  
  while (generated < target && attempts < max_attempts) {
    attempts <- attempts + 1L
    
    protein  <- proteins[[sample.int(n_proteins, 1L)]]
    prot_len <- nchar(protein)
    if (prot_len < length) next
    
    start     <- sample.int(prot_len - length + 1L, 1L)
    candidate <- substr(protein, start, start + length - 1L)
    
    ## Reject non-canonical residues
    if (!grepl(CANONICAL_AA, candidate)) next
    
    ## Reject if confirmed positive at any Class II locus
    if (exists(candidate, envir = rejection_env, inherits = FALSE)) next
    
    ## Reject duplicates within negative set
    if (exists(candidate, envir = negatives_env, inherits = FALSE)) next
    
    negatives <- c(negatives, candidate)
    assign(candidate, TRUE, envir = negatives_env)
    generated <- generated + 1L
  }
  
  if (generated < target) {
    cat(sprintf(
      "  Warning: only generated %d/%d negatives for length %d\n",
      generated, target, length
    ))
  } else {
    cat(sprintf("  Length %3d: %d negatives generated\n", length, generated))
  }
}

cat(sprintf("\nTotal negatives generated: %s\n",
            format(length(negatives), big.mark = ",")))

## ============================================================
## SECTION 7: Write positives and raw negatives lists
## ============================================================
## Writes each set as a plain text file, one peptide per line,
## no header.
## ============================================================

pos_out <- file.path(OUT_DIR, "DQ_benchmark_positives.txt")
neg_out <- file.path(OUT_DIR, "DQ_benchmark_negatives.txt")

writeLines(dq_pos_sample$peptide, pos_out)
cat(sprintf("Positives written: %d peptides -> %s\n",
            nrow(dq_pos_sample), pos_out))

writeLines(negatives, neg_out)
cat(sprintf("Negatives written: %d peptides -> %s\n",
            length(negatives), neg_out))

## ============================================================
## SECTION 8: Write FASTA with deduplicated indexed headers
## ============================================================
## Writes negatives to FASTA format. Headers use a sequential
## base name (neg_1, neg_2, ...) which is unique by design.
## The deduplication counter also handles the general case of
## protein-name-based headers by appending an incrementing
## suffix (e.g. Albumin_1, Albumin_2) on collision.
## ============================================================

cat("\nWriting FASTA with deduplicated headers...\n")

faa_out      <- file.path(OUT_DIR, "DQ_benchmark_negatives.faa")
name_counts  <- list()
fasta_lines  <- character(length(negatives) * 2L)
dupes_found  <- 0L

`%||%` <- function(a, b) if (!is.null(a)) a else b

for (i in seq_along(negatives)) {
  base              <- paste0("neg_", i)
  name_counts[[base]] <- (name_counts[[base]] %||% 0L) + 1L
  n                 <- name_counts[[base]]
  
  if (n > 1L) {
    header     <- paste0(base, "_", n)
    dupes_found <- dupes_found + 1L
  } else {
    header <- base
  }
  
  fasta_lines[2L * i - 1L] <- paste0(">", header)
  fasta_lines[2L * i]      <- negatives[i]
}

writeLines(fasta_lines, faa_out)

if (dupes_found > 0L) {
  cat(sprintf("  Disambiguated %d duplicate header(s)\n", dupes_found))
} else {
  cat("  No duplicate headers detected\n")
}
cat(sprintf("FASTA written: %d entries -> %s\n", length(negatives), faa_out))

## ============================================================
## SECTION 9: Build and write shuffled benchmark dataset
## ============================================================
## Combines DQ positives (label 1) and decoy negatives (label 0)
## into a single data.table, shuffles rows, and writes a
## tab-separated file with no column header.
## Format: <peptide>\t<label>
## ============================================================

cat("\nBuilding shuffled benchmark dataset...\n")

benchmark <- data.table(
  peptide = c(dq_pos_sample$peptide, negatives),
  label   = c(rep(1L, nrow(dq_pos_sample)), rep(0L, length(negatives)))
)
benchmark <- benchmark[sample(.N)]

bench_out <- file.path(OUT_DIR, "DQ_benchmark_dataset.txt")
fwrite(benchmark, bench_out, sep = "\t", col.names = FALSE)

cat(sprintf("Benchmark dataset written -> %s\n", bench_out))
cat(sprintf("  DQ positives: %d\n", nrow(dq_pos_sample)))
cat(sprintf("  Negatives:    %d\n", length(negatives)))
cat(sprintf("  Total:        %d\n", nrow(benchmark)))
cat("\nDONE\n")

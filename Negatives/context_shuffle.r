## ============================================================
## Protein-derived negative peptide generation — DQ locus
## ============================================================
##
## PURPOSE
## -------
## Generates biologically grounded negative peptides by:
##   1. Masking all known DQ binding regions per source protein,
##      extended by FLANK_BUFFER residues on each side to exclude
##      immediate anchor-position context
##   2. Concatenating residues outside masked regions into a pool
##   3. Shuffling that residue pool (preserves AA composition,
##      destroys positional binding signal)
##   4. Sampling random 12-25 AA windows from the shuffled pool
##   5. Validating candidates against known DQ positives
##
## INPUT
##   context_iedb_full.csv — Annotated IEDB table (upstream pipeline)
##                           Must contain: locus, label, uniprot_id,
##                           sequence, peptide, start, end
##
## OUTPUT
##   DQ_derived_negatives.txt — Plain peptide list, no header
##   DQ_derived_negatives.faa — FASTA with DQ_0_{i} headers
##
## DEPENDENCIES
##   data.table
## ============================================================

library(data.table)

## ── Paths ──────────────────────────────────────────────────
FULL_CSV   <- "/gpfs/home4/skleynhans/report/iedb_data/official_files/context_iedb_full.csv"
OUTPUT_TXT <- "/gpfs/home4/skleynhans/report/iedb_data/rm_binding_shuffle/DQ_derived_negatives.txt"
OUTPUT_FAA <- "/gpfs/home4/skleynhans/report/iedb_data/rm_binding_shuffle/DQ_derived_negatives.faa"

## ── Sampling parameters ────────────────────────────────────
N_PER_PROTEIN <- 10L    ## negatives to generate per source protein
MAX_ATTEMPTS  <- 500L   ## max sampling tries per shuffled pool
FLANK_BUFFER  <- 2L     ## residues to mask on each side of a binding region
                        ## captures immediate anchor-position context flanking
                        ## the binding core; increase for more conservative masking
set.seed(42)

## ============================================================
## HELPERS
## ============================================================

## ── Merge overlapping binding intervals ─────────────────────
## Sorts intervals by start and merges any that overlap or are
## adjacent. Returns a data.table with columns start and end.
merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) return(data.table(start = integer(0), end = integer(0)))

  ord    <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]

  cur_s <- starts[1]
  cur_e <- ends[1]

  m_starts <- c()
  m_ends   <- c()

  if (length(starts) >= 2) {
    for (i in seq(2, length(starts))) {
      if (starts[i] <= cur_e) {
        cur_e <- max(cur_e, ends[i])
      } else {
        m_starts <- c(m_starts, cur_s)
        m_ends   <- c(m_ends,   cur_e)
        cur_s    <- starts[i]
        cur_e    <- ends[i]
      }
    }
  }

  m_starts <- c(m_starts, cur_s)
  m_ends   <- c(m_ends,   cur_e)

  data.table(start = m_starts, end = m_ends)
}

## ── Extract unmasked residues from a protein ────────────────
## Merges binding intervals, extends each by FLANK_BUFFER on
## both sides (clamped to sequence boundaries), then concatenates
## all residues outside every masked interval into a single string.
## The buffer excludes immediate anchor-position context flanking
## the binding core, producing cleaner negative candidates.
get_unmasked_residues <- function(sequence, starts, ends, flank = FLANK_BUFFER) {
  prot_len  <- nchar(sequence)
  intervals <- merge_intervals(starts, ends)

  ## Build a logical mask: TRUE = unmasked position
  mask <- rep(TRUE, prot_len)
  for (i in seq_len(nrow(intervals))) {
    s <- max(1L,        intervals$start[i] - flank)
    e <- min(prot_len,  intervals$end[i]   + flank)
    mask[s:e] <- FALSE
  }

  aa <- strsplit(sequence, "")[[1]]
  paste0(aa[mask], collapse = "")
}

## ── Write peptide-only FASTA ────────────────────────────────
## Headers follow the format: >{prefix}_{i}
write_peptide_fasta <- function(dt, prefix, path) {
  headers     <- paste0(">", prefix, "_", seq_len(nrow(dt)))
  fasta_lines <- as.vector(rbind(headers, dt$peptide))
  writeLines(fasta_lines, path)
  cat("FASTA written:", nrow(dt), "entries ->", path, "\n")
}

## ============================================================
## SECTION 1: Load annotated IEDB data
## ============================================================
cat("Loading annotated IEDB data...\n")
df <- fread(FULL_CSV)

## ============================================================
## SECTION 2: Subset to DQ positives
## ============================================================
dq_pos <- df[locus == "DQ" & label == 1]

cat("DQ positive records:", nrow(dq_pos), "\n")
cat("Unique source proteins:", uniqueN(dq_pos$uniprot_id), "\n\n")

## ============================================================
## SECTION 3: Build known DQ positive peptide set
## ============================================================
dq_positives <- unique(dq_pos$peptide)
cat("Known DQ positive peptides:", length(dq_positives), "\n\n")

## ============================================================
## SECTION 4: Build protein sequence lookup
## ============================================================
prot_seqs <- unique(dq_pos[, .(uniprot_id, sequence)])
prot_seqs <- prot_seqs[!is.na(sequence) & sequence != ""]

cat("DQ proteins available for sampling:", nrow(prot_seqs), "\n\n")

## ============================================================
## SECTION 5: Build per-protein binding region map
## ============================================================
binding_map <- dq_pos[, .(
  starts = list(start),
  ends   = list(end)
), by = uniprot_id]

## ============================================================
## SECTION 6: Build shuffled unmasked residue pool per protein
## ============================================================
## For each protein:
##   1. Retrieve its binding intervals
##   2. Extend intervals by FLANK_BUFFER on each side
##   3. Concatenate all residues outside those extended intervals
##   4. Shuffle the resulting residue pool
##
## The flank buffer removes the immediate N- and C-terminal
## anchor context around each binding region, ensuring sampled
## negatives carry no partial binding signal from groove-
## contacting flanking positions.
## Proteins where the entire sequence is masked are dropped.
## ============================================================

cat("Building shuffled unmasked residue pools (FLANK_BUFFER =", FLANK_BUFFER, ")...\n")

prot_seqs[, shuffled_pool := {
  row_map <- binding_map[uniprot_id == .BY$uniprot_id]
  seq1    <- sequence[1]

  if (nrow(row_map) == 0) {
    aa <- strsplit(seq1, "")[[1]]
    paste0(sample(aa), collapse = "")
  } else {
    unmasked <- get_unmasked_residues(
      seq1,
      row_map$starts[[1]],
      row_map$ends[[1]]
    )
    if (nchar(unmasked) == 0) {
      NA_character_
    } else {
      aa <- strsplit(unmasked, "")[[1]]
      paste0(sample(aa), collapse = "")
    }
  }
}, by = uniprot_id]

## Drop proteins with no unmasked residues remaining
prot_seqs <- prot_seqs[!is.na(shuffled_pool) & nchar(shuffled_pool) >= 12L]
cat("Proteins with usable shuffled pools:", nrow(prot_seqs), "\n\n")

## ============================================================
## SECTION 7: Sample negatives from shuffled pools
## ============================================================
cat("Sampling negatives from shuffled unmasked pools...\n")

negatives     <- character(0)
negatives_set <- character(0)

for (i in seq_len(nrow(prot_seqs))) {

  uid      <- prot_seqs$uniprot_id[i]
  pool     <- prot_seqs$shuffled_pool[i]
  pool_len <- nchar(pool)

  generated <- 0L
  attempts  <- 0L

  while (generated < N_PER_PROTEIN && attempts < MAX_ATTEMPTS) {
    attempts <- attempts + 1L

    pep_len <- sample(12L:25L, 1L)
    if (pool_len < pep_len) break

    s         <- sample(1L:(pool_len - pep_len + 1L), 1L)
    e         <- s + pep_len - 1L
    candidate <- substr(pool, s, e)

    if (candidate %in% dq_positives)  next
    if (candidate %in% negatives_set) next

    negatives     <- c(negatives, candidate)
    negatives_set <- c(negatives_set, candidate)
    generated     <- generated + 1L
  }

  if (generated < N_PER_PROTEIN) {
    cat(sprintf(
      "Warning: Only generated %d/%d negatives for protein %s\n",
      generated, N_PER_PROTEIN, uid
    ))
  }
}

cat("Total DQ protein-derived negatives generated:", length(negatives), "\n")

## ============================================================
## SECTION 8: Write outputs
## ============================================================
neg_dt <- data.table(peptide = negatives)

writeLines(neg_dt$peptide, OUTPUT_TXT)
cat("TXT written:", nrow(neg_dt), "peptides ->", OUTPUT_TXT, "\n")

write_peptide_fasta(neg_dt, "DQ_0", OUTPUT_FAA)

cat("\nDONE\n")

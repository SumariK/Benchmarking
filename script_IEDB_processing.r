## ============================================================
## IEDB MHC Class II Ligand Pipeline
## ============================================================
##
## PURPOSE
## -------
## Processes the IEDB MHC ligand dataset to produce per-locus
## peptide input files and protein FASTA files for use with
## MixMHC2pred. Handles both UniProt and NCBI-sourced sequences,
## recovers misrouted accessions, and produces two sets of output:
## (1) full deduplicated dataset
## (2) downsampled to a 1:5 positive-to-negative ratio
## (3) CAPTAn-ready FASTA files with guaranteed unique headers
##
## INPUT
##   mhc_ligand_full.csv — Full IEDB MHC ligand export (skip=1)
##
## DEPENDENCIES
##   data.table, httr, Biostrings, rentrez
## ============================================================

library(data.table)
library(httr)
library(Biostrings)
library(rentrez)

## ============================================================
## SECTION 1: Load and filter IEDB data
## ============================================================

setwd("/gpfs/home4/skleynhans/examples_captan")

df <- fread(
  "/gpfs/home4/skleynhans/examples_captan/mhc_ligand_full.csv",
  skip = 1
)
setnames(df, make.unique(names(df)))

df <- df[, .(
  peptide     = Name,
  allele      = `Name.6`,
  start       = as.integer(`Starting Position`),
  source_iri  = `Source Molecule IRI`,
  organism    = `Source Organism`,
  class       = Class,
  qualitative = `Qualitative Measurement`
)]

df <- df[
  class    == "II"           &
  organism == "Homo sapiens" &
  !is.na(start)
]

df <- df[!grepl("\\+", peptide)]
df <- df[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", peptide)]
df <- df[nchar(peptide) >= 12 & nchar(peptide) <= 25]
df[, end := start + nchar(peptide) - 1L]
df <- df[grepl("^HLA-(DR|DQ|DP)", allele)]

cat("Qualitative measurement breakdown:\n")
print(table(df$qualitative, useNA = "ifany"))
cat("Negatives after allele filter:", nrow(df[qualitative == "Negative"]), "\n\n")

## ============================================================
## SECTION 2: Protein ID extraction and validation
## ============================================================

df[, db_source := fcase(
  grepl("uniprot", source_iri, ignore.case = TRUE), "uniprot",
  grepl("ncbi",    source_iri, ignore.case = TRUE), "ncbi",
  default = "other"
)]

df[, raw_id := sub("\\..*", "", sub(".*/", "", source_iri))]
df[db_source == "uniprot", uniprot_id := raw_id]
df[db_source == "ncbi",    ncbi_id    := raw_id]

valid_uniprot <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$"

df[, id_valid := fcase(
  db_source == "uniprot", grepl(valid_uniprot, uniprot_id),
  db_source == "ncbi",    !is.na(ncbi_id) & ncbi_id != "",
  default = FALSE
)]

## Retain negatives even without a valid ID
df <- df[qualitative == "Negative" | id_valid == TRUE]

cat("Negatives retained after ID validation:\n")
cat("  Total:          ", nrow(df[qualitative == "Negative"]), "\n")
cat("  with UniProt ID:", nrow(df[qualitative == "Negative" & db_source == "uniprot"]), "\n")
cat("  with NCBI ID:  ", nrow(df[qualitative == "Negative" & db_source == "ncbi"]), "\n\n")

## ============================================================
## SECTION 3: Sequence retrieval
## ============================================================

fetch_ncbi_seqs <- function(ids, batch_size = 200) {
  results <- list()
  batches  <- split(ids, ceiling(seq_along(ids) / batch_size))

  for (i in seq_along(batches)) {
    cat("Fetching NCBI batch", i, "of", length(batches), "\n")
    tryCatch({
      raw   <- entrez_fetch(
        db      = "protein",
        id      = batches[[i]],
        rettype = "fasta",
        retmode = "text"
      )
      lines   <- strsplit(raw, "\n")[[1]]
      headers <- grep("^>", lines)

      for (j in seq_along(headers)) {
        start_l <- headers[j] + 1L
        end_l   <- ifelse(j < length(headers), headers[j + 1L] - 1L, length(lines))
        seq     <- paste(lines[start_l:end_l], collapse = "")
        hdr     <- lines[headers[j]]
        acc     <- sub("^>\\S+\\|(\\S+)\\|.*", "\\1", hdr)
        if (acc == hdr) acc <- sub("^>(\\S+).*", "\\1", hdr)
        results[[acc]] <- seq
      }
      Sys.sleep(0.4)
    }, error = function(e) cat("Batch", i, "failed:", conditionMessage(e), "\n"))
  }
  results
}

fetch_fasta_chunk <- function(chunk_ids) {
  tryCatch({
    res <- GET(
      "https://rest.uniprot.org/uniprotkb/accessions",
      query = list(
        accessions = paste(chunk_ids, collapse = ","),
        format     = "fasta"
      )
    )
    if (status_code(res) != 200) {
      warning("UniProt fetch returned status: ", status_code(res))
      return(NULL)
    }
    content(res, "text", encoding = "UTF-8")
  }, error = function(e) {
    warning("UniProt fetch error: ", conditionMessage(e))
    NULL
  })
}

## Fetch NCBI sequences
ncbi_ids_to_fetch <- unique(df[db_source == "ncbi" & !is.na(ncbi_id), ncbi_id])
cat("Fetching", length(ncbi_ids_to_fetch), "unique NCBI IDs...\n")

ncbi_seq_list <- fetch_ncbi_seqs(ncbi_ids_to_fetch)
ncbi_seq_dt   <- data.table(
  ncbi_id  = names(ncbi_seq_list),
  sequence = unlist(ncbi_seq_list)
)
cat("NCBI sequences retrieved:", nrow(ncbi_seq_dt), "\n")

## Recover misrouted UniProt IDs filed under NCBI IRIs
missing_ncbi  <- ncbi_ids_to_fetch[!ncbi_ids_to_fetch %in% ncbi_seq_dt$ncbi_id]
misrouted_ids <- missing_ncbi[grepl(valid_uniprot, missing_ncbi)]
cat("Misrouted UniProt-format IDs to recover:", length(misrouted_ids), "\n\n")

## Fetch UniProt sequences (including recovered IDs)
uniprot_ids <- unique(c(
  df[db_source == "uniprot" & !is.na(uniprot_id), uniprot_id],
  misrouted_ids
))

chunk_size  <- 100
id_chunks   <- split(uniprot_ids, ceiling(seq_along(uniprot_ids) / chunk_size))

cat("Fetching", length(uniprot_ids), "unique UniProt IDs (incl. misrouted)...\n")
fasta_chunks <- lapply(id_chunks, function(chunk) {
  Sys.sleep(0.2)
  fetch_fasta_chunk(chunk)
})

fasta_text <- paste(Filter(Negate(is.null), fasta_chunks), collapse = "\n")
tmp        <- tempfile(fileext = ".fasta")
writeLines(trimws(fasta_text), tmp)

fa           <- readAAStringSet(tmp)
uniprot_seqs <- as.character(fa)
names(uniprot_seqs) <- sapply(strsplit(names(fa), "\\|"), `[`, 2)
cat("UniProt sequences retrieved:", length(uniprot_seqs), "\n\n")

## Merge sequences back onto df
df[db_source == "ncbi" & ncbi_id %in% misrouted_ids, uniprot_id := ncbi_id]

df[db_source == "uniprot" | (db_source == "ncbi" & ncbi_id %in% misrouted_ids),
   sequence := uniprot_seqs[uniprot_id]]
df[db_source == "ncbi" & !ncbi_id %in% misrouted_ids,
   sequence := ncbi_seq_dt$sequence[match(ncbi_id, ncbi_seq_dt$ncbi_id)]]

df <- df[!is.na(sequence) & sequence != ""]

cat("Sequence coverage diagnostics:\n")
cat("  Negatives with sequence:", nrow(df[qualitative == "Negative"]), "\n")
cat("  Positives with sequence:", nrow(df[qualitative != "Negative"]), "\n\n")

## ============================================================
## SECTION 4: Flanking context computation
## ============================================================
## Extracts 12 characters of context per peptide:
##   3 upstream + 3 N-terminal + 3 C-terminal + 3 downstream
## Terminus positions are padded with PAD ("-" by default).
## Change PAD to "X" if your MixMHC2pred version requires it.
## ============================================================

PAD <- "-"

get_context_12 <- function(prot, start, end) {
  aa <- strsplit(prot, "")[[1]]
  n  <- length(aa)
  get <- function(i) if (i < 1L || i > n) PAD else aa[i]

  up   <- c(get(start - 3L), get(start - 2L), get(start - 1L))
  nter <- c(get(start),      get(start + 1L), get(start + 2L))
  cter <- c(get(end   - 2L), get(end   - 1L), get(end))
  down <- c(get(end   + 1L), get(end   + 2L), get(end   + 3L))

  paste0(c(up, nter, cter, down), collapse = "")
}

cat("Computing flanking context for", nrow(df), "rows...\n")

seq_vec   <- df$sequence
start_vec <- df$start
end_vec   <- df$end

df[, context12 := mapply(
  function(s, st, en) {
    if (is.null(s) || length(s) == 0L || is.na(s)) return(NA_character_)
    get_context_12(s, st, en)
  },
  seq_vec, start_vec, end_vec,
  USE.NAMES = FALSE,
  SIMPLIFY  = TRUE
)]
df[, context12 := as.character(context12)]

## ============================================================
## SECTION 5: Annotation
## ============================================================

df[, locus := sub("HLA-(DR|DQ|DP).*", "\\1", allele)]

df[, label := fifelse(
  grepl("^Positive", qualitative, ignore.case = TRUE), 1L, 0L
)]

df[, haplotype := {
  mol <- sub(".*\\*(\\d+):?.*",    "\\1", allele)
  ser <- sub(".*[A-Za-z](\\d+)$",    "\\1", allele)
  dpw <- sub(".*[A-Za-z]w(\\d+)$",   "\\1", allele)
  ifelse(grepl("\\*",     allele), mol,
  ifelse(grepl("w\\d+$",  allele), dpw,
  ifelse(grepl("\\d+$",   allele), ser, "unknown")))
}]

cat("Alleles with unknown haplotype (retained):\n")
print(unique(df[haplotype == "unknown", allele]))

df[, fasta_header := paste0(locus, "_", label, "_", haplotype)]

## ============================================================
## SECTION 6: Deduplication
## ============================================================
## Removes duplicate peptides within each locus+label group,
## keeping one representative entry per unique peptide sequence.
## ============================================================

df <- df[, .SD[!duplicated(peptide)], by = .(locus, label)]

cat("\nMissing contexts: ",  sum(is.na(df$context12)), "\n")
cat("Final rows (full): ",    nrow(df), "\n")
cat("Peptide counts per locus and label (full dataset):\n")
print(df[, .N, by = .(locus, label)][order(locus, label)])

df_full <- copy(df)

## ============================================================
## SECTION 7: Downsampling (1:5 positive-to-negative ratio)
## ============================================================
## Positives are downsampled to floor(n_neg / 5) per locus.
## Negatives are never downsampled.
## set.seed(42) ensures reproducibility.
## ============================================================

set.seed(42)

df_ds <- rbindlist(lapply(c("DR", "DQ", "DP"), function(loc) {
  negs      <- df_full[locus == loc & label == 0L]
  pos       <- df_full[locus == loc & label == 1L]
  n_pos_max <- floor(nrow(negs) / 5L)

  cat(loc,
      "— negatives:", nrow(negs),
      "| positives before:", nrow(pos),
      "| positives after:", n_pos_max, "\n")

  if (nrow(pos) > n_pos_max) pos <- pos[sample(.N, n_pos_max)]
  rbind(negs, pos)
}))

cat("\nFinal counts after downsampling:\n")
print(df_ds[, .N, by = .(locus, label)][order(locus, label)])
cat("Downsampled rows: ", nrow(df_ds), "\n\n")

## ============================================================
## SECTION 8: Build combined sequence lookup
## ============================================================
## Merges UniProt and NCBI sequence vectors into a single named
## lookup used when writing FASTA files in Sections 9 and 10.
## ============================================================

df_full[, seq_key := fifelse(!is.na(uniprot_id), uniprot_id, ncbi_id)]
df_ds[,  seq_key := fifelse(!is.na(uniprot_id), uniprot_id, ncbi_id)]

all_seqs <- c(
  uniprot_seqs,
  setNames(ncbi_seq_dt$sequence, ncbi_seq_dt$ncbi_id)
)
all_seqs <- all_seqs[!duplicated(names(all_seqs))]

## ============================================================
## SECTION 9: Save outputs
## ============================================================

out_dir <- "/gpfs/home4/skleynhans/report/iedb_data"

write_locus_files <- function(data, loc, lbl, tag, suffix = "") {
  ## --- Peptide + context text file ---
  pep_dt <- data[locus == loc & label == lbl, .(peptide, context12)]
  if (nrow(pep_dt) > 0) {
    fwrite(
      pep_dt,
      file.path(out_dir, paste0(loc, "_", tag, suffix, ".txt")),
      sep = "\t", col.names = FALSE
    )
    cat(loc, tag, suffix, "written:", nrow(pep_dt), "peptides\n")
  } else {
    cat(loc, tag, suffix, ": no rows found, file not written\n")
  }

  ## --- Source protein FASTA ---
  df_sub <- unique(data[locus == loc & label == lbl], by = "seq_key")
  if (nrow(df_sub) > 0) {
    fa_lines <- unlist(mapply(
      function(header, key) {
        prot <- all_seqs[[key]]
        if (is.null(prot) || prot == "") return(NULL)
        c(paste0(">", header), prot)
      },
      df_sub$fasta_header, df_sub$seq_key,
      SIMPLIFY = FALSE
    ))
    writeLines(fa_lines, file.path(out_dir, paste0(loc, "_", tag, suffix, ".faa")))
    cat(loc, tag, suffix, "FASTA written:", nrow(df_sub), "unique proteins\n")
  } else {
    cat(loc, tag, suffix, ": no rows found, FASTA not written\n")
  }
}

## ---- Full dataset ----------------------------------------
cat("\n--- Writing full dataset outputs ---\n")
fwrite(df_full, file.path(out_dir, "context_iedb_full.csv"))
cat("Full annotated table saved\n")

fwrite(
  df_full[, .(peptide, context12)],
  file.path(out_dir, "mixmhc_context_all.txt"),
  sep = "\t", col.names = FALSE
)

df_fasta    <- unique(df_full, by = "seq_key")
fasta_lines <- unlist(mapply(
  function(header, key) {
    prot <- all_seqs[[key]]
    if (is.null(prot) || prot == "") return(NULL)
    c(paste0(">", header), prot)
  },
  df_fasta$fasta_header, df_fasta$seq_key,
  SIMPLIFY = FALSE
))
writeLines(fasta_lines, file.path(out_dir, "iedb_proteins_all.faa"))
cat("Full protein FASTA written:", nrow(df_fasta), "unique proteins\n")

for (loc in c("DR", "DQ", "DP")) {
  write_locus_files(df_full, loc, 1L, "positives", suffix = "")
  write_locus_files(df_full, loc, 0L, "negatives", suffix = "")
}

## ---- Downsampled dataset ---------------------------------
cat("\n--- Writing downsampled dataset outputs ---\n")
fwrite(df_ds, file.path(out_dir, "context_iedb_downsampled.csv"))
cat("Downsampled annotated table saved\n")

fwrite(
  df_ds[, .(peptide, context12)],
  file.path(out_dir, "mixmhc_context_all_downsampled.txt"),
  sep = "\t", col.names = FALSE
)

df_fasta_ds    <- unique(df_ds, by = "seq_key")
fasta_lines_ds <- unlist(mapply(
  function(header, key) {
    prot <- all_seqs[[key]]
    if (is.null(prot) || prot == "") return(NULL)
    c(paste0(">", header), prot)
  },
  df_fasta_ds$fasta_header, df_fasta_ds$seq_key,
  SIMPLIFY = FALSE
))
writeLines(fasta_lines_ds, file.path(out_dir, "iedb_proteins_all_downsampled.faa"))
cat("Downsampled protein FASTA written:", nrow(df_fasta_ds), "unique proteins\n")

for (loc in c("DR", "DQ", "DP")) {
  write_locus_files(df_ds, loc, 1L, "positives", suffix = "_downsampled")
  write_locus_files(df_ds, loc, 0L, "negatives", suffix = "_downsampled")
}

## ============================================================
## SECTION 10: CAPTAn-ready FASTA files with unique headers
## ============================================================
## Generates per-locus protein FASTA files from the downsampled
## dataset with guaranteed unique headers for CAPTAn submission.
## Also saves a header mapping CSV for use in evaluation.
##
## Header format: {locus}_{label}_{haplotype}_{i}
##   e.g. DQ_1_03_1, DQ_0_unknown_12
## The trailing index i guarantees uniqueness within each
## locus+label group regardless of haplotype collisions.
## ============================================================

cat("\n--- Writing CAPTAn-ready FASTA files ---\n")

write_captan_fasta <- function(data, loc, lbl, out_dir) {
  lbl_str <- ifelse(lbl == 1L, "1", "0")
  prefix  <- paste0(loc, "_", lbl_str)

  ## Deduplicate by sequence — one FASTA entry per unique protein
  sub <- unique(
    data[locus == loc & label == lbl & !is.na(sequence) & sequence != ""],
    by = "sequence"
  )

  if (nrow(sub) == 0) {
    cat(prefix, ": no rows found, CAPTAn FASTA not written\n")
    return(NULL)
  }

  ## Unique header: {locus}_{label}_{haplotype}_{i}
  sub[, captan_header := paste0(prefix, "_", haplotype, "_", seq_len(.N))]
  stopifnot(uniqueN(sub$captan_header) == nrow(sub))

  ## Write FASTA
  fasta_path <- file.path(out_dir, paste0(prefix, "_captan.faa"))
  fa_lines   <- unlist(mapply(
    function(header, seq) c(paste0(">", header), seq),
    sub$captan_header,
    sub$sequence,
    SIMPLIFY = FALSE
  ))
  writeLines(fa_lines, fasta_path)
  cat(prefix, "CAPTAn FASTA written:", nrow(sub),
      "unique proteins ->", fasta_path, "\n")

  ## Return old -> new header mapping
  return(sub[, .(
    peptide,
    sequence,
    locus,
    label,
    haplotype,
    old_header    = fasta_header,
    captan_header
  )])
}

## Generate CAPTAn FASTA and header mappings for all loci
all_maps <- list()

for (loc in c("DR", "DQ", "DP")) {
  pos_map <- write_captan_fasta(df_ds, loc, 1L, out_dir)
  neg_map <- write_captan_fasta(df_ds, loc, 0L, out_dir)
  loc_map <- rbind(pos_map, neg_map)

  if (!is.null(loc_map) && nrow(loc_map) > 0) {
    map_path <- file.path(out_dir, paste0(loc, "_header_mapping.csv"))
    fwrite(loc_map, map_path)
    cat(loc, "header mapping saved ->", map_path, "\n")
    all_maps[[loc]] <- loc_map
  }
}

## Save combined mapping across all loci
combined_map <- rbindlist(all_maps)
fwrite(combined_map, file.path(out_dir, "all_loci_header_mapping.csv"))
cat("\nCombined header mapping saved ->",
    file.path(out_dir, "all_loci_header_mapping.csv"), "\n")

## Final uniqueness check
cat("\nCAPTAn header uniqueness check:\n")
for (loc in c("DR", "DQ", "DP")) {
  loc_map <- all_maps[[loc]]
  if (!is.null(loc_map))
    cat(loc, "- unique headers:",
        uniqueN(loc_map$captan_header) == nrow(loc_map), "\n")
}

cat("\nALL DONE\n")
cat("Submit *_captan.faa files to CAPTAn\n")
cat("Use *_header_mapping.csv files to match predictions back to peptides\n")

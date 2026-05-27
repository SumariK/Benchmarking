## ============================================================
## IEDB MHC Class II — Negatives Processing Pipeline
## ============================================================
##
## PURPOSE
## -------
## Filters IEDB MHC ligand data to human Class II negatives,
## retrieves source protein sequences from UniProt, builds an
## annotated FASTA, and produces cleaned output CSVs.
##
## INPUT
##   mhc_ligand_full.csv — Full IEDB MHC ligand export (skip=1)
##
## OUTPUT
##   protein_with_peptides.fasta — Annotated protein FASTA
##   peptide_mapping.csv         — UniProt ID + peptide positions
##   sel_cols.csv                — Column-selected filtered table
##
## DEPENDENCIES
##   data.table, httr, Biostrings, readr
## ============================================================

library(data.table)
library(httr)
library(Biostrings)
library(readr)

## ============================================================
## SECTION 1: Load and filter data
## ============================================================
## Keeps only: human · Class II · Negative · valid start position
## Deduplicates immediately after loading.
## ============================================================

df <- fread(
  "/gpfs/home4/skleynhans/examples_captan/mhc_ligand_full.csv",
  skip = 1
)

df <- df[
  Class                  == "II"           &
  `Source Organism`      == "Homo sapiens" &
  `Qualitative Measurement` == "Negative"  &
  !is.na(`Starting Position`)
]

df <- unique(df)

## ============================================================
## SECTION 2: Extract and validate UniProt IDs
## ============================================================
## Parses accessions from the Source Molecule IRI column.
## Strips the base URL and any version suffixes (e.g. ".1").
## Only IDs matching the canonical UniProt format are kept.
## ============================================================

df[, uniprot_id := sub(".*/",  "", `Source Molecule IRI`)]
df[, uniprot_id := sub("\..*", "", uniprot_id)]
df <- df[!is.na(uniprot_id) & uniprot_id != ""]

ids <- unique(df$uniprot_id)
ids <- trimws(sub("\..*", "", ids))

## Keep only valid UniProt accession format
valid_uniprot <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$"
ids <- ids[grepl(valid_uniprot, ids)]

## ============================================================
## SECTION 3: Clean peptide info
## ============================================================
## Renames the peptide column, removes PTM annotations (+),
## non-standard amino acid characters, and empty entries.
## Computes start/end positions and deduplicates.
## ============================================================

df[, peptide := Name]
df[, allele  := `Name.6`]

## Drop empty peptides
df <- df[nchar(peptide) > 0]

## Remove PTM-modified peptides
df <- df[!grepl("\+", peptide)]

## Keep only standard amino acid characters
df <- df[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", peptide)]

## Compute positions
df[, start := as.integer(`Starting Position`)]
df[, end   := start + nchar(peptide) - 1L]

## Deduplicate by protein + peptide + position
df <- unique(df, by = c("uniprot_id", "peptide", "start", "end"))

## ============================================================
## SECTION 4: Fetch protein sequences from UniProt (bulk)
## ============================================================
## Splits UniProt IDs into chunks of 100 and queries the
## UniProt REST API for FASTA sequences. Responses are
## concatenated, written to a temp file, and parsed with
## Biostrings. UniProt IDs are extracted from |ID| in headers.
## ============================================================

chunk_size <- 100
id_chunks  <- split(ids, ceiling(seq_along(ids) / chunk_size))

fetch_fasta_chunk <- function(chunk_ids) {
  res <- GET(
    "https://rest.uniprot.org/uniprotkb/accessions",
    query = list(
      accessions = paste(chunk_ids, collapse = ","),
      format     = "fasta"
    )
  )
  if (status_code(res) != 200) {
    cat("Chunk failed:", status_code(res), "\n")
    cat(content(res, "text"), "\n")
    return(NULL)
  }
  content(res, "text", encoding = "UTF-8")
}

cat("Fetching FASTA in chunks...\n")

fasta_chunks <- lapply(seq_along(id_chunks), function(i) {
  cat("Chunk", i, "of", length(id_chunks), "\n")
  Sys.sleep(0.2)
  fetch_fasta_chunk(id_chunks[[i]])
})

fasta_text <- paste(
  unlist(fasta_chunks[!sapply(fasta_chunks, is.null)]),
  collapse = "\n"
)

tmp <- tempfile(fileext = ".fasta")
writeLines(fasta_text, tmp)

fa          <- readAAStringSet(tmp)
seqs        <- as.list(as.character(fa))
clean_names <- sub(".*\|([^|]+)\|.*", "\1", names(fa))
names(seqs) <- clean_names

## ============================================================
## SECTION 5: Build protein lookup table
## ============================================================
## Converts parsed sequences into a named list for fast lookup.
## Prints a sanity check: how many df UniProt IDs were found.
## ============================================================

proteins <- lapply(names(seqs), function(id) {
  list(header = id, seq = seqs[[id]])
})
names(proteins) <- names(seqs)

cat("Example protein IDs:\n")
print(head(names(proteins)))

cat("Example df IDs:\n")
print(head(df$uniprot_id))

cat("Matches:\n")
print(sum(df$uniprot_id %in% names(proteins)))

## ============================================================
## SECTION 6: Build annotated FASTA output
## ============================================================
## For each retrieved protein, writes:
##   >UNIPROT_ID
##   <full sequence>
##   # peptide=<seq> start=<n> end=<n>   (one line per peptide)
## ============================================================

out_file <- "protein_with_peptides.fasta"
file.create(out_file)

for (id in names(proteins)) {

  prot <- proteins[[id]]
  seq  <- prot$seq

  pep_df <- df[uniprot_id == id]
  pep_df <- unique(pep_df, by = c("peptide", "start", "end"))

  ## Protein header and sequence
  cat(paste0(">", id, "\n"),  file = out_file, append = TRUE)
  cat(paste0(seq,  "\n"),     file = out_file, append = TRUE)

  ## Peptide position annotations
  if (nrow(pep_df) > 0) {
    for (i in seq_len(nrow(pep_df))) {
      cat(
        paste0(
          "# peptide=", pep_df$peptide[i],
          " start=",    pep_df$start[i],
          " end=",      pep_df$end[i], "\n"
        ),
        file = out_file, append = TRUE
      )
    }
  }
}

cat("Annotated FASTA written to:", out_file, "\n")

## ============================================================
## SECTION 7: Save peptide mapping table
## ============================================================
## Writes a four-column CSV:
##   uniprot_id | peptide | start | end
## ============================================================

fwrite(
  df[, .(uniprot_id, peptide, start, end)],
  "peptide_mapping.csv"
)

cat("Peptide mapping table written to: peptide_mapping.csv\n")
cat("DONE\n")
getwd()

## ============================================================
## SECTION 8: Column selection — save filtered output CSV
## ============================================================
## Selects only the columns needed downstream and writes a
## clean CSV. The helper function validates column existence
## before writing and raises a clear error if any are missing.
## ============================================================

## Helper: read a CSV, validate columns, write a subset CSV
create_csv_with_selected_columns <- function(input_file, output_file, selected_columns) {
  tryCatch({
    if (!file.exists(input_file))
      stop(paste("File not found:", input_file))

    data <- read_csv(input_file, show_col_types = FALSE)

    missing_cols <- setdiff(selected_columns, colnames(data))
    if (length(missing_cols) > 0)
      stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

    write_csv(data[selected_columns], output_file)
    message("CSV created: ", output_file)

  }, error = function(e) message("Error: ", e$message))
}

## Select and save columns directly from the in-memory df
columns_to_keep <- c("Name.6", "uniprot_id", "peptide")

missing_cols <- setdiff(columns_to_keep, colnames(df))
if (length(missing_cols) > 0)
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

df_selected <- df[, ..columns_to_keep]

fwrite(
  df_selected,
  "/gpfs/home4/skleynhans/examples_captan/sel_cols.csv"
)

cat("Column-selected CSV written to: sel_cols.csv\n")

 ## ============================================================
## UTILITY: merge_ranges(starts, ends)
## ============================================================
## Merges overlapping integer ranges into non-overlapping
## intervals. Useful for masking peptide-covered regions
## before shuffling or sampling remaining sequence.
##
## Example:
##   merge_ranges(c(10, 15, 30), c(20, 25, 40))
##   # Returns merged intervals: [10,25], [30,40]
## ============================================================

merge_ranges <- function(starts, ends) {

  ## Remove NAs
  valid  <- !is.na(starts) & !is.na(ends)
  starts <- starts[valid]
  ends   <- ends[valid]

  if (length(starts) == 0)
    return(data.table(start = integer(), end = integer()))

  ## Sort by start position
  ord    <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]

  merged_starts <- c()
  merged_ends   <- c()
  cur_start     <- starts[1]
  cur_end       <- ends[1]

  for (i in 2:length(starts)) {
    if (starts[i] <= cur_end) {
      cur_end <- max(cur_end, ends[i])
    } else {
      merged_starts <- c(merged_starts, cur_start)
      merged_ends   <- c(merged_ends,   cur_end)
      cur_start     <- starts[i]
      cur_end       <- ends[i]
    }
  }

  merged_starts <- c(merged_starts, cur_start)
  merged_ends   <- c(merged_ends,   cur_end)

  data.table(start = merged_starts, end = merged_ends)
}          

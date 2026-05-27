## ============================================================
## MixMHC2pred Performance Assessment — DQ / DR / DP
## ============================================================
##
## PURPOSE
## -------
## Evaluates MixMHC2pred binding predictions per locus.
## Labels are assigned from explicit row ranges — no label
## file needed. Blacklisted peptides are removed before any
## metrics are computed.
##
## Per locus:
##   1. Load predictions + assign labels from row ranges
##   2. Remove blacklisted peptides + clean %Rank_best
##   3. %Rank_best distribution histogram (pos vs neg)
##   4. ROC curve + AUC
##   5. PR  curve + AP
##   6. Summary table + save
##
## OUTPUT (written to OUT_DIR):
##   {locus}_labelled.csv   — predictions + Label column
##   {locus}_Rank_Dist.png  — %Rank_best histogram
##   {locus}_ROC.png        — ROC curve
##   {locus}_PR.png         — PR curve
##   allele_metrics_summary.csv — AUC + AP + N table
##
## DEPENDENCIES
##   data.table, ggplot2, pROC, PRROC
## ============================================================

library(data.table)
library(ggplot2)
library(pROC)
library(PRROC)

## ============================================================
## CONFIG: Output directory and locus definitions
## ============================================================
## For each locus, specify:
##   pred_file — path to MixMHC2pred output TSV
##   neg_rows  — 1-based row indices for negative peptides
##   pos_rows  — 1-based row indices for positive peptides
##
## Row ranges accept: continuous (1:1209), vectors (c(1,2,5)),
## or mixed (c(1:500, 600:700)).
##
## DR and DP are commented out — uncomment and fill in paths
## and row ranges to activate them.
## ============================================================

OUT_DIR <- "/zfs/omics/personal/15843580/Internship/mixmhc_results/no_ctx_iedb"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOCI <- list(

  DQ = list(
    pred_file = "/zfs/omics/personal/15843580/Internship/mixmhc_results/dq_noctx_iedb_results.txt",
    neg_rows  = 1:228,    ## <-- edit these numbers
    pos_rows  = 229:273   ## <-- edit these numbers
  )

  # DR = list(
  #   pred_file = "/gpfs/home4/skleynhans/benchmarking/mix_mhc/mixmhc_hr_res.txt",
  #   neg_rows  = 1:1209,
  #   pos_rows  = 1210:2418
  # ),
  #
  # DP = list(
  #   pred_file = "/gpfs/home4/skleynhans/benchmarking/mix_mhc/mixmhc_dp_res.txt",
  #   neg_rows  = 1:600,
  #   pos_rows  = 601:1200
  # )

)

## ============================================================
## CONFIG: Blacklisted peptides
## ============================================================
## Negative peptides confirmed to overlap with positives in
## the UniProt-derived dataset. Excluded before any metrics
## are computed. Add sequences here as new overlaps are found.
## ============================================================

BLACKLIST <- list(

  DP = c(
    "ANITPFGADQVKDRGPE",
    "ADLLDYIKALNRNSDR",
    "KDPEGLFLQDNIVA",
    "VILAKGAEEMETVIPVD",
    "FPDSLIVKGFNVVSAW",
    "TSLKDYIKAEEDKLEQ",
    "VENSFFLNVNSQ",
    "IMNSFVNDIFERIASE"
  ),

  DR = c(
    "IDSVIVVDNVPQVG",
    "FDASKIKIEFTPEQIEEF",
    "ADLLDYIKALNRNSDR",
    "QQRLIFAGKQLED",
    "KSVLQAVQKTDEGHPF",
    "KLEFSIYPAPQVSTA",
    "IASVAGLTAAAYR",
    "LDDFHVNGGELILIH",
    "LKAMLVFAEHRYYGE",
    "AVLSLYASGRTTG",
    "VPIYEGYALPHAILRLD",
    "GEALGRLLVVYPWTQRF",
    "KETWKALEALVAKG",
    "TARYMSINTHLGKEQ"
  ),

  DQ = c(
    "VPIYEGYALPHAILRLD",
    "IGLDAAALPGQPME",
    "YRGAAGALLVYDITRR"
  )

)

## ============================================================
## HELPERS
## ============================================================

## ── Save a ggplot to OUT_DIR at 150 dpi ─────────────────────
save_plot <- function(p, filename, width = 7, height = 5) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, plot = p, width = width, height = height, dpi = 150)
  cat("  Plot saved ->", path, "\n")
}

## ── Load MixMHC2pred TSV, skipping # comment lines ──────────
load_pred <- function(path) {
  raw   <- readLines(path)
  lines <- raw[!grepl("^#", raw)]
  fread(text = paste(lines, collapse = "\n"), sep = "\t", header = TRUE)
}

## ============================================================
## PER-LOCUS STEP 1: Load predictions and assign labels
## ============================================================
## Reads the prediction TSV (skipping # lines), assigns
## Label = 0 to neg_rows and Label = 1 to pos_rows.
## Validates that all row numbers are within file bounds and
## stops with a clear error message if not.
## Rows not covered by either range are flagged and dropped.
## ============================================================

metrics_rows <- list()   ## AUC + AP per locus

for (locus in names(LOCI)) {

  cfg <- LOCI[[locus]]

  cat("\n", strrep("=", 60), "\n")
  cat("Locus:", locus, "\n")
  cat(strrep("=", 60), "\n")

  ## Check file exists
  if (!file.exists(cfg$pred_file)) {
    cat("  WARNING: prediction file not found —", cfg$pred_file, "\n")
    cat("  Skipping locus", locus, "\n")
    next
  }

  pred <- load_pred(cfg$pred_file)
  cat("  Rows loaded:", nrow(pred), "\n")

  ## Validate row numbers are within bounds
  bad_neg <- cfg$neg_rows[cfg$neg_rows < 1L | cfg$neg_rows > nrow(pred)]
  bad_pos <- cfg$pos_rows[cfg$pos_rows < 1L | cfg$pos_rows > nrow(pred)]
  if (length(bad_neg) > 0 || length(bad_pos) > 0) {
    stop(sprintf(
      "Row numbers out of bounds for locus %s (file has %d rows).\n  Bad neg: %s\n  Bad pos: %s",
      locus, nrow(pred),
      paste(bad_neg, collapse = ", "),
      paste(bad_pos, collapse = ", ")
    ))
  }

  pred[, Label := NA_integer_]
  pred[cfg$neg_rows, Label := 0L]
  pred[cfg$pos_rows, Label := 1L]

  unlabelled <- sum(is.na(pred$Label))
  if (unlabelled > 0)
    cat("  NOTE:", unlabelled,
        "row(s) not covered by neg_rows or pos_rows — will be dropped\n")

  cat("  Label distribution:\n")
  print(pred[, .N, by = Label])

  fwrite(pred, file.path(OUT_DIR, paste0(locus, "_labelled.csv")), sep = "\t")

  ## ============================================================
  ## PER-LOCUS STEP 2: Remove blacklisted peptides and clean
  ## ============================================================
  ## Removes any peptide in BLACKLIST[[locus]] from the table.
  ## Reports sequences found and removed, and flags any
  ## blacklisted sequences not present in the file.
  ##
  ## Then coerces %Rank_best to numeric and drops rows where
  ## it is NA or Label is NA.
  ## ============================================================

  blacklist_locus <- BLACKLIST[[locus]]

  if (!is.null(blacklist_locus) && length(blacklist_locus) > 0) {
    n_before  <- nrow(pred)
    found     <- pred$Peptide[pred$Peptide %in% blacklist_locus]
    pred      <- pred[!Peptide %in% blacklist_locus]
    n_removed <- n_before - nrow(pred)
    cat("  Blacklist — removed:", n_removed, "peptide(s)\n")
    if (length(found) > 0) {
      cat("  Removed sequences:\n")
      for (pep in found) cat("   ", pep, "\n")
    }
    not_found <- blacklist_locus[!blacklist_locus %in% found]
    if (length(not_found) > 0) {
      cat("  NOTE:", length(not_found),
          "blacklisted peptide(s) not found in file:\n")
      for (pep in not_found) cat("   ", pep, "\n")
    }
  } else {
    cat("  Blacklist — no entries for this locus\n")
  }

  pred[, `%Rank_best` := suppressWarnings(as.numeric(`%Rank_best`))]
  df_clean <- pred[!is.na(`%Rank_best`) & !is.na(Label)]
  cat("  Rows after cleaning:", nrow(df_clean), "\n")

  df_clean[, Label_f := factor(Label, levels = c(0, 1),
                               labels = c("Negative", "Positive"))]

  ## ============================================================
  ## PER-LOCUS STEP 3: %Rank_best distribution histogram
  ## ============================================================
  ## Overlaid histogram of %Rank_best for positives and negatives.
  ## Lower %Rank_best = stronger predicted binder.
  ## Saved as {locus}_Rank_Dist.png
  ## ============================================================

  p_hist <- ggplot(df_clean, aes(x = `%Rank_best`, fill = Label_f)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("Negative" = "#4C72B0",
                                 "Positive" = "#DD8452")) +
    labs(
      title = paste0("MixMHC2pred %Rank Distribution: HLA-", locus),
      x     = "%Rank_best",
      y     = "Count",
      fill  = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")

  save_plot(p_hist, paste0(locus, "_Rank_Dist.png"))

   ## ============================================================
  ## PER-LOCUS STEP 4: ROC curve
  ## ============================================================
  ## %Rank_best is negated so lower rank = higher score = better
  ## binder, as required by pROC::roc().
  ## Saved as {locus}_ROC.png
  ## ============================================================

  roc_obj <- roc(
    response  = df_clean$Label,
    predictor = -df_clean$`%Rank_best`,
    quiet     = TRUE
  )
  roc_auc <- as.numeric(auc(roc_obj))
  cat("  AUC:", round(roc_auc, 4), "\n")

  roc_dt <- data.table(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )

  p_roc <- ggplot(roc_dt, aes(x = FPR, y = TPR)) +
    geom_line(colour = "#4C72B0", linewidth = 0.9) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    annotate("text", x = 0.65, y = 0.08,
             label = paste0("AUC = ", round(roc_auc, 3)), size = 4) +
    labs(
      title = paste0("ROC Curve – MixMHC2pred: HLA-", locus),
      x     = "False Positive Rate",
      y     = "True Positive Rate"
    ) +
    theme_bw()

  save_plot(p_roc, paste0(locus, "_ROC.png"))

    ## ============================================================
  ## PER-LOCUS STEP 5: Precision-Recall curve
  ## ============================================================
  ## %Rank_best is negated so lower rank = higher score.
  ## Saved as {locus}_PR.png
  ## ============================================================

  pr_obj <- pr.curve(
    scores.class0 = -df_clean[Label == 1]$`%Rank_best`,
    scores.class1 = -df_clean[Label == 0]$`%Rank_best`,
    curve         = TRUE
  )
  pr_auc <- pr_obj$auc.integral
  cat("  PR-AUC:", round(pr_auc, 4), "\n")

  pr_dt <- as.data.table(pr_obj$curve)
  setnames(pr_dt, c("Recall", "Precision", "threshold"))

  p_pr <- ggplot(pr_dt, aes(x = Recall, y = Precision)) +
    geom_line(colour = "#DD8452", linewidth = 0.9) +
    annotate("text", x = 0.65, y = 0.92,
             label = paste0("AP = ", round(pr_auc, 3)), size = 4) +
    labs(
      title = paste0("Precision-Recall Curve – MixMHC2pred: HLA-", locus),
      x     = "Recall (Sensitivity)",
      y     = "Precision (PPV)"
    ) +
    theme_bw()

  save_plot(p_pr, paste0(locus, "_PR.png"))

  ## Accumulate metrics
  metrics_rows[[locus]] <- data.table(
    Locus = locus,
    AUC   = round(roc_auc, 4),
    AP    = round(pr_auc,  4),
    N_pos = sum(df_clean$Label == 1L),
    N_neg = sum(df_clean$Label == 0L)
  )

}  ## end per-locus loop

## ============================================================
## SUMMARY TABLE
## ============================================================
## Prints and saves AUC, AP, N_pos, N_neg per locus.
## Saved as allele_metrics_summary.csv
## ============================================================

metrics_dt <- rbindlist(metrics_rows)

cat("\n", strrep("=", 60), "\n")
cat("PERFORMANCE SUMMARY\n")
cat(strrep("=", 60), "\n")
print(metrics_dt)

fwrite(
  metrics_dt,
  file.path(OUT_DIR, "allele_metrics_summary.csv"),
  sep = "\t"
)
cat("\nSummary table ->",
    file.path(OUT_DIR, "allele_metrics_summary.csv"), "\n")

cat("\nDONE\n")

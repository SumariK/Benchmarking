## ============================================================
## CAPTAn DQ Evaluation Pipeline
## ============================================================
##
## PURPOSE
## -------
## Evaluates CAPTAn MHC Class II binding predictions against
## IEDB ground truth labels for the DQ locus.
## Recovers original fasta_headers by stripping the trailing
## index from CAPTAn's unique headers — no mapping file needed.
##
## INPUT
##   DQ/summary_CAPTAn_core.csv       — CAPTAn core predictions
##   context_iedb_downsampled.csv      — IEDB ground truth labels
##
## OUTPUT
##   captan_dq_evaluation_results.csv  — Full matched results
##   captan_dq_score_dist.png          — Score density by label
##   captan_dq_scorecore_dist.png      — ScoreCore density by label
##   captan_dq_roc.png                 — ROC curve
##   captan_dq_pr_curve.png            — PR curve
##   captan_dq_score_vs_scorecore.png  — Score vs ScoreCore scatter
##
## DEPENDENCIES
##   data.table, ggplot2, pROC, PRROC
## ============================================================

library(data.table)
library(ggplot2)
library(pROC)
library(PRROC)

## ── Paths ──────────────────────────────────────────────────
BASE_DIR   <- "/zfs/omics/personal/15843580/Internship/captan_iedb/"
CAPTAN_CSV <- file.path(BASE_DIR, "DQ/summary_CAPTAn_core.csv")
GT_CSV     <- file.path(BASE_DIR, "context_iedb_downsampled.csv")

## ============================================================
## SECTION 1: Load CAPTAn output
## ============================================================
## Loads the CAPTAn CSV and classifies each row by allele type.
## Filters to DQ proteins scored against DQ alleles only:
##   Protein must start with "DQ_"
##   Allele must contain "DQA" or "DQB"
## ============================================================

dt <- fread(CAPTAN_CSV)

dt[, allele_type := fcase(
  grepl("DQA|DQB", Allele), "DQ",
  grepl("DPA|DPB", Allele), "DP",
  grepl("DRB",     Allele), "DR",
  default = "other"
)]

dq <- dt[
  grepl("^DQ_",    Protein) &
  grepl("DQA|DQB", Allele)
]

cat("DQ rows after filtering:", nrow(dq), "\n")
cat("Unique DQ proteins:     ", uniqueN(dq$Protein), "\n\n")

## ============================================================
## SECTION 2: Recover original fasta_header from captan_header
## ============================================================
## CAPTAn appends a trailing _{i} index to guarantee unique
## FASTA headers. Stripping it recovers the original header:
##
##   DQ_1_03_7  →  DQ_1_03
##
## No separate header mapping file is required.
## ============================================================

dq[, fasta_header := sub("_\\d+$", "", Protein)]

cat("Example header recovery:\n")
print(head(unique(dq[, .(Protein, fasta_header)]), 10))

## ============================================================
## SECTION 3: Load ground truth
## ============================================================
## Loads the downsampled IEDB CSV, filters to DQ locus, and
## retains the columns needed for matching and evaluation.
## Prints both header sets for a quick visual sanity check.
## ============================================================

gt    <- fread(GT_CSV)
gt_dq <- gt[locus == "DQ", .(peptide, start, end, label, fasta_header)]

cat("\nGround truth DQ peptides:", nrow(gt_dq), "\n")
print(gt_dq[, .N, by = label])

cat("\nGround truth fasta headers:\n")
print(unique(gt_dq$fasta_header))

cat("\nRecovered fasta headers from CAPTAn:\n")
print(unique(dq$fasta_header))

## ============================================================
## SECTION 4: Collapse to best DQ allele score per window
## ============================================================
## For each unique (fasta_header, Start, End) combination,
## takes the maximum Score, ScoreCore, and ScoreContext across
## all DQ alleles and records which allele achieved the best.
## Collapses the many-alleles-per-window structure to one row
## per protein window for downstream matching.
## ============================================================

dq_best <- dq[, .(
  Score        = max(Score),
  ScoreCore    = max(ScoreCore),
  ScoreContext = max(ScoreContext),
  best_allele  = Allele[which.max(Score)]
), by = .(fasta_header, Start, End)]

cat("\nUnique protein-window combinations:", nrow(dq_best), "\n\n")

## ============================================================
## SECTION 5: Match ground truth peptides to CAPTAn windows
## ============================================================
## For each ground truth peptide, finds all CAPTAn windows on
## the same protein (fasta_header match) that positionally
## overlap the peptide's [start, end] interval.
## Selects the highest-scoring overlapping window.
##
## Three match outcomes:
##   matched               — overlapping window found
##   no_overlapping_window — protein found, no overlap
##   protein_not_found     — no CAPTAn rows for this protein
##
## Unmatched peptides (NA score) are set to Score = 0.
## ============================================================

matches <- gt_dq[, {
  captan_rows <- dq_best[fasta_header == fasta_header]

  if (nrow(captan_rows) == 0) {
    list(
      Score        = NA_real_,
      ScoreCore    = NA_real_,
      ScoreContext = NA_real_,
      best_allele  = NA_character_,
      match_type   = "protein_not_found"
    )
  } else {
    overlapping <- captan_rows[Start <= end & End >= start]

    if (nrow(overlapping) == 0) {
      list(
        Score        = NA_real_,
        ScoreCore    = NA_real_,
        ScoreContext = NA_real_,
        best_allele  = NA_character_,
        match_type   = "no_overlapping_window"
      )
    } else {
      best <- overlapping[which.max(Score)]
      list(
        Score        = best$Score,
        ScoreCore    = best$ScoreCore,
        ScoreContext = best$ScoreContext,
        best_allele  = best$best_allele,
        match_type   = "matched"
      )
    }
  }
}, by = .(peptide, start, end, label, fasta_header)]

## Unmatched peptides → Score = 0
matches[is.na(Score),        Score        := 0]
matches[is.na(ScoreCore),    ScoreCore    := 0]
matches[is.na(ScoreContext), ScoreContext := 0]

cat("──────────────── MATCH REPORT ────────────────\n")
cat("Ground truth peptides:  ", nrow(gt_dq), "\n")
print(matches[, .N, by = .(label, match_type)][order(label, match_type)])

## ============================================================
## SECTION 6: Evaluation
## ============================================================
## Computes ROC-AUC with 95% CI using pROC::roc() and
## PR-AUC using PRROC::pr.curve(). Prints both to console.
## ============================================================

roc_obj <- roc(matches$label, matches$Score, quiet = TRUE)
pr_obj  <- pr.curve(
  scores.class0 = matches$Score[matches$label == 1],
  scores.class1 = matches$Score[matches$label == 0],
  curve         = TRUE
)

cat("\n──────────────── PERFORMANCE ────────────────\n")
cat("ROC-AUC: ", round(auc(roc_obj), 3), "\n")
cat("95% CI:  ", round(ci(roc_obj)[1], 3), "-", round(ci(roc_obj)[3], 3), "\n")
cat("PR-AUC:  ", round(pr_obj$auc.integral, 3), "\n\n")

## ============================================================
## SECTION 7: Score summary table
## ============================================================
## Prints per-label descriptive statistics for Score:
##   N, Min, Median, Mean, Max, SD
## ============================================================

cat("──────────────── SCORE SUMMARY ────────────────\n")
score_summary <- matches[, .(
  N      = .N,
  Min    = round(min(Score),    3),
  Median = round(median(Score), 3),
  Mean   = round(mean(Score),   3),
  Max    = round(max(Score),    3),
  SD     = round(sd(Score),     3)
), by = label]
print(score_summary)

## ============================================================
## SECTION 8: Score distribution plot
## ============================================================
## Density plot of CAPTAn Score split by label.
##   Red   (#e74c3c) = Negative (label 0)
##   Green (#2ecc71) = Positive (label 1)
## Saved to: captan_dq_score_dist.png
## ============================================================

p1 <- ggplot(matches, aes(x = Score, fill = factor(label))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c("0" = "#e74c3c", "1" = "#2ecc71"),
    labels = c("Negative", "Positive")
  ) +
  theme_minimal() +
  labs(
    title    = "CAPTAn Score Distribution — DQ",
    subtitle = "Best DQ allele score per peptide position",
    x        = "CAPTAn Score",
    y        = "Density",
    fill     = "Label"
  )
print(p1)
ggsave(file.path(BASE_DIR, "captan_dq_score_dist.png"),
       p1, width = 7, height = 5, dpi = 300)

## ============================================================
## SECTION 9: ScoreCore distribution plot
## ============================================================
## Density plot of CAPTAn ScoreCore split by label.
## Saved to: captan_dq_scorecore_dist.png
## ============================================================

p2 <- ggplot(matches, aes(x = ScoreCore, fill = factor(label))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c("0" = "#e74c3c", "1" = "#2ecc71"),
    labels = c("Negative", "Positive")
  ) +
  theme_minimal() +
  labs(
    title    = "CAPTAn ScoreCore Distribution — DQ",
    subtitle = "Binding core score per peptide position",
    x        = "CAPTAn ScoreCore",
    y        = "Density",
    fill     = "Label"
  )
print(p2)
ggsave(file.path(BASE_DIR, "captan_dq_scorecore_dist.png"),
       p2, width = 7, height = 5, dpi = 300)

## ============================================================
## SECTION 10: ROC curve
## ============================================================
## ROC curve with AUC and 95% CI in the title.
## Diagonal dashed line = random classifier baseline.
## coord_equal() enforces a square aspect ratio.
## Saved to: captan_dq_roc.png
## ============================================================

roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

p3 <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#2ecc71", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  labs(
    title = paste0(
      "ROC Curve — CAPTAn DQ",
      "\nAUC = ", round(auc(roc_obj), 3),
      "  95% CI: ", round(ci(roc_obj)[1], 3),
      "–", round(ci(roc_obj)[3], 3)
    ),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_minimal() +
  coord_equal()
print(p3)
ggsave(file.path(BASE_DIR, "captan_dq_roc.png"),
       p3, width = 6, height = 6, dpi = 300)

## ============================================================
## SECTION 11: Precision-recall curve
## ============================================================
## PR curve with AUC in the title.
## Dashed horizontal line = positive class prevalence baseline.
## Saved to: captan_dq_pr_curve.png
## ============================================================

pr_df    <- data.frame(
  recall    = pr_obj$curve[, 1],
  precision = pr_obj$curve[, 2]
)
baseline <- sum(matches$label == 1) / nrow(matches)

p4 <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(color = "#3498db", linewidth = 1) +
  geom_hline(yintercept = baseline,
             linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.7, y = baseline + 0.02,
           label = paste0("Baseline (", round(baseline, 3), ")"),
           color = "grey40", size = 3.5) +
  labs(
    title = paste0(
      "Precision-Recall Curve — CAPTAn DQ",
      "\nAUC = ", round(pr_obj$auc.integral, 3)
    ),
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal()
print(p4)
ggsave(file.path(BASE_DIR, "captan_dq_pr_curve.png"),
       p4, width = 7, height = 5, dpi = 300)

## ============================================================
## SECTION 12: Score vs ScoreCore scatter
## ============================================================
## Scatter plot of Score vs ScoreCore coloured by label.
## Useful for visualising whether core and full scores
## separate labels similarly, and how correlated they are.
## Saved to: captan_dq_score_vs_scorecore.png
## ============================================================

p5 <- ggplot(matches, aes(x = Score, y = ScoreCore,
                           color = factor(label))) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(
    values = c("0" = "#e74c3c", "1" = "#2ecc71"),
    labels = c("Negative", "Positive")
  ) +
  theme_minimal() +
  labs(
    title = "CAPTAn Score vs ScoreCore — DQ",
    x     = "Score",
    y     = "ScoreCore",
    color = "Label"
  )
print(p5)
ggsave(file.path(BASE_DIR, "captan_dq_score_vs_scorecore.png"),
       p5, width = 7, height = 5, dpi = 300)

## ============================================================
## SECTION 13: Save results
## ============================================================
## Writes the full matched evaluation table to CSV.
## ============================================================

fwrite(matches, file.path(BASE_DIR, "captan_dq_evaluation_results.csv"))
cat("\nSaved results\nDONE\n")

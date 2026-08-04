# ============================================================
# Project : CPTAC Multi-Omics MOFA
# Cancer  : Pancreatic Ductal Adenocarcinoma (PDAC)
# Script  : 03_Proteomics_Preprocessing.R
# Purpose : Preprocess PDAC proteomics data and save checkpoint
# ============================================================

rm(list = ls())

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed/PDAC", showWarnings = FALSE, recursive = TRUE)

clinical <- read.csv(
  "data/raw/PDAC/clinical.csv",
  row.names = 1,
  check.names = FALSE
)

proteomics <- read.csv(
  "data/raw/PDAC/proteomics_umich.csv",
  row.names = 1,
  check.names = FALSE
)

proteomics_data <- proteomics[
  !rownames(proteomics) %in% c("Database_ID", "Patient_ID"),
  ,
  drop = FALSE
]

proteomics_data <- as.data.frame(
  lapply(proteomics_data, as.numeric),
  check.names = FALSE
)

rownames(proteomics_data) <- rownames(proteomics)[
  !rownames(proteomics) %in% c("Database_ID", "Patient_ID")
]

missing_rate <- colMeans(is.na(proteomics_data)) * 100

missing_summary <- data.frame(
  Protein = names(missing_rate),
  Missing_Percent = as.numeric(missing_rate),
  row.names = NULL
)

write.csv(
  missing_summary,
  "results/PDAC_proteomics_missingness_summary.csv",
  row.names = FALSE
)

missing_thresholds <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)

proteins_remaining <- sapply(
  missing_thresholds,
  function(threshold) sum(missing_rate <= threshold)
)

missingness_threshold_summary <- data.frame(
  Threshold_Percent = missing_thresholds,
  Proteins_Remaining = as.numeric(proteins_remaining)
)

write.csv(
  missingness_threshold_summary,
  "results/PDAC_proteomics_missingness_threshold_summary.csv",
  row.names = FALSE
)

proteins_keep <- names(missing_rate)[missing_rate <= 50]
proteins_removed <- names(missing_rate)[missing_rate > 50]

proteomics_filtered <- proteomics_data[, proteins_keep, drop = FALSE]

missingness_filter_summary <- data.frame(
  Total_Proteins = ncol(proteomics_data),
  Retained_Proteins = length(proteins_keep),
  Removed_Proteins = length(proteins_removed),
  Retained_Percent = round(100 * length(proteins_keep) / ncol(proteomics_data), 2),
  Removed_Percent = round(100 * length(proteins_removed) / ncol(proteomics_data), 2)
)

write.csv(
  missingness_filter_summary,
  "results/PDAC_proteomics_missingness_filter_summary.csv",
  row.names = FALSE
)

write.csv(
  proteomics_filtered,
  "data/processed/PDAC/proteomics_umich_filtered_50pct.csv"
)

save(
  clinical,
  proteomics,
  proteomics_data,
  missing_rate,
  missing_summary,
  missing_thresholds,
  proteins_remaining,
  missingness_threshold_summary,
  proteins_keep,
  proteins_removed,
  proteomics_filtered,
  missingness_filter_summary,
  file = "results/PDAC_proteomics_preprocessing_checkpoint.RData"
)

cat("\n==============================\n")
cat("Proteomics Preprocessing Done\n")
cat("==============================\n")
cat("Raw proteomics:", dim(proteomics), "\n")
cat("Numeric data:", dim(proteomics_data), "\n")
cat("Filtered proteomics:", dim(proteomics_filtered), "\n")
cat("Retained proteins:", length(proteins_keep), "\n")
cat("Removed proteins:", length(proteins_removed), "\n")

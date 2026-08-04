# ============================================
# Step 02: Clinical Data Summary for PDAC
# Project: CPTAC Multi-Omics MOFA
# ============================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- getwd()

clinical_file <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data/clinical.csv"
output_dir <- file.path(project_dir, "results")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(clinical_file)) {
  stop("File not found: ", clinical_file)
}

clinical <- fread(clinical_file)

cat("=== Clinical data loaded successfully ===\n")
cat("Dimensions:", nrow(clinical), "rows x", ncol(clinical), "columns\n\n")

cat("=== Column names ===\n")
print(colnames(clinical))

stage_col <- "tumor_stage_pathological"


if (!stage_col %in% colnames(clinical)) {
  stop("Stage column not found: ", stage_col)
}

cat("\n=== Stage column detected ===\n")
cat(stage_col, "\n\n")

clinical_summary <- data.frame(
  metric = c("n_samples", "n_variables"),
  value = c(nrow(clinical), ncol(clinical))
)

fwrite(clinical_summary, file.path(output_dir, "PDAC_clinical_summary.csv"))

tumor_stage_summary <- as.data.frame(table(clinical[[stage_col]], useNA = "ifany"))
colnames(tumor_stage_summary) <- c("tumor_stage_pathological", "count")

fwrite(tumor_stage_summary, file.path(output_dir, "PDAC_tumor_stage_summary.csv"))

cat("=== Tumor stage distribution ===\n")
print(tumor_stage_summary)

cat("\n=== Output files written ===\n")
cat(file.path(output_dir, "PDAC_clinical_summary.csv"), "\n")
cat(file.path(output_dir, "PDAC_tumor_stage_summary.csv"), "\n")

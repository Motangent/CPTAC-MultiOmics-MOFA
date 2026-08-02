# ============================================================
# Project : CPTAC Multi-Omics MOFA
# Cancer  : Pancreatic Ductal Adenocarcinoma (PDAC)
# Script  : 01_Load_PDAC_Data.R
# Purpose : Load exported CPTAC PDAC datasets from CSV files
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Load Clinical Data
# ------------------------------------------------------------

clinical <- read.csv(
  "data/raw/PDAC/clinical.csv",
  row.names = 1,
  check.names = FALSE
)

# ------------------------------------------------------------
# Quick inspection
# ------------------------------------------------------------

str(clinical)

table(clinical$tumor_stage_pathological)

# ------------------------------------------------------------
# Dataset dimensions
# ------------------------------------------------------------

proteomics <- read.csv(
  "data/raw/PDAC/proteomics_umich.csv",
  row.names = 1,
  check.names = FALSE
)

transcriptomics <- read.csv(
  "data/raw/PDAC/transcriptomics_broad.csv",
  row.names = 1,
  check.names = FALSE
)

phosphoproteomics <- read.csv(
  "data/raw/PDAC/phosphoproteomics_umich.csv",
  row.names = 1,
  check.names = FALSE
)

cat("\n==============================\n")
cat("PDAC Dataset Summary\n")
cat("==============================\n")

cat("Clinical:", dim(clinical), "\n")
cat("Proteomics:", dim(proteomics), "\n")
cat("Transcriptomics:", dim(transcriptomics), "\n")
cat("Phosphoproteomics:", dim(phosphoproteomics), "\n")


# ============================================================
# Script  : 05_Transcriptomics_Preprocessing.R
# Purpose : Transcriptomics preprocessing outputs (NO rerun)
# Input   : data/processed/PDAC_Step2_Workspace.RData
# Output  : data/processed/PDAC/transcriptomics_bcm_overlap_140x59286.csv
#           data/processed/PDAC/transcriptomics_bcm_hvg_5000_140x5000.csv
#           results/PDAC_transcriptomics_preprocessing_checkpoint.RData
# ============================================================

rm(list = ls())

load("data/processed/PDAC_Step2_Workspace.RData")

stopifnot(exists("trans_overlap"), exists("trans_hvg"))
stopifnot(all(dim(trans_overlap) == c(140, 59286)))
stopifnot(all(dim(trans_hvg) == c(140, 5000)))

dir.create("data/processed/PDAC", recursive = TRUE, showWarnings = FALSE)

write.csv(
  trans_overlap,
  "data/processed/PDAC/transcriptomics_bcm_overlap_140x59286.csv",
  row.names = TRUE
)

write.csv(
  trans_hvg,
  "data/processed/PDAC/transcriptomics_bcm_hvg_5000_140x5000.csv",
  row.names = TRUE
)

save(
  trans_overlap,
  trans_hvg,
  file = "results/PDAC_transcriptomics_preprocessing_checkpoint.RData"
)

cat("\nTranscriptomics preprocessing checkpoint created (no rerun).\n")

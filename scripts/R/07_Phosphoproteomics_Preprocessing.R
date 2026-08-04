# Step 07: Phosphoproteomics preprocessing outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("phospho_filtered"))
stopifnot(is.matrix(phospho_filtered) || is.data.frame(phospho_filtered))

dir.create("data/processed/PDAC", recursive = TRUE, showWarnings = FALSE)

write.csv(
  phospho_filtered,
  "data/processed/PDAC/phosphoproteomics_filtered_50pct.csv",
  row.names = TRUE
)

save(
  phospho_filtered,
  file = "results/PDAC_phospho_preprocessing_checkpoint.RData"
)

cat("\nPhosphoproteomics preprocessing checkpoint created (no rerun).\n")

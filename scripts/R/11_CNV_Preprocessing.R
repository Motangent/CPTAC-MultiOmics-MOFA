# Step 11: CNV preprocessing outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("cnv_bcm_final"))

dir.create("data/processed/PDAC", recursive = TRUE, showWarnings = FALSE)

write.csv(
  cnv_bcm_final,
  "data/processed/PDAC/cnv_bcm_final.csv",
  row.names = TRUE
)

save(
  cnv_bcm_final,
  file = "results/PDAC_cnv_preprocessing_checkpoint.RData"
)

cat("\nCNV preprocessing checkpoint created (no rerun).\n")

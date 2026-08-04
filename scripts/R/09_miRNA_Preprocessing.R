# Step 09: miRNA preprocessing outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("mirna_bcm_final"))

dir.create("data/processed/PDAC", recursive = TRUE, showWarnings = FALSE)

write.csv(
  mirna_bcm_final,
  "data/processed/PDAC/mirna_bcm_final.csv",
  row.names = TRUE
)

save(
  mirna_bcm_final,
  file = "results/PDAC_mirna_preprocessing_checkpoint.RData"
)

cat("\nmiRNA preprocessing checkpoint created (no rerun).\n")

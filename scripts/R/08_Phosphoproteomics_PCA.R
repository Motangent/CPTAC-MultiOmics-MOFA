# Step 08: Phosphoproteomics PCA outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("pca_phospho"))
stopifnot(inherits(pca_phospho, "prcomp"))

variance_explained <- data.frame(
  PC = paste0("PC", seq_along(pca_phospho$sdev)),
  Standard_Deviation = pca_phospho$sdev,
  Variance_Explained = pca_phospho$sdev^2 / sum(pca_phospho$sdev^2),
  Cumulative_Variance = cumsum(pca_phospho$sdev^2 / sum(pca_phospho$sdev^2))
)

dir.create("results/PCA", recursive = TRUE, showWarnings = FALSE)

write.csv(
  variance_explained,
  "results/PCA/PDAC_Phosphoproteomics_PCA_variance_explained.csv",
  row.names = FALSE
)

save(
  pca_phospho,
  variance_explained,
  file = "results/PCA/PDAC_Phosphoproteomics_PCA_Checkpoint.RData"
)

cat("\nPhosphoproteomics PCA checkpoint created from existing object (no rerun).\n")

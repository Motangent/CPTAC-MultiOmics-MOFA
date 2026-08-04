# Step 10: miRNA PCA outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("pca_mirna"))
stopifnot(inherits(pca_mirna, "prcomp"))

variance_explained <- data.frame(
  PC = paste0("PC", seq_along(pca_mirna$sdev)),
  Standard_Deviation = pca_mirna$sdev,
  Variance_Explained = pca_mirna$sdev^2 / sum(pca_mirna$sdev^2),
  Cumulative_Variance = cumsum(pca_mirna$sdev^2 / sum(pca_mirna$sdev^2))
)

dir.create("results/PCA", recursive = TRUE, showWarnings = FALSE)

write.csv(
  variance_explained,
  "results/PCA/PDAC_miRNA_PCA_variance_explained.csv",
  row.names = FALSE
)

save(
  pca_mirna,
  variance_explained,
  file = "results/PCA/PDAC_miRNA_PCA_Checkpoint.RData"
)

cat("\nmiRNA PCA checkpoint created from existing object (no rerun).\n")

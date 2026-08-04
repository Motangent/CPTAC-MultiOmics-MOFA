# Step 12: CNV PCA outputs (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

stopifnot(exists("pca_cnv"))
stopifnot(inherits(pca_cnv, "prcomp"))

variance_explained <- data.frame(
  PC = paste0("PC", seq_along(pca_cnv$sdev)),
  Standard_Deviation = pca_cnv$sdev,
  Variance_Explained = pca_cnv$sdev^2 / sum(pca_cnv$sdev^2),
  Cumulative_Variance = cumsum(pca_cnv$sdev^2 / sum(pca_cnv$sdev^2))
)

dir.create("results/PCA", recursive = TRUE, showWarnings = FALSE)

write.csv(
  variance_explained,
  "results/PCA/PDAC_CNV_PCA_variance_explained.csv",
  row.names = FALSE
)

save(
  pca_cnv,
  variance_explained,
  file = "results/PCA/PDAC_CNV_PCA_Checkpoint.RData"
)

cat("\nCNV PCA checkpoint created from existing object (no rerun).\n")

# Script 06: Transcriptomics PCA output from existing checkpoint

rm(list = ls())

load("data/processed/PDAC_Step2_Workspace.RData")

stopifnot(exists("pca_trans"))
stopifnot(inherits(pca_trans, "prcomp"))

variance_explained <- data.frame(
  PC = paste0("PC", seq_along(pca_trans$sdev)),
  Standard_Deviation = pca_trans$sdev,
  Variance_Explained = pca_trans$sdev^2 / sum(pca_trans$sdev^2),
  Cumulative_Variance = cumsum(pca_trans$sdev^2 / sum(pca_trans$sdev^2))
)

dir.create("results/PCA", recursive = TRUE, showWarnings = FALSE)

write.csv(
  variance_explained,
  "results/PCA/PDAC_Transcriptomics_PCA_variance_explained.csv",
  row.names = FALSE
)

save(
  pca_trans,
  variance_explained,
  file = "results/PCA/PDAC_Transcriptomics_PCA_Checkpoint.RData"
)

cat("\nTranscriptomics PCA checkpoint created from existing object (no rerun).\n")

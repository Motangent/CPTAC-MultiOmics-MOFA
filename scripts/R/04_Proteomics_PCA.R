# ============================================================
# Project : CPTAC Multi-Omics MOFA
# Cancer  : Pancreatic Ductal Adenocarcinoma (PDAC)
# Script  : 04_Proteomics_PCA.R
# Purpose : Perform complete-case PCA on proteomics data
# ============================================================

rm(list = ls())

# Load data and filter
load("results/PDAC_proteomics_preprocessing_checkpoint.RData")

# Ensure complete-case for PCA (as per Lab Book, page 23)
proteomics_complete <- proteomics_filtered[, colSums(is.na(proteomics_filtered)) == 0, drop = FALSE]

# Execute PCA
pca_res <- prcomp(proteomics_complete, scale. = TRUE)

# Calculate variance explained
pca_var <- pca_res$sdev^2
pca_var_exp <- round(pca_var / sum(pca_var) * 100, 2)

# Prepare variance summary
variance_summary <- data.frame(
  PC = paste0("PC", 1:length(pca_var_exp)),
  Variance_Explained = pca_var_exp
)

write.csv(
  variance_summary,
  "results/PDAC_proteomics_PCA_variance_explained.csv",
  row.names = FALSE
)

# Save PCA checkpoint
save(
  pca_res,
  proteomics_complete,
  variance_summary,
  file = "results/PDAC_Proteomics_PCA_Final_Checkpoint.RData"
)

cat("\n==============================\n")
cat("Proteomics PCA Done\n")
cat("==============================\n")
cat("Complete proteins used:", ncol(proteomics_complete), "\n")
cat("PC1 Variance:", pca_var_exp[1], "%\n")
cat("PC2 Variance:", pca_var_exp[2], "%\n")

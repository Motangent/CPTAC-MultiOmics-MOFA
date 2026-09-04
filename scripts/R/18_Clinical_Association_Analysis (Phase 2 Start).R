# ==============================================================================
# CPTAC PDAC: Complete Clinical Associations Analysis (Zero-to-Hero)
# ==============================================================================

# 1. Clear Environment & Close Graphics Devices
rm(list = ls())
while (!is.null(dev.list())) {
  dev.off()
}

# 2. Set Working Directory & Load Required Package
setwd("C:/Users/darvi/Desktop/MOFA/CPTAC-MultiOmics-MOFA-Project/CPTAC-MultiOmics-MOFA")
library(pheatmap)

# 3. Load & Reshape MOFA Factors (105 samples × 15 factors)
factors_long <- read.csv("MOFA_factors_all.csv", check.names = FALSE, stringsAsFactors = FALSE)
mofa_pdac_factors <- reshape(
  factors_long[, c("sample", "factor", "value")],
  idvar = "sample",
  timevar = "factor",
  direction = "wide"
)

rownames(mofa_pdac_factors) <- mofa_pdac_factors$sample
mofa_pdac_factors$sample <- NULL
colnames(mofa_pdac_factors) <- sub("^value\\.", "", colnames(mofa_pdac_factors))

fac_mat <- as.data.frame(lapply(mofa_pdac_factors, as.numeric))
rownames(fac_mat) <- rownames(mofa_pdac_factors)
factors_list <- grep("^Factor[0-9]+$", colnames(fac_mat), value = TRUE)

# 4. Load Clinical Data & Align Samples
load("data/processed/PDAC/PDAC_MOFA_Inputs_105_Samples.RData")

common_ids <- intersect(rownames(fac_mat), clinical_data$Patient_ID)
fac_mat <- fac_mat[common_ids, factors_list, drop = FALSE]
cli_sub <- clinical_data[match(common_ids, clinical_data$Patient_ID), , drop = FALSE]

stopifnot(
  identical(rownames(fac_mat), cli_sub$Patient_ID),
  nrow(fac_mat) == 105,
  ncol(fac_mat) == 15
)

# 5. Statistical Tests for Categorical Variables (Kruskal-Wallis)
clin_vars <- c("tumor_stage_pathological", "histologic_grade", "sex", "tobacco_smoking_history")

p_mat <- matrix(
  NA_real_,
  nrow = length(clin_vars),
  ncol = length(factors_list),
  dimnames = list(clin_vars, factors_list)
)

for (v in clin_vars) {
  for (f in factors_list) {
    x <- cli_sub[[v]]
    y <- fac_mat[[f]]
    valid <- !is.na(x) & x != "" & is.finite(y)
    if (sum(valid) >= 3 && length(unique(x[valid])) > 1 && length(unique(y[valid])) > 1) {
      p_mat[v, f] <- kruskal.test(y[valid] ~ factor(x[valid]))$p.value
    }
  }
}

# 6. Statistical Test for Continuous Variable (Age - Spearman Correlation)
age_p <- sapply(factors_list, function(f) {
  age <- suppressWarnings(as.numeric(as.character(cli_sub$age)))
  y <- fac_mat[[f]]
  valid <- is.finite(age) & is.finite(y)
  if (sum(valid) >= 3 && length(unique(age[valid])) > 1 && length(unique(y[valid])) > 1) {
    cor.test(age[valid], y[valid], method = "spearman", exact = FALSE)$p.value
  } else {
    NA_real_
  }
})

# 7. Combine Results & Multi-Testing Correction (BH-FDR)
p_mat_full <- rbind(p_mat, age = age_p)

p_adj <- matrix(
  NA_real_,
  nrow = nrow(p_mat_full),
  ncol = ncol(p_mat_full),
  dimnames = dimnames(p_mat_full)
)

valid_p <- is.finite(p_mat_full)
p_adj[valid_p] <- p.adjust(p_mat_full[valid_p], method = "BH")

# 8. Render Heatmap Directly in Plots Pane (No Automatic File Saving)
pheatmap(
  mat = p_adj,
  display_numbers = TRUE,
  number_format = "%.3f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "CPTAC PDAC Clinical Associations: BH-FDR-adjusted p-values",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  na_col = "grey90",
  border_color = "grey80",
  fontsize = 10,
  fontsize_number = 8
)

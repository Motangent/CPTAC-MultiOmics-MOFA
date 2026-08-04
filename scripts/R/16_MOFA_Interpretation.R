# ============================================================
# 16_MOFA_Interpretation.R
# CPTAC PDAC Multi-Omics Integration with MOFA2
# Purpose: Extract and register MOFA interpretation outputs
# No model retraining
# ============================================================

rm(list = ls())

library(MOFA2)

load("results/MOFA/PDAC_MOFA_Training_Checkpoint.RData")

stopifnot(
  exists("mofa_trained"),
  inherits(mofa_trained, "MOFA"),
  mofa_trained@status == "trained"
)

dir.create("results/MOFA/Interpretation", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Model dimensions
# -----------------------------

model_dimensions <- data.frame(
  Dimension = c("Views", "Groups", "Samples", "Factors"),
  Value = c(
    mofa_trained@dimensions$M,
    mofa_trained@dimensions$G,
    mofa_trained@dimensions$N[[1]],
    mofa_trained@dimensions$K
  )
)

write.csv(
  model_dimensions,
  "results/MOFA/Interpretation/PDAC_MOFA_Model_Dimensions.csv",
  row.names = FALSE
)

# -----------------------------
# 2. Variance explained
# -----------------------------

variance_explained <- get_variance_explained(mofa_trained)

total_variance <- data.frame(
  View = names(variance_explained$r2_total[[1]]),
  Variance_Explained_Percent = as.numeric(
    variance_explained$r2_total[[1]]
  )
)

write.csv(
  total_variance,
  "results/MOFA/Interpretation/PDAC_MOFA_Total_Variance_Explained.csv",
  row.names = FALSE
)

per_factor_variance <- as.data.frame(
  variance_explained$r2_per_factor[[1]]
)

per_factor_variance$Factor <- rownames(per_factor_variance)
per_factor_variance <- per_factor_variance[, c(
  "Factor",
  setdiff(names(per_factor_variance), "Factor")
)]

write.csv(
  per_factor_variance,
  "results/MOFA/Interpretation/PDAC_MOFA_Variance_Explained_Per_Factor.csv",
  row.names = FALSE
)

# -----------------------------
# 3. Factor scores
# -----------------------------

factor_scores <- get_factors(
  mofa_trained,
  factors = "all",
  groups = "all",
  as.data.frame = TRUE
)

write.csv(
  factor_scores,
  "results/MOFA/Interpretation/PDAC_MOFA_Factor_Scores.csv",
  row.names = FALSE
)

# -----------------------------
# 4. Feature weights
# -----------------------------

feature_weights <- get_weights(
  mofa_trained,
  views = "all",
  factors = "all",
  scale = TRUE,
  as.data.frame = TRUE
)

write.csv(
  feature_weights,
  "results/MOFA/Interpretation/PDAC_MOFA_Feature_Weights.csv",
  row.names = FALSE
)

# -----------------------------
# 5. Standard plots
# -----------------------------

pdf(
  "results/MOFA/Interpretation/PDAC_MOFA_Variance_Explained_Plot.pdf",
  width = 9,
  height = 6
)

plot_variance_explained(
  mofa_trained,
  x = "view",
  y = "factor",
  plot_total = TRUE
)


dev.off()

pdf(
  "results/MOFA/Interpretation/PDAC_MOFA_Variance_Explained_Per_Factor_Plot.pdf",
  width = 10,
  height = 7
)

plot_variance_explained(
  mofa_trained,
  x = "factor",
  y = "view"
)

dev.off()

# -----------------------------
# 6. Final interpretation checkpoint
# -----------------------------

save(
  mofa_trained,
  variance_explained,
  total_variance,
  per_factor_variance,
  factor_scores,
  feature_weights,
  model_dimensions,
  file = "results/MOFA/Interpretation/PDAC_MOFA_Interpretation_Checkpoint.RData"
)

cat(
  "\nMOFA interpretation completed successfully.\n",
  "Model: 5 views | 105 samples | 15 factors\n",
  "Outputs saved in: results/MOFA/Interpretation/\n",
  sep = ""
)


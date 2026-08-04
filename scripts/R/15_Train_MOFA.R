# Step 15: MOFA trained-model checkpoint (NO retraining)

rm(list = ls())

load("PDAC_Phase6_PostMOFA_SAFE.RData")

stopifnot(exists("mofa_trained"))

dir.create("results/MOFA", recursive = TRUE, showWarnings = FALSE)

save(
  mofa_trained,
  file = "results/MOFA/PDAC_MOFA_Training_Checkpoint.RData"
)

cat(
  "\nMOFA trained-model checkpoint registered successfully (no retraining).\n"
)

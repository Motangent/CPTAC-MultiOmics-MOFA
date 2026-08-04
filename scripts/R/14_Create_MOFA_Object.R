# Step 14: Create MOFA object (NO rerun of data processing)

rm(list = ls())

# Load all five final omics matrices
load("results/PDAC_common_five_omics_checkpoint.RData")

# Load clinical data for sample metadata
# We assume clinical object is available from earlier scripts' environment or needs to be loaded separately.
# For simplicity and speed, we will load the main five-omics checkpoint again to ensure access to prot_aligned for clinical info
load("results/PDAC_aligned_five_omics_checkpoint.RData")
clinical_data <- clinical_aligned[Patient_ID %in% common_samples]

# --- Create MOFA Data Structure ---
# Data must be a list of matrices where samples are rows and features are columns (MOFA2 requires samples in rows)

data_list <- list(
  Proteomics = prot_final,
  Transcriptomics = trans_final,
  Phosphoproteomics = phospho_final,
  miRNA = mirna_final,
  CNV = cnv_final
)

# --- Create MOFA Model Options (Simplified) ---

# This step does NOT train the model, it just creates the object for training.
# MOFA2 package is needed here, though we avoid full training.

# We will export the aligned data and the clinical data for the next step (15_Train_MOFA.R) to use.
# Since running MOFA installation and object creation might add overhead, we simply save the final inputs needed.

save(
  data_list,
  clinical_data,
  file = "data/processed/PDAC/PDAC_MOFA_Inputs_105_Samples.RData"
)

# This is the checkpoint for the training step
save(
  data_list,
  clinical_data,
  file = "results/PDAC_MOFA_initial_checkpoint.RData"
)

cat(
  "\nMOFA data input structure created for",
  length(data_list),
  "views with",
  nrow(data_list[[1]]),
  "samples.\n"
)

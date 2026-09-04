# ==============================================================================
# Script 19: Survival Data Preparation and Quality Control
# Project: CPTAC-MultiOmics-MOFA-Project (PDAC)
# Purpose: Prepare and validate overall-survival data aligned with MOFA factors
# ==============================================================================

# 1. Verification of Required Workspace Objects
stopifnot(
  exists("clinical_data"),
  exists("mofa_pdac_factors"),
  exists("common_ids")
)

# 2. Define the 15 MOFA factor names
factor_cols <- paste0("Factor", 1:15)

# 3. Extract and Standardize Clinical Survival Variables
clinical_survival <- data.frame(
  Sample_ID = as.character(common_ids),
  time      = as.numeric(clinical_data$`Overall survival, days`),
  status    = as.numeric(clinical_data$`Survival status (1, dead; 0, alive)`),
  stringsAsFactors = FALSE
)

# 4. Prepare MOFA Factor Data Frame
mofa_factors_df <- as.data.frame(mofa_pdac_factors)
mofa_factors_df$Sample_ID <- rownames(mofa_pdac_factors)
mofa_factors_df <- mofa_factors_df[, c("Sample_ID", factor_cols), drop = FALSE]

# 5. Validate the Clinical Survival Table Structure
stopifnot(
  nrow(clinical_survival) == 105,
  identical(colnames(clinical_survival), c("Sample_ID", "time", "status")),
  !anyDuplicated(clinical_survival$Sample_ID),
  all(!is.na(clinical_survival$Sample_ID))
)

# 6. Merge Clinical Survival Data with MOFA Factors
survival_merged <- merge(
  clinical_survival,
  mofa_factors_df,
  by = "Sample_ID",
  all = FALSE,
  sort = FALSE
)

# 7. Filter for Valid Survival Analysis Cases
valid_survival <- !is.na(survival_merged$time) &
  survival_merged$time > 0 &
  !is.na(survival_merged$status) &
  survival_merged$status %in% c(0, 1) &
  complete.cases(survival_merged[, factor_cols, drop = FALSE])

survival_ready <- survival_merged[valid_survival, c("Sample_ID", "time", "status", factor_cols), drop = FALSE]

# 8. Identify Excluded Samples and Exclusion Reasons
all_ids     <- as.character(clinical_survival$Sample_ID)
kept_ids    <- as.character(survival_ready$Sample_ID)
removed_ids <- setdiff(all_ids, kept_ids)

removed_report <- clinical_survival[
  match(removed_ids, clinical_survival$Sample_ID),
  c("Sample_ID", "time", "status"),
  drop = FALSE
]

removed_reasons <- data.frame(
  Sample_ID         = removed_report$Sample_ID,
  missing_time      = is.na(removed_report$time),
  non_positive_time = (!is.na(removed_report$time) & removed_report$time <= 0),
  missing_status    = is.na(removed_report$status),
  invalid_status    = (!is.na(removed_report$status) & !removed_report$status %in% c(0, 1)),
  stringsAsFactors  = FALSE
)

# 9. Final Quality Control Assertions
stopifnot(
  nrow(clinical_survival) == 105,
  nrow(survival_ready) == 102,
  ncol(survival_ready) == 18,
  length(removed_ids) == 3,
  identical(colnames(survival_ready), c("Sample_ID", "time", "status", factor_cols)),
  !anyDuplicated(survival_ready$Sample_ID),
  all(!is.na(survival_ready$Sample_ID)),
  all(!is.na(survival_ready$time)),
  all(survival_ready$time > 0),
  all(!is.na(survival_ready$status)),
  all(survival_ready$status %in% c(0, 1)),
  sum(is.na(survival_ready)) == 0,
  all(colSums(is.na(survival_ready[, factor_cols, drop = FALSE])) == 0)
)

# 10. Summary Output
cat("=== Step 19: Survival Data Preparation Completed ===\n")
cat("Clinical samples:", nrow(clinical_survival), "\n")
cat("Retained samples:", nrow(survival_ready), "\n")
cat("Removed samples:", length(removed_ids), "\n")
cat("Deaths (status = 1):", sum(survival_ready$status == 1), "\n")
cat("Censored (status = 0):", sum(survival_ready$status == 0), "\n")
cat("MOFA factors:", length(factor_cols), "\n")
cat("Missing values in survival_ready:", sum(is.na(survival_ready)), "\n\n")

cat("Excluded samples report:\n")
print(removed_reasons, row.names = FALSE)

cat("\nSurvival time summary (days):\n")
print(summary(survival_ready$time))

cat("\nFinal quality control: PASSED\n")

# ==============================================================================
# Script 22: Multivariable Cox Regression Adjusted for Clinical Variables
# ==============================================================================

library(survival)

# 1. Load datasets
rds_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC-MultiOmics-MOFA-Project/data/processed/PDAC/survival_ready.rds"
clinical_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data/clinical.csv"

stopifnot(file.exists(rds_path))
stopifnot(file.exists(clinical_path))

survival_ready <- readRDS(rds_path)

clinical <- read.csv(
  clinical_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 2. Select clinical variables
clinical_sub <- clinical[, c(
  "Patient_ID",
  "age",
  "sex",
  "tumor_stage_pathological"
)]

# 3. Remove duplicate patient records
clinical_sub <- clinical_sub[!duplicated(clinical_sub$Patient_ID), ]

# 4. Merge MOFA survival data with clinical data
merged <- merge(
  survival_ready,
  clinical_sub,
  by.x = "Sample_ID",
  by.y = "Patient_ID",
  all.x = TRUE,
  sort = FALSE
)

cat("Survival-ready samples:", nrow(survival_ready), "\n")
cat("Merged samples:", nrow(merged), "\n")
cat("Matched age:", sum(!is.na(merged$age)), "\n")
cat("Matched sex:", sum(!is.na(merged$sex)), "\n")
cat("Matched stage:", sum(!is.na(merged$tumor_stage_pathological)), "\n\n")

# 5. Keep complete cases for the multivariable model
model_data <- merged[
  complete.cases(
    merged[, c(
      "time",
      "status",
      "Factor2",
      "Factor9",
      "age",
      "sex",
      "tumor_stage_pathological"
    )]
  ),
]

cat("Complete cases used in model:", nrow(model_data), "\n\n")

# 6. Prepare model variables
model_data$sex <- factor(model_data$sex)

model_data$tumor_stage_pathological <- factor(
  model_data$tumor_stage_pathological
)

# 7. Fit multivariable Cox regression
cox_multi <- coxph(
  Surv(time, status) ~
    Factor2 +
    Factor9 +
    age +
    sex +
    tumor_stage_pathological,
  data = model_data
)

# 8. Display model summary
cat("=== Multivariable Cox Regression Summary ===\n\n")
print(summary(cox_multi))

# 9. Create formatted results table
cox_summary <- summary(cox_multi)

results_multi <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  CI_Lower = cox_summary$conf.int[, "lower .95"],
  CI_Upper = cox_summary$conf.int[, "upper .95"],
  P_Value = cox_summary$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)

cat("\n=== Formatted Hazard Ratios ===\n")
print(results_multi, digits = 4)

# 10. Check proportional-hazards assumption
cat("\n=== Proportional-Hazards Test ===\n\n")
print(cox.zph(cox_multi))

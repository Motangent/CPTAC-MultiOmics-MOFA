# ==============================================================================
# Script 20: Univariate Cox Proportional Hazards Regression
# Project: CPTAC-MultiOmics-MOFA-Project (PDAC)
# Purpose: Evaluate prognostic value of all 15 MOFA factors on Overall Survival
# ==============================================================================

# 1. Load Survival Package
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
library(survival)

# 2. Load the Prepared Survival Dataset
rds_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC-MultiOmics-MOFA-Project/data/processed/PDAC/survival_ready.rds"

if (!file.exists(rds_path)) {
  stop("Survival file not found at: ", rds_path)
}

survival_ready <- readRDS(rds_path)

# Verify loaded data integrity
stopifnot(
  nrow(survival_ready) == 102,
  ncol(survival_ready) == 18,
  sum(is.na(survival_ready)) == 0
)

# 3. Univariate Cox Regression across all 15 MOFA Factors
factor_cols <- paste0("Factor", 1:15)

cox_results_list <- lapply(factor_cols, function(f) {
  formula_str <- as.formula(paste("Surv(time, status) ~", f))
  fit <- coxph(formula_str, data = survival_ready)
  fit_summary <- summary(fit)
  
  coef_val <- fit_summary$coefficients[f, "coef"]
  hr_val   <- fit_summary$coefficients[f, "exp(coef)"]
  se_val   <- fit_summary$coefficients[f, "se(coef)"]
  p_val    <- fit_summary$coefficients[f, "Pr(>|z|)"]
  ci_lower <- fit_summary$conf.int[f, "lower .95"]
  ci_upper <- fit_summary$conf.int[f, "upper .95"]
  
  data.frame(
    Factor   = f,
    Coef     = round(coef_val, 4),
    HR       = round(hr_val, 4),
    CI_95_LL = round(ci_lower, 4),
    CI_95_UL = round(ci_upper, 4),
    P_Value  = p_val,
    stringsAsFactors = FALSE
  )
})

cox_summary <- do.call(rbind, cox_results_list)

# 4. Multiple Testing Correction (FDR / Benjamini-Hochberg)
cox_summary$FDR <- p.adjust(cox_summary$P_Value, method = "BH")

# Sort by raw P-value
cox_summary_sorted <- cox_summary[order(cox_summary$P_Value), ]

# Prepare display table
cox_display <- cox_summary_sorted
cox_display$P_Value <- formatC(cox_display$P_Value, format = "e", digits = 3)
cox_display$FDR     <- formatC(cox_display$FDR, format = "e", digits = 3)

# 5. Display Results
cat("=== Step 20: Univariate Cox Regression Results ===\n\n")
print(cox_display, row.names = FALSE)

cat("\nNominally significant factors (P < 0.05):\n")
sig_p <- cox_summary_sorted$Factor[cox_summary_sorted$P_Value < 0.05]
if (length(sig_p) > 0) {
  print(sig_p)
} else {
  cat("None\n")
}

cat("\nFDR significant factors (FDR < 0.05):\n")
sig_fdr <- cox_summary_sorted$Factor[cox_summary_sorted$FDR < 0.05]
if (length(sig_fdr) > 0) {
  print(sig_fdr)
} else {
  cat("None\n")
}

# ==============================================================================
# Script 21: Proportional Hazards Assumption & Kaplan-Meier Survival Curves
# Project: CPTAC-MultiOmics-MOFA-Project (PDAC)
# ==============================================================================

# 1. Load Required Packages
library(survival)
library(survminer)

# 2. Load Survival Dataset
rds_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC-MultiOmics-MOFA-Project/data/processed/PDAC/survival_ready.rds"
stopifnot(file.exists(rds_path))
survival_ready <- readRDS(rds_path)

stopifnot(
  nrow(survival_ready) == 102,
  all(c("time", "status", "Factor2", "Factor9") %in% colnames(survival_ready))
)

# 3. Test Proportional Hazards Assumption (Schoenfeld Residuals)
cox_f2 <- coxph(Surv(time, status) ~ Factor2, data = survival_ready)
cox_f9 <- coxph(Surv(time, status) ~ Factor9, data = survival_ready)

zph_f2 <- cox.zph(cox_f2)
zph_f9 <- cox.zph(cox_f9)

cat("=== Proportional Hazards Assumption Test (Schoenfeld Residuals) ===\n")
cat("\n--- Factor 2 PH Test ---\n")
print(zph_f2)

cat("\n--- Factor 9 PH Test ---\n")
print(zph_f9)

# 4. Stratify Patients by Median Split (High vs. Low)
survival_ready$Factor2_Group <- factor(
  ifelse(survival_ready$Factor2 >= median(survival_ready$Factor2), "High", "Low"),
  levels = c("Low", "High")
)

survival_ready$Factor9_Group <- factor(
  ifelse(survival_ready$Factor9 >= median(survival_ready$Factor9), "High", "Low"),
  levels = c("Low", "High")
)

# 5. Fit Kaplan-Meier Curves
km_fit_f2 <- survfit(Surv(time, status) ~ Factor2_Group, data = survival_ready)
km_fit_f9 <- survfit(Surv(time, status) ~ Factor9_Group, data = survival_ready)

cat("\n=== Kaplan-Meier Fit Summary (Median Survival in Days) ===\n")
cat("\n--- Factor 2 Groups ---\n")
print(km_fit_f2)

cat("\n--- Factor 9 Groups ---\n")
print(km_fit_f9)

# 6. Log-Rank Tests
lr_f2 <- survdiff(Surv(time, status) ~ Factor2_Group, data = survival_ready)
lr_f9 <- survdiff(Surv(time, status) ~ Factor9_Group, data = survival_ready)

p_lr_f2 <- 1 - pchisq(lr_f2$chisq, length(lr_f2$n) - 1)
p_lr_f9 <- 1 - pchisq(lr_f9$chisq, length(lr_f9$n) - 1)

cat(sprintf("\nLog-Rank Test P-value (Factor 2): %.4e\n", p_lr_f2))
cat(sprintf("Log-Rank Test P-value (Factor 9): %.4e\n", p_lr_f9))

# 7. Generate and Display Kaplan-Meier Plots
p1 <- ggsurvplot(
  km_fit_f2,
  data = survival_ready,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Overall Survival by MOFA Factor 2 (Median Split)",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  legend.title = "Factor 2",
  legend.labs = c("Low", "High"),
  ggtheme = theme_classic(base_size = 12)
)

p2 <- ggsurvplot(
  km_fit_f9,
  data = survival_ready,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#E64B35", "#00A087"),
  title = "Overall Survival by MOFA Factor 9 (Median Split)",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  legend.title = "Factor 9",
  legend.labs = c("Low", "High"),
  ggtheme = theme_classic(base_size = 12)
)

print(p1)
print(p2)

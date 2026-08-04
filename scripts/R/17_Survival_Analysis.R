# ============================================================
# Project : CPTAC Multi-Omics MOFA
# Cancer  : Pancreatic Ductal Adenocarcinoma (PDAC)
# Script  : 02_Survival_Analysis.R
# Purpose : Survival analysis for MOFA Factor1 and Factor3
# ============================================================

rm(list = ls())

library(survival)
library(survminer)

raw_data_dir <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data"
survival_dir <- file.path("results", "Survival")

if (!dir.exists(survival_dir)) {
  dir.create(survival_dir, recursive = TRUE)
}

clinical <- read.csv(
  file.path(raw_data_dir, "clinical.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

factor_df <- read.csv(
  "PDAC_MOFA_results/factor_values_wide.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

factor_df <- factor_df[, c("sample", "Factor1", "Factor3")]
clinical$sample <- clinical$Patient_ID

merged_survival <- merge(factor_df, clinical, by = "sample")

cox_data <- merged_survival[, c(
  "sample",
  "Factor1",
  "Factor3",
  "Overall survival, days",
  "Survival status (1, dead; 0, alive)"
)]

missing_survival <- cox_data[
  !complete.cases(cox_data),
  c(
    "sample",
    "Overall survival, days",
    "Survival status (1, dead; 0, alive)"
  )
]

write.csv(
  missing_survival,
  file.path(
    survival_dir,
    "Excluded_Samples_Missing_Survival_Data.csv"
  ),
  row.names = FALSE
)

cox_data <- cox_data[complete.cases(cox_data), ]

cox_f1 <- coxph(
  Surv(
    `Overall survival, days`,
    `Survival status (1, dead; 0, alive)`
  ) ~ Factor1,
  data = cox_data
)

cox_f3 <- coxph(
  Surv(
    `Overall survival, days`,
    `Survival status (1, dead; 0, alive)`
  ) ~ Factor3,
  data = cox_data
)

cox_model <- coxph(
  Surv(
    `Overall survival, days`,
    `Survival status (1, dead; 0, alive)`
  ) ~ Factor1 + Factor3,
  data = cox_data
)

cox_models <- list(
  Factor1_univariable = cox_f1,
  Factor3_univariable = cox_f3,
  Factor1_Factor3_multivariable = cox_model
)

cox_results <- do.call(
  rbind,
  lapply(names(cox_models), function(model_name) {
    fit <- cox_models[[model_name]]
    fit_summary <- summary(fit)
    
    data.frame(
      model = model_name,
      variable = rownames(fit_summary$coefficients),
      n = fit$n,
      events = fit$nevent,
      beta = fit_summary$coefficients[, "coef"],
      hazard_ratio = fit_summary$coefficients[, "exp(coef)"],
      ci_lower_95 = fit_summary$conf.int[, "lower .95"],
      ci_upper_95 = fit_summary$conf.int[, "upper .95"],
      p_value = fit_summary$coefficients[, "Pr(>|z|)"],
      row.names = NULL,
      check.names = FALSE
    )
  })
)

ph_results <- as.data.frame(cox.zph(cox_model)$table)
ph_results$term <- rownames(ph_results)
rownames(ph_results) <- NULL

km_data <- cox_data

km_data$Factor1_group <- ifelse(
  km_data$Factor1 >= median(km_data$Factor1),
  "High",
  "Low"
)

km_data$Factor3_group <- ifelse(
  km_data$Factor3 >= median(km_data$Factor3),
  "High",
  "Low"
)

fit_km_f1 <- survfit(
  Surv(
    `Overall survival, days`,
    `Survival status (1, dead; 0, alive)`
  ) ~ Factor1_group,
  data = km_data
)

fit_km_f3 <- survfit(
  Surv(
    `Overall survival, days`,
    `Survival status (1, dead; 0, alive)`
  ) ~ Factor3_group,
  data = km_data
)

plot_f1 <- ggsurvplot(
  fit_km_f1,
  data = km_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by MOFA Factor 1",
  legend.title = "Factor 1",
  legend.labs = c("High", "Low")
)

plot_f3 <- ggsurvplot(
  fit_km_f3,
  data = km_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by MOFA Factor 3",
  legend.title = "Factor 3",
  legend.labs = c("High", "Low")
)

png(
  file.path(survival_dir, "KM_Overall_Survival_MOFA_Factor1.png"),
  width = 3000,
  height = 2250,
  res = 300
)
print(plot_f1)
dev.off()

pdf(
  file.path(survival_dir, "KM_Overall_Survival_MOFA_Factor1.pdf"),
  width = 10,
  height = 7.5
)
print(plot_f1)
dev.off()

png(
  file.path(survival_dir, "KM_Overall_Survival_MOFA_Factor3.png"),
  width = 3000,
  height = 2250,
  res = 300
)
print(plot_f3)
dev.off()

pdf(
  file.path(survival_dir, "KM_Overall_Survival_MOFA_Factor3.pdf"),
  width = 10,
  height = 7.5
)
print(plot_f3)
dev.off()

write.csv(
  cox_data,
  file.path(survival_dir, "Survival_Analysis_Data.csv"),
  row.names = FALSE
)

write.csv(
  cox_results,
  file.path(survival_dir, "Cox_Regression_Results.csv"),
  row.names = FALSE
)

write.csv(
  ph_results,
  file.path(survival_dir, "Cox_Proportional_Hazards_Test.csv"),
  row.names = FALSE
)

capture.output(
  summary(cox_f1),
  summary(cox_f3),
  summary(cox_model),
  cox.zph(cox_model),
  file = file.path(
    survival_dir,
    "Survival_Analysis_Console_Output.txt"
  )
)

print(cox_results)
print(ph_results)

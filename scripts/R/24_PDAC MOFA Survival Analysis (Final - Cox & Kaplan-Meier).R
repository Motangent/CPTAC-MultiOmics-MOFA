# ==============================================================================
# Script 24: PDAC MOFA Survival Analysis (Final - Cox & Kaplan-Meier)
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(survival)
  library(survminer)
})

# ۱. بارگذاری مدل MOFA اصلی (۱۰۵ نمونه، شناسه‌های C3L)
all_mofa_files <- list.files("C:/Users/darvi/Desktop/MOFA", pattern = "mofa.*\\.(hdf5|rds)$", full.names = TRUE, recursive = TRUE)
selected_factors <- NULL

for (f in all_mofa_files) {
  obj <- tryCatch({
    if (grepl("\\.rds$", f)) readRDS(f) else load_model(f)
  }, error = function(e) NULL)
  
  if (!is.null(obj)) {
    fmat <- tryCatch({
      if (inherits(obj, "MOFA")) get_factors(obj, factors = "all")[[1]] else as.matrix(obj)
    }, error = function(e) NULL)
    
    if (!is.null(fmat) && any(grepl("^C3L-", rownames(fmat)))) {
      selected_factors <- fmat
      cat(sprintf("مدل اصلی لود شد: %s (%d نمونه)\n", basename(f), nrow(fmat)))
      break
    }
  }
}
if (is.null(selected_factors)) stop("مدل MOFA با شناسه‌های C3L یافت نشد.")

factors_df <- as.data.frame(selected_factors)
factors_df$Patient_ID <- rownames(factors_df)

# ۲. بارگذاری متغیرهای بقا از فایل بالینی
clin_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data/clinical.csv"
clin_raw <- read.csv(clin_path, stringsAsFactors = FALSE, check.names = FALSE)
colnames(clin_raw) <- make.unique(colnames(clin_raw))

surv_df <- data.frame(
  Patient_ID = as.character(clin_raw$Patient_ID),
  OS_days    = suppressWarnings(as.numeric(clin_raw[["Overall survival, days"]])),
  OS_event   = suppressWarnings(as.numeric(clin_raw[["Survival status (1, dead; 0, alive)"]])),
  stringsAsFactors = FALSE
)

# ۳. ادغام و اعتبارسنجی
valid_surv <- merge(factors_df, surv_df, by = "Patient_ID") %>%
  filter(!is.na(OS_days) & !is.na(OS_event) & OS_days > 0)
valid_surv$OS_months <- valid_surv$OS_days / 30.4375

cat(sprintf("\nنمونه‌های منطبق: %d | رخدادهای فوت: %d\n\n",
            nrow(valid_surv), sum(valid_surv$OS_event == 1)))

# ۴. Cox تک‌متغیره برای تمامی فاکتورها
factor_cols <- grep("^Factor", colnames(valid_surv), value = TRUE)
cox_results <- data.frame(
  Factor       = factor_cols,
  Hazard_Ratio = NA_real_,
  CI_Lower     = NA_real_,
  CI_Upper     = NA_real_,
  P_Value      = NA_real_
)

for (i in seq_along(factor_cols)) {
  fn <- factor_cols[i]
  fit_c <- tryCatch(coxph(as.formula(paste("Surv(OS_months, OS_event) ~", fn)), data = valid_surv), error = function(e) NULL)
  if (!is.null(fit_c)) {
    s_fit <- summary(fit_c)
    cox_results$Hazard_Ratio[i] <- s_fit$coefficients[1, "exp(coef)"]
    cox_results$CI_Lower[i]     <- s_fit$conf.int[1, "lower .95"]
    cox_results$CI_Upper[i]     <- s_fit$conf.int[1, "upper .95"]
    cox_results$P_Value[i]      <- s_fit$coefficients[1, "Pr(>|z|)"]
  }
}

cox_results$FDR <- p.adjust(cox_results$P_Value, method = "fdr")
cox_results <- cox_results[order(cox_results$P_Value), ]

cat("=== Cox Proportional Hazards (Univariate) ===\n")
print(cox_results, digits = 4, row.names = FALSE)

# ۵. Kaplan-Meier برای فاکتورهای معنادار (Cox p < 0.05)
sig_factors <- cox_results$Factor[which(cox_results$P_Value < 0.05)]
cat("\nفاکتورهای معنادار:", paste(sig_factors, collapse = ", "), "\n")

km_plot_list <- list()

for (fn in sig_factors) {
  valid_surv$Factor_Group <- ifelse(valid_surv[[fn]] >= median(valid_surv[[fn]], na.rm = TRUE), "High", "Low")
  km_fit <- survfit(Surv(OS_months, OS_event) ~ Factor_Group, data = valid_surv)
  
  p_km <- ggsurvplot(
    km_fit,
    data = valid_surv,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    legend.title = paste(fn, "Level"),
    legend.labs = c("High", "Low"),
    palette = c("#E64B35", "#4DBBD5"),
    title = paste("Overall Survival Analysis by", fn),
    xlab = "Time (Months)",
    ylab = "Overall Survival Probability",
    ggtheme = theme_classic()
  )
  
  km_plot_list[[fn]] <- p_km
  print(p_km)
}

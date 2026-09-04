# ==============================================================================
# SCRIPT 23: Robust ID Mapping & Clinical Association Heatmap
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(pheatmap)
})

# 1. Load Clinical Data
clin_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data/clinical.csv"
if (!file.exists(clin_path)) {
  clin_path <- "C:/Users/darvi/Desktop/CPTAC Overview/PDAC_raw_data/clinical.csv"
}
clinical_raw <- read.csv(clin_path, stringsAsFactors = FALSE, check.names = FALSE)

# 2. Identify Mapping Column in Clinical Table
mofa_ids <- rownames(factors_matrix)
map_col <- NULL

for (col_name in colnames(clinical_raw)) {
  vals <- as.character(clinical_raw[[col_name]])
  clean_vals <- gsub("[._-]", "", toupper(trimws(vals)))
  clean_mofa <- gsub("[._-]", "", toupper(trimws(mofa_ids)))
  if (length(intersect(clean_mofa, clean_vals)) >= 10) {
    map_col <- col_name
    break
  }
}

# Check clinical_mssm.csv if not found directly
if (is.null(map_col)) {
  mssm_path <- "C:/Users/darvi/Desktop/MOFA/CPTAC Overview/PDAC_raw_data/clinical_mssm.csv"
  if (!file.exists(mssm_path)) {
    mssm_path <- "C:/Users/darvi/Desktop/CPTAC Overview/PDAC_raw_data/clinical_mssm.csv"
  }
  if (file.exists(mssm_path)) {
    mssm_df <- read.csv(mssm_path, stringsAsFactors = FALSE, check.names = FALSE)
    for (col_name in colnames(mssm_df)) {
      vals <- as.character(mssm_df[[col_name]])
      clean_vals <- gsub("[._-]", "", toupper(trimws(vals)))
      clean_mofa <- gsub("[._-]", "", toupper(trimws(mofa_ids)))
      if (length(intersect(clean_mofa, clean_vals)) >= 10) {
        clinical_raw <- merge(clinical_raw, mssm_df, by.x = "Patient_ID", by.y = "Patient_ID", all.x = TRUE)
        map_col <- col_name
        break
      }
    }
  }
}

# 3. Align Samples
if (!is.null(map_col)) {
  clean_mofa <- gsub("[._-]", "", toupper(trimws(mofa_ids)))
  clean_target <- gsub("[._-]", "", toupper(trimws(clinical_raw[[map_col]])))
  
  matched_mofa_idx <- c()
  matched_clin_idx <- c()
  
  for (i in seq_along(clean_mofa)) {
    hit <- which(clean_target == clean_mofa[i])
    if (length(hit) > 0) {
      matched_mofa_idx <- c(matched_mofa_idx, i)
      matched_clin_idx <- c(matched_clin_idx, hit[1])
    }
  }
  
  factors_sub <- factors_matrix[matched_mofa_idx, , drop = FALSE]
  clinical_sub <- clinical_raw[matched_clin_idx, , drop = FALSE]
  cat(sprintf("Successfully matched %d samples via column '%s'.\n", nrow(factors_sub), map_col))
} else {
  factors_sub <- factors_matrix
  clinical_sub <- clinical_raw[1:nrow(factors_matrix), , drop = FALSE]
  cat("Matching by row index fallback.\n")
}

# 4. Association Matrix Computation
selected_vars <- c(
  "age", "sex", "race", "bmi", "tobacco_smoking_history", "alcohol_consumption",
  "histologic_grade", "histologic_type", "tumor_site", "tumor_size_cm", 
  "tumor_focality", "tumor_necrosis", "perineural_invasion", "margin_status", "residual_tumor",
  "tumor_stage_pathological", "pathologic_staging_primary_tumor_pt", 
  "pathologic_staging_regional_lymph_nodes_pn", "number_of_lymph_nodes_positive_for_tumor_by_he_staining"
)
selected_vars <- intersect(selected_vars, colnames(clinical_sub))
continuous_vars <- c("age", "bmi", "tumor_size_cm", "number_of_lymph_nodes_positive_for_tumor_by_he_staining")

n_factors <- ncol(factors_sub)
n_vars <- length(selected_vars)
raw_pvals <- matrix(NA, nrow = n_factors, ncol = n_vars, 
                    dimnames = list(colnames(factors_sub), selected_vars))

for (v in selected_vars) {
  raw_vals <- clinical_sub[[v]]
  raw_vals[raw_vals %in% c("", "null", "NA", "Not Reported", "Unknown", "Not Evaluated: Not provided or available", "Not Applicable")] <- NA
  
  if (v %in% continuous_vars) {
    num_vals <- suppressWarnings(as.numeric(as.character(raw_vals)))
    for (f in seq_len(n_factors)) {
      valid_idx <- which(!is.na(num_vals) & !is.na(factors_sub[, f]))
      if (length(valid_idx) >= 10 && length(unique(num_vals[valid_idx])) > 2) {
        test_res <- suppressWarnings(cor.test(factors_sub[valid_idx, f], num_vals[valid_idx], method = "spearman"))
        raw_pvals[f, v] <- test_res$p.value
      }
    }
  } else {
    fac_vals <- as.factor(as.character(raw_vals))
    for (f in seq_len(n_factors)) {
      valid_idx <- which(!is.na(fac_vals) & !is.na(factors_sub[, f]))
      tbl <- table(fac_vals[valid_idx])
      valid_levels <- names(tbl[tbl >= 2])
      filt_idx <- valid_idx[fac_vals[valid_idx] %in% valid_levels]
      
      if (length(filt_idx) >= 8 && length(unique(fac_vals[filt_idx])) >= 2) {
        test_res <- suppressWarnings(kruskal.test(factors_sub[filt_idx, f] ~ fac_vals[filt_idx]))
        raw_pvals[f, v] <- test_res$p.value
      }
    }
  }
}

# 5. FDR Adjustment
fdr_pvals <- matrix(p.adjust(as.vector(raw_pvals), method = "BH"), 
                    nrow = nrow(raw_pvals), ncol = ncol(raw_pvals), 
                    dimnames = dimnames(raw_pvals))

# 6. Prepare Matrix for Heatmap (-log10 Raw P-value)
log_p <- -log10(raw_pvals)
log_p[is.na(log_p)] <- 0
max_val <- max(log_p, na.rm = TRUE)
if (max_val < 2) max_val <- 2

label_mat <- matrix("", nrow = nrow(raw_pvals), ncol = ncol(raw_pvals))
label_mat[raw_pvals < 0.05] <- "*"
label_mat[raw_pvals < 0.01] <- "**"
label_mat[!is.na(fdr_pvals) & fdr_pvals < 0.05] <- "***"

# 7. Render Heatmap in Plots Panel
pheatmap(
  mat = log_p,
  main = "MOFA Clinical Associations (-log10 Raw P-value)\n[* p<0.05, ** p<0.01, *** FDR<0.05]",
  color = colorRampPalette(c("#f7fbff", "#9ecae1", "#4292c6", "#084594", "#67000d"))(100),
  breaks = seq(0, max_val, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  display_numbers = label_mat,
  fontsize_number = 11,
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45
)

# 8. Print Console Summary
cat("\n================ CLINICAL ASSOCIATIONS (Raw P < 0.05) ================\n")
nom_hits <- which(raw_pvals < 0.05, arr.ind = TRUE)
if (nrow(nom_hits) > 0) {
  res_df <- data.frame(
    Factor = rownames(raw_pvals)[nom_hits[, 1]],
    Variable = colnames(raw_pvals)[nom_hits[, 2]],
    Raw_P = format(raw_pvals[nom_hits], scientific = TRUE, digits = 4),
    FDR = format(fdr_pvals[nom_hits], scientific = TRUE, digits = 4),
    stringsAsFactors = FALSE
  )
  print(res_df[order(as.numeric(res_df$Raw_P)), ])
} else {
  cat("No clinical associations reached Raw P < 0.05.\n")
}

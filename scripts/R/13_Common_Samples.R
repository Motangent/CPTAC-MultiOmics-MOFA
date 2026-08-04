# Step 13: Common samples across five omics (NO rerun)

rm(list = ls())

load("results/PDAC_aligned_five_omics_checkpoint.RData")

omics_objects <- list(
  Proteomics = prot_aligned,
  Transcriptomics = trans_aligned,
  Phosphoproteomics = phospho_filtered,
  miRNA = mirna_bcm_final,
  CNV = cnv_bcm_final
)

common_samples <- Reduce(
  intersect,
  lapply(omics_objects, rownames)
)

stopifnot(length(common_samples) > 0)

prot_final <- prot_aligned[common_samples, , drop = FALSE]
trans_final <- trans_aligned[common_samples, , drop = FALSE]
phospho_final <- phospho_filtered[common_samples, , drop = FALSE]
mirna_final <- mirna_bcm_final[common_samples, , drop = FALSE]
cnv_final <- cnv_bcm_final[common_samples, , drop = FALSE]

dir.create("results", recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(Patient_ID = common_samples),
  "results/PDAC_common_samples_five_omics.csv",
  row.names = FALSE
)

save(
  common_samples,
  prot_final,
  trans_final,
  phospho_final,
  mirna_final,
  cnv_final,
  file = "results/PDAC_common_five_omics_checkpoint.RData"
)

cat(
  "\nFive-omics common-sample checkpoint created:",
  length(common_samples),
  "samples.\n"
)

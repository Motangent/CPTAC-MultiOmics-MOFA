
# ============================================================
# Project: CPTAC Multi-Omics MOFA
# Cancer : Pancreatic Ductal Adenocarcinoma (PDAC)
# Step   : Load PDAC datasets
# ============================================================

library(cptac)

# Load PDAC dataset
pdac <- cptac::Pdac()

# Clinical data
clinical <- pdac$get_clinical(source = "mssm")

# Proteomics
proteomics_bcm <- pdac$get_proteomics(source = "bcm")
proteomics_umich <- pdac$get_proteomics(source = "umich")

# Phosphoproteomics
phospho_bcm <- pdac$get_phosphoproteomics(source = "bcm")
phospho_umich <- pdac$get_phosphoproteomics(source = "umich")

# Transcriptomics
transcriptomics_bcm <- pdac$get_transcriptomics(source = "bcm")
transcriptomics_broad <- pdac$get_transcriptomics(source = "broad")
transcriptomics_washu <- pdac$get_transcriptomics(source = "washu")

# CNV
cnv_bcm <- pdac$get_CNV(source = "bcm")
cnv_washu <- pdac$get_CNV(source = "washu")

# miRNA
mirna_bcm <- pdac$get_miRNA(source = "bcm")
mirna_washu <- pdac$get_miRNA(source = "washu")
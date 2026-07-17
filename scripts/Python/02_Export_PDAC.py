from pathlib import Path
import cptac

# ============================================================
# Project: CPTAC Multi-Omics MOFA
# Cancer : PDAC
# Step   : Export downloaded datasets to CSV
# ============================================================

# Load dataset
pdac = cptac.Pdac()

# Output folder
project_root = Path(__file__).resolve().parents[2]
output_dir = project_root / "data" / "raw" / "PDAC"
output_dir.mkdir(parents=True, exist_ok=True)

datasets = {
    "clinical_mssm": pdac.get_clinical(source="mssm"),
    "medical_history": pdac.get_medical_history(source="mssm"),

    "proteomics_bcm": pdac.get_proteomics(source="bcm"),
    "proteomics_umich": pdac.get_proteomics(source="umich"),

    "phosphoproteomics_bcm": pdac.get_phosphoproteomics(source="bcm"),
    "phosphoproteomics_umich": pdac.get_phosphoproteomics(source="umich"),

    "transcriptomics_bcm": pdac.get_transcriptomics(source="bcm"),
    "transcriptomics_broad": pdac.get_transcriptomics(source="broad"),
    "transcriptomics_washu": pdac.get_transcriptomics(source="washu"),

    "cnv_bcm": pdac.get_CNV(source="bcm"),
    "cnv_washu": pdac.get_CNV(source="washu"),

    "mirna_bcm": pdac.get_miRNA(source="bcm"),
    "mirna_washu": pdac.get_miRNA(source="washu"),

    "circular_rna": pdac.get_circular_RNA(source="bcm"),

    "somatic_mutation_washu": pdac.get_somatic_mutation(source="washu"),
    "somatic_mutation_harmonized": pdac.get_somatic_mutation(source="harmonized"),

    "ancestry": pdac.get_ancestry_prediction(source="harmonized"),
    "cibersort": pdac.get_cibersort(source="washu"),
    "hla_typing": pdac.get_hla_typing(source="washu"),
    "tumor_purity": pdac.get_tumor_purity(source="washu"),
    "xcell": pdac.get_xcell(source="washu"),
}

for name, df in datasets.items():
    outfile = output_dir / f"{name}.csv"
    df.to_csv(outfile)
    print(f"Exported: {outfile.name}")

def export_dataframe(name, df):
    file_path = output_dir / f"{name}.csv"
    df.to_csv(file_path)
    print(f"✓ {name} -> {file_path}")

datasets = {
    "clinical": pdac.get_clinical(),
    "medical_history": pdac.get_medical_history(),

    "proteomics_bcm": pdac.get_proteomics(source="bcm"),
    "proteomics_umich": pdac.get_proteomics(source="umich"),

    "phosphoproteomics_bcm": pdac.get_phosphoproteomics(source="bcm"),
    "phosphoproteomics_umich": pdac.get_phosphoproteomics(source="umich"),

    "transcriptomics_bcm": pdac.get_transcriptomics(source="bcm"),
    "transcriptomics_broad": pdac.get_transcriptomics(source="broad"),
    "transcriptomics_washu": pdac.get_transcriptomics(source="washu"),

    "cnv_bcm": pdac.get_CNV(source="bcm"),
    "cnv_washu": pdac.get_CNV(source="washu"),

    "mirna_bcm": pdac.get_miRNA(source="bcm"),
    "mirna_washu": pdac.get_miRNA(source="washu"),

    "circular_rna": pdac.get_circular_RNA(source="bcm"),

    "somatic_mutation_washu": pdac.get_somatic_mutation(source="washu"),
    "somatic_mutation_harmonized": pdac.get_somatic_mutation(source="harmonized"),

    "ancestry": pdac.get_ancestry_prediction(source="harmonized"),
    "cibersort": pdac.get_cibersort(source="washu"),
    "hla_typing": pdac.get_hla_typing(source="washu"),
    "tumor_purity": pdac.get_tumor_purity(source="washu"),
    "xcell": pdac.get_xcell(source="washu"),
}

for name, df in datasets.items():
    export_dataframe(name, df)

print("\nAll PDAC datasets exported successfully.")



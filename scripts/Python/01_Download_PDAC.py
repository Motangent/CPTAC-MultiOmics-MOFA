"""
Project:
CPTAC Multi-Omics MOFA

Script:
01_Download_PDAC.py

Purpose:
Download (or load if already downloaded) the CPTAC PDAC dataset
and verify that the dataset is available.
"""

import cptac

print("Loading CPTAC PDAC dataset ...")

pdac = cptac.Pdac()

print("PDAC dataset loaded successfully.")

print("\nAvailable data types:")

print(pdac.list_data())
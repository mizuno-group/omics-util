#!/usr/bin/env python3
"""
Created on 2025-03-11 (Tue) 14:08:30

0. define file paths
1. load count data
2. normalize
3. concat (lymphocyte and myeloid)
4. visualize

@author: I.Azuma
"""
import os
import numpy as np
import pandas as pd

BASE_DIR = "C:/github/omics-util/fcs_analysis"
os.chdir(BASE_DIR)

import analyzer_utils as au
import plot_utils as pu


def main():
    # 0. define file paths
    file_lymph = "./data/20250307_lymph.csv"
    file_mye = "./data/20250307_mye.csv"
    lymph_targets=["main_population","NK","B","CD4+T","CD8+T","CD4-CD8-"]
    mye_targets=["CD45+Dump-","eosinophil","neutrophil","monocyte","Kupffer","MonoMac"]
    samples=["APAP-PBS_1", "APAP-PBS_2", "APAP-PBS_3", "APAP-Ly6G_1", "APAP-Ly6G_2", "APAP-Ly6G_3", "PBS-PBS_1", "PBS-PBS_2"]

    # 1. load count data
    df_lymph = au.load_counts(filein=file_lymph, targets=lymph_targets)
    df_mye = au.load_counts(filein=file_mye, targets=mye_targets)
    df_lymph.index = samples
    df_mye.index = samples

    # 2. normalize
    df_norm_lymph = au.norm_counts(df_lymph, base_col="main_population")
    df_norm_mye = au.norm_counts(df_mye, base_col="CD45+Dump-")

    # 3. concat (lymphocyte and myeloid)
    merged_df = pd.concat([df_norm_lymph, df_norm_mye], axis=1).T

    # 4. visualize
    pu.plot_box(merged_df, immune="NK")
    pu.plot_bar(merged_df, immune="NK", sort_sample=["PBS-PBS","APAP-PBS","APAP-Ly6G"])
    pu.plot_group_scatter(merged_df, cell1='neutrophil', cell2='NK')

if __name__ == "__main__":
    main()

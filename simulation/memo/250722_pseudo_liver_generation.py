#!/usr/bin/env python3
"""
Created on 2025-07-22 (Tue) 22:58:25

@author: I.Azuma
"""
# %%
BASE_DIR = '/workspace/cluster/HDD/azuma/TopicModel_Deconv'

import pandas as pd
import scanpy as sc

import sys
sys.path.append(BASE_DIR+'/github/TRIAD')
from _utils import simulation

"""
adata = sc.read_h5ad(BASE_DIR+'/datasource/scRNASeq/LiverCellAtlas/mouseStStAll/processed/liver_adata_148202x19052.h5ad')
summary_df = pd.read_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/gs_8000x12.csv', index_col=0)

dat = simulation.LiverCellAtlas_Simulator(sample_size=8000, method='dirichlet')
dat.set_data(adata=adata)
dat.split_cell_idx(save_dir=BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/cell_idx')
#cell_idx_dict = dat.cell_idx_dict
#pd.to_pickle(cell_idx_dict, BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/cell_idx/cell_idx_dict.pkl')
"""

# create train (8000 samples) dataset
summary_df = pd.read_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/gs_8000x12.csv', index_col=0)
cell_idx_dict = pd.read_pickle(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/cell_idx/cell_idx_dict.pkl')
dat = simulation.LiverCellAtlas_Simulator(sample_size=8000, method='dirichlet')
dat.set_data(summary_df=summary_df, cell_idx_dict=cell_idx_dict)
bulk_df = dat.create_sim_bulk(pool_size=500, mode='train')
bulk_df.to_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/bulk_19052_8000.csv')

# create test (100 samples) dataset
summary_df = pd.read_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/gs_100x12.csv', index_col=0)
cell_idx_dict = pd.read_pickle(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/cell_idx/cell_idx_dict.pkl')
dat = simulation.LiverCellAtlas_Simulator(sample_size=100, method='dirichlet')
dat.set_data(summary_df=summary_df, cell_idx_dict=cell_idx_dict)
bulk_df = dat.create_sim_bulk(pool_size=500, mode='test')
bulk_df.to_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/bulk_19052_100.csv')

# %% Save
import anndata as ad
import pandas as pd

df1 = pd.read_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/bulk_19052_8000.csv', index_col=0)
df2 = pd.read_csv(BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/gs_8000x12.csv', index_col=0)
df2['ds'] = 'LiverCellAtlas'
df2['batch'] = 0

adata_x = df1.T
adata = ad.AnnData(X=adata_x)
adata.obs = df2
adata.var = pd.DataFrame(index=adata_x.columns)
OUTPUT_FILENAME = BASE_DIR+'/datasource/Simulated_Data/LiverCellAtlas/processed/dirichlet/livercellatlas_simulated.h5ad'
adata.write_h5ad(OUTPUT_FILENAME)

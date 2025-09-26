# -*- coding: utf-8 -*-
"""
Created on 2024-09-22 (Sun) 14:23:13

Reproduction of the scaden pipeline for the PBMC data.

Reference:
- pbmc_processing_10x6K.ipynb downloaded from https://figshare.com/articles/software/Publication_Figures/8234030?file=17855789

@author: I.Azuma
"""
# %%
BASE_DIR = "/workspace/mnt/cluster/HDD/azuma/TopicModel_Deconv/"

import gc
import pandas as pd
import scanpy as sc
import numpy as np
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, color_map='viridis')
sc.logging.print_versions()

# Parameters
num_top_genes = 100
out_path = BASE_DIR+"/datasource/scRNASeq/PBMCs/6k_processed/"

# %% Load the data
# Paths to three different datasets
data_path = BASE_DIR+"/datasource/scRNASeq/Scaden/pbmc6k/filtered_matrices_mex/hg19/"
# Load PBMC 6k dataset
adata = sc.read(data_path + 'matrix.mtx').T
adata.var_names = pd.read_csv(data_path + 'genes.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv(data_path + 'barcodes.tsv', header=None)[0]
adata.var_names_make_unique()
adata  # n_obs × n_vars = 5419 × 32738

# %% Preprocessing for count matrix
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=5)
adata.obs['n_counts'] = adata.X.sum(axis=1)
sc.pl.scatter(adata, x='n_counts', y='n_genes')

adata = adata[adata.obs['n_genes'] < 2000, :]
adata = adata[adata.obs['n_counts'] < 6000, :]
adata.raw = sc.pp.log1p(adata, copy=True)

# Normalize
sc.pp.normalize_per_cell(adata)

# Save the complete matrix
df = pd.DataFrame(adata.X.todense())
df.columns = adata.var.index
df.to_csv(out_path + "data6k_norm_counts_all.txt", sep="\t")

gc.collect()
# %%
filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.filter_genes_dispersion(filter_result)

adata = adata[:, filter_result.gene_subset]
sc.pp.normalize_per_cell(adata)
adata.X_noscale = adata.X
sc.pp.log1p(adata)
sc.pp.regress_out(adata, ['n_counts'])
sc.pp.scale(adata)

gc.collect()

# calculate UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)  # NOTE: this is slow
sc.tl.louvain(adata, resolution=0.6)
sc.pl.umap(adata, color='louvain')

# visualization
sc.pl.umap(adata, color=['IL7R', 'LYZ', 'MS4A1', 'GNLY', 'FCER1A', 'FCGR3A', 'CST3', 'CD8A', 'CCL5'])
new_celltypes = ['CD4Tcells','Bcells', 'Monocytes', 'CD8Tcells','Monocytes2', 'NK', 'Dendritic']
adata.rename_categories('louvain', new_celltypes)
sc.pl.umap(adata, color='louvain', legend_loc='on data')

# Make new category
celltypes = pd.DataFrame(adata.obs['louvain'])
celltypes.louvain.replace(['Monocytes2'],
                          ['Monocytes'], inplace=True)
adata.obs['celltype'] = celltypes.louvain
sc.pl.umap(adata, color='celltype', legend_loc='on data')

# Save celltypes
celltypes = pd.DataFrame(adata.obs['celltype'])
celltypes.columns = ['Celltype']
celltypes.to_csv(out_path + "data6k_celltypes.txt", sep="\t")

# %% Rank the genes
sc.tl.rank_genes_groups(adata, 'celltype')
top_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
top_genes.head(num_top_genes)

# Save top genes
top_genes = top_genes.head(num_top_genes)
top_genes_flat = pd.DataFrame(top_genes.values.flatten())
top_genes_flat.to_csv(out_path + "data6k_top_genes_" + str(num_top_genes) +".txt", sep="\t")

# Save adata
adata.write(out_path + "data6k_processed.h5ad")

# %% reference data processing
celltype_info = pd.read_table(BASE_DIR+'/datasource/scRNASeq/PBMCs/6k_processed/data6k_celltypes.txt',index_col=0)
celltype_dict = dict(zip(celltype_info.index, celltype_info['Celltype']))

# Paths to three different datasets
data_path = BASE_DIR+"/datasource/scRNASeq/Scaden/pbmc6k/filtered_matrices_mex/hg19/"
# Load PBMC 6k dataset
adata = sc.read(data_path + 'matrix.mtx').T
adata.var_names = pd.read_csv(data_path + 'genes.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv(data_path + 'barcodes.tsv', header=None)[0]
adata.var_names_make_unique()
adata  # n_obs × n_vars = 5419 × 32738

# reflect celtype information
adata.obs['celltype'] = [celltype_dict[x] if x in celltype_dict.keys() else 'Unknown' for x in adata.obs_names]

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=5)
adata.obs['n_counts'] = adata.X.sum(axis=1)
sc.pl.scatter(adata, x='n_counts', y='n_genes')

adata = adata[adata.obs['n_genes'] < 2000, :]
adata = adata[adata.obs['n_counts'] < 6000, :]
#adata.raw = sc.pp.log1p(adata, copy=True)

# Normalize
sc.pp.normalize_per_cell(adata)

# Cell type specific gene expression
mean_expression_per_celltype = pd.DataFrame(adata.X.toarray(),index=adata.obs.index, columns=adata.var.index).groupby(adata.obs['celltype']).mean()
mean_expression_per_celltype.to_csv(out_path + "data6k_6x13137_reference.csv")

# %%

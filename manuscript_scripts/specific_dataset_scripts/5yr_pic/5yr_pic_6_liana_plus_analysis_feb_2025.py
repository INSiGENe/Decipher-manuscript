import numpy as np
import pandas as pd
import scanpy as sc

import plotnine as p9

import liana as li
import decoupler as dc
import omnipath as op

# Import DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import corneto as cn

import random

import os


# Set the seed
random.seed(123)  # You can use any number as the seed
# Set the seed
np.random.seed(123)  # Use the same or a different seed as needed

dataset_path = "results/5yr_pic"
# Specify the path to the folder you want to loop through
output_data_path = os.path.join(dataset_path,"liana/data")
h5ad_seurat_file = os.path.join(output_data_path,"seurat_object_oi.h5ad")
figures_data_path = os.path.join(dataset_path,"figures")

adata = sc.read_h5ad(filename = h5ad_seurat_file)
adata

sample_key = "orig.ident"
groupby = 'cluster'
condition_key = 'condition'

# filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


sc.pp.neighbors(adata)
sc.tl.umap(adata)
# Show pre-computed UMAP
#sc.pl.umap(adata, color=[condition_key, sample_key, groupby], frameon=False, ncols=2)

#pseudobulk samples (sample-pair)
pdata = dc.get_pseudobulk(
    adata,
    sample_col=sample_key,
    groups_col=groupby,
    mode='sum',
    min_cells=10,
    min_counts=10000
)


dc.plot_psbulk_samples(pdata, groupby=[sample_key, groupby], figsize=(11, 4))

dea_results = {}
for cell_group in pdata.obs[groupby].unique():
    # Select cell profiles
    ctdata = pdata[pdata.obs[groupby] == cell_group].copy()
    # Obtain genes that pass the edgeR-like thresholds
    # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
    genes = dc.filter_by_expr(ctdata,
                              group=condition_key,
                              min_count=5, # a minimum number of counts in a number of samples
                              min_total_count=10 # a minimum total number of reads across samples
                              )
    # Filter by these genes
    ctdata = ctdata[:, genes].copy()
    if ctdata.obs['condition'].nunique() >= 2:
        pass  
    else:
        continue
    # Build DESeq2 object
    # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
    dds = DeseqDataSet(
        adata=ctdata,
        design_factors=condition_key,
        ref_level=[condition_key, 'ctrl'], # set control as reference
        refit_cooks=True,
        n_cpus=None,
    )
    # Compute LFCs
    dds.deseq2()
    # Contrast between stim and ctrl
    stat_res = DeseqStats(dds, contrast=[condition_key, 'stim', 'ctrl'])
    # Compute Wald test
    stat_res.summary()
    # Shrink LFCs
    stat_res.lfc_shrink(coeff='condition_stim_vs_ctrl') # {condition_key}_cond_vs_ref
    dea_results[cell_group] = stat_res.results_df

# concat results across cell types
dea_df = pd.concat(dea_results)

dea_df = dea_df.reset_index().rename(columns={'level_0': groupby,"level_1":"index"})
dea_df.head()

adata = adata[adata.obs[condition_key]=='stim'].copy()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

dea_df = dea_df.set_index('index')

arr = dea_df.index.isin(adata.var_names.values.tolist())

filtered_dea_df = dea_df[arr]

lr_res = li.multi.df_to_lr(adata,
                           dea_df=filtered_dea_df,
                           resource_name='connectomedb2020',
                           expr_prop=0.1, # calculated for adata as passed - used to filter interactions
                           groupby=groupby,
                           stat_keys=['stat', 'pvalue', 'padj'],
                           use_raw=False,
                           complex_col='stat', # NOTE: we use the Wald Stat to deal with complexes
                           verbose=True,
                           return_all_lrs=False,
                           )

lr_res = lr_res.sort_values("interaction_stat", ascending=False, key=abs)
lr_res.head()


lr_res.to_csv(os.path.join(output_data_path,"liana_p_interaction_results.csv"))
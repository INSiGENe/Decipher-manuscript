import anndata as ad
from scipy.io import mmwrite
import gc
import pandas as pd
import random
import numpy as np

# Set the random seed
random.seed(123)
np.random.seed(123)


#first decompress the .h5ad.gz file
# Load the .h5ad file with memory-mapping
#backed='r'
adata = ad.read_h5ad('manuscript_analysis/lupus/pre_processing/GSE174188_CLUES1_adjusted.h5ad')
#adata.obs['Status'].unique().tolist()


del adata.uns['umap']
del adata.uns['rank_genes_groups']
del adata.uns['neighbors']
del adata.uns['pca']
del adata.obsp['distances']
del adata.obsp['connectivities']
del adata.obsm['X_pca']
del adata.obsm['X_umap']

subset_mask = (adata.obs['pop_cov'] == 'Asian') & \
              (adata.obs['Sex'] == 'Female') & \
              (adata.obs['Status'].isin(['Healthy', 'Managed']) & \
               (adata.obs['Processing_Cohort'] == '4.0')).copy()

# Subset data
adata_subset = adata[subset_mask].copy()
del adata
gc.collect()

this_raw = adata_subset.raw.X.copy()

non_zero_rows = np.array((this_raw != 0).sum(axis=0)).flatten() > 0
filtered_matrix = this_raw[:,non_zero_rows]

non_zero_cols = np.array((filtered_matrix != 0).sum(axis=1)).flatten() >= 200
filtered_matrix = filtered_matrix[non_zero_cols,:]


# Extract gene names for filtered columns
gene_names = adata_subset.raw.var_names[non_zero_rows]


# Extract metadata for filtered rows
filtered_obs = adata_subset.obs.iloc[non_zero_cols]

# Save gene names and metadata
pd.Series(gene_names).to_csv('manuscript_analysis/lupus/pre_processing/gene_names.csv', index=False)
filtered_obs.to_csv('manuscript_analysis/lupus/pre_processing/filtered_metadata.csv', index=True)

mmwrite(f'manuscript_analysis/lupus/pre_processing/output_matrix.mtx', filtered_matrix)

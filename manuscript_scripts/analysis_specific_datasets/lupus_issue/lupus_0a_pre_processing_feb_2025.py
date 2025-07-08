import anndata as ad
from scipy.io import mmwrite
import gc
import pandas as pd
import random
import numpy as np
import os

# Set the random seed
random.seed(123)
np.random.seed(123)

#raw data download from here https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2


#first decompress the .h5ad.gz file
# Load the .h5ad file with memory-mapping
#backed='r'
adata = ad.read_h5ad('data/lupus/lupus_sc_data.h5ad')
#adata.obs['Status'].unique().tolist()


subset_mask = (adata.obs['self_reported_ethnicity'] == 'Asian') & \
              (adata.obs['sex'] == 'female') & \
              (adata.obs['disease_state'].isin(['na', 'managed']) & \
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
gene_names = adata_subset.var['feature_name'].values[non_zero_rows]

# Extract cell names (row indices from adata.obs)
cell_names = adata_subset.obs_names[non_zero_cols]


# Extract metadata for filtered rows
filtered_obs = adata_subset.obs.iloc[non_zero_cols]

# Define the directory path
output_dir = "results/lupus/pre_processing"

# Create the directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Save gene names and metadata
pd.Series(gene_names).to_csv('results/lupus/pre_processing/gene_names.csv', index=False)
filtered_obs.to_csv('results/lupus/pre_processing/filtered_metadata.csv', index=True)
# Save cell names to a file
pd.Series(cell_names).to_csv(f"{output_dir}/cell_names.csv", index=False)

mmwrite(f'results/lupus/pre_processing/output_matrix.mtx', filtered_matrix)

import anndata as ad
from scipy.io import mmwrite
import gc
import pandas as pd
import random
import numpy as np
import os
import json
import argparse
import sys

#functions
def is_integer_subset(matrix, n=100):
    # Subset the matrix
    sub = matrix[:n, :n]
    arr = sub.A if hasattr(sub, "A") else sub
    # Check if values are close to integers
    return np.all(np.isclose(arr, np.round(arr)))

# Set the random seed
random.seed(123)
np.random.seed(123)

# ---- Argument parser to select dataset key ----
parser = argparse.ArgumentParser(description="Preprocess h5ad file using config.json")
parser.add_argument("dataset_key", type=str, help="Key for the dataset in config.json (e.g., cz_rcc)")
args = parser.parse_args()

# ---- Load config ----
with open("scripts/config.json") as f:
    config = json.load(f)

cfg = config[args.dataset_key]
preproc = cfg["pre_processing"]
input_path = os.path.join(preproc["input_path"], "dataset.h5ad")
output_dir = os.path.join(preproc["output_path"], preproc["step"])
os.makedirs(output_dir, exist_ok=True)

# ---- Load data ----
adata = ad.read_h5ad(input_path)

# ---- Flexible subsetting logic ----
subset_mask = pd.Series(True, index=adata.obs.index)

if "subset_logic" in preproc:
    for field, values in preproc["subset_logic"].items():
        if isinstance(values, list):
            subset_mask &= adata.obs[field].isin(values)
        else:
            subset_mask &= adata.obs[field] == values

adata_subset = adata[subset_mask].copy()
del adata
gc.collect()

# ---- Extract raw counts and filter ----
# Check if adata.raw exists
# Determine which count matrix to use
used_raw = False
if adata_subset.raw is not None and hasattr(adata_subset.raw, "X"):
    matrix = adata_subset.raw.X
    if is_integer_subset(matrix):
        print("Using adata.raw.X (appears integer in subset)")
        this_raw = matrix.copy()
        used_raw = True
    else:
        print("Warning: adata.raw.X contains non-integer values in subset. Exiting.")
        sys.exit(1)
else:
    matrix = adata_subset.X
    if is_integer_subset(matrix):
        print("Using adata.X (appears integer in subset)")
        this_raw = matrix.copy()
    else:
        print("Warning: adata.X contains non-integer values in subset. Exiting.")
        sys.exit(1)

non_zero_genes = np.array((this_raw != 0).sum(axis=0)).flatten() > 0
filtered_matrix = this_raw[:, non_zero_genes]

min_genes = preproc.get("min_genes_per_cell", 0)
non_zero_cells = np.array((filtered_matrix != 0).sum(axis=1)).flatten() >= min_genes
filtered_matrix = filtered_matrix[non_zero_cells, :]

# ---- Extract and save metadata ----
if used_raw:
    gene_names = adata_subset.raw.var_names[non_zero_genes]
else:
    gene_names = adata_subset.var_names[non_zero_genes]
cell_names = adata_subset.obs_names[non_zero_cells]
filtered_obs = adata_subset.obs.iloc[non_zero_cells]

pd.Series(gene_names).to_csv(f"{output_dir}/gene_names.csv", index=False)
filtered_obs.to_csv(f"{output_dir}/filtered_metadata.csv", index=True)
pd.Series(cell_names).to_csv(f"{output_dir}/cell_names.csv", index=False)
mmwrite(f"{output_dir}/output_matrix.mtx", filtered_matrix)

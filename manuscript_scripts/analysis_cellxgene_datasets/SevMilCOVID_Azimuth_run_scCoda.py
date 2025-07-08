import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import anndata


import random
import numpy as np
import tensorflow as tf

random.seed(123)
np.random.seed(123)
tf.random.set_seed(123)


# Read the counts and metadata
counts = pd.read_csv("counts.csv", index_col=0)  # shape: genes x cells
obs = pd.read_csv("obs.csv", index_col=0)
var = pd.read_csv("var.csv")

# Transpose counts to (cells x genes)
counts = counts.T

# Create AnnData object
adata = anndata.AnnData(X=counts.values, obs=obs, var=var.set_index("gene_id"))

# Summarize counts per sample and cell type
# Adjust 'celltype' and 'sample_column' to your metadata fields
celltype_column = "predicted.celltype.l2"
sample_column = "sample_id"  # or replace with your column name, e.g., "orig.ident"
condition_column = "severity_group"  # adjust as needed

# Create count table per sample × cell type
counts = adata.obs.groupby(["sample_id", "predicted.celltype.l2"]).size().unstack(fill_value=0).reset_index()

# Add severity_group per sample
counts["severity_group"] = adata.obs.groupby("sample_id")["severity_group"].first().values

#Severe vs. Healthy
subset_severe = counts[counts["severity_group"].isin(["Severe", "Healthy"])]
data_severe = dat.from_pandas(subset_severe.drop(columns=["sample_id"]), covariate_columns=["severity_group"])
data_severe.obs["sample_id"] = subset_severe["sample_id"].values

model_severe = mod.CompositionalAnalysis(data_severe, formula="severity_group", reference_cell_type="automatic")
results_severe = model_severe.sample_hmc()
results_severe.set_fdr(est_fdr=0.05)

print("Severe vs. Healthy")
print(results_severe.summary())
results_severe.save("sccoda_severe_vs_healthy.pkl")

effects_df = results_severe.effect_df
effects_df.to_csv("data/SevMilCOVID/results_sccoda_severe_vs_healthy.csv")


#Moderate vs. Healthy
subset_moderate = counts[counts["severity_group"].isin(["Moderate", "Healthy"])]
data_moderate = dat.from_pandas(subset_moderate.drop(columns=["sample_id"]), covariate_columns=["severity_group"])
data_moderate.obs["sample_id"] = subset_moderate["sample_id"].values

model_moderate = mod.CompositionalAnalysis(data_moderate, formula="severity_group", reference_cell_type="automatic")
results_moderate = model_moderate.sample_hmc()
results_moderate.set_fdr(est_fdr=0.05)

print("Moderate vs. Healthy")
print(results_moderate.summary())
results_moderate.save("sccoda_moderate_vs_healthy.pkl")

effects_df_mod = results_moderate.effect_df
effects_df_mod.to_csv("data/SevMilCOVID/results_sccoda_moderate_vs_healthy.csv")




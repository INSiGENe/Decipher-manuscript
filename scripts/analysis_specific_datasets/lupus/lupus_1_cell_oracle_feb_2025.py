import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import random
import celloracle as co
import glob
co.__version__

# Set the seed
random.seed(123)  # You can use any number as the seed
# Set the seed
np.random.seed(123)  # Use the same or a different seed as needed

#parameters
data_folder = "results/lupus"
this_species = "human"

celloracle_folder = os.path.join(data_folder,"cellOracle")
os.makedirs(celloracle_folder,exist_ok=True)

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

save_folder = os.path.join(celloracle_folder,"figures")
os.makedirs(save_folder, exist_ok=True)
os.makedirs(os.path.join(save_folder,"pca_figures"), exist_ok=True)

# Specify the path to the folder you want to loop through
output_data_path = os.path.join(data_folder,"cellOracle/data")
os.makedirs(output_data_path,exist_ok=True)
os.makedirs(os.path.join(output_data_path,"GRN"),exist_ok=True)
h5ad_folder = os.path.join(data_folder,"pre_processing/h5ad_by_cluster")

# Check the 'species' parameter and load the corresponding base GRN
if this_species == "human":
    base_GRN = co.data.load_human_promoter_base_GRN()
elif this_species == "mouse":
    base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()
else:
    raise ValueError("Error: 'species' parameter must be 'human' or 'mouse'.")


for file in glob.glob(os.path.join(h5ad_folder,"*.h5ad")):
    print(file)
    adata = sc.read_h5ad(file)

    cluster_name = adata.obs["cluster"].unique()[0]
  
    # Only consider genes with more than 1 count
    sc.pp.filter_genes(adata, min_counts=1)

    # Normalize gene expression matrix with total UMI count per cell
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

    # Select top 2000 highly-variable genes
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                #flavor='cell_ranger',
                                                n_top_genes=3000,
                                                log=False)

    # Subset the genes
    adata = adata[:, filter_result.gene_subset]

    # Renormalize after filtering
    sc.pp.normalize_per_cell(adata)
    # keep raw cont data before log transformation
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()

    # Log transformation and scaling
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    #filename = f"{cluster_name}.png"
    #sc.pl.pca(adata,color = "condition", save=filename,title=cluster_name)

    # Instantiate Oracle object
    oracle = co.Oracle()

    # In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
    adata.X = adata.layers["raw_count"].copy()

    # Instantiate Oracle object.
    oracle.import_anndata_as_raw_count(adata=adata,
                                    cluster_column_name="cluster",
                                    embedding_name="X_pca")
    
    # You can load TF info dataframe with the following code.
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    # Perform PCA
    oracle.perform_PCA()

    # Select important PCs
    plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    plt.axvline(n_comps, c="k")
    #plt.show()
    
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = int(0.025*n_cell)
    print(f"Auto-selected k is :{k}")

    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=6)
    
    links = oracle.get_links(
        cluster_name_for_GRN_unit="cluster", 
        alpha=10,
        verbose_level=10)
    
    filename = cluster_name + ".csv"
 
    # Save as csv
    links.links_dict[cluster_name].to_csv(os.path.join(output_data_path,"GRN",filename))
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
    links.plot_scores_as_rank(cluster=cluster_name, n_gene=30, save=f"{save_folder}/ranked_score")



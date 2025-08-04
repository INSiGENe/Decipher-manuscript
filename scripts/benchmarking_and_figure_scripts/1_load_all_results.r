library(devtools)
load_all()
library(data.table)

set.seed(123)

# Define datasets to benchmark/visualize
datasets <- list(
  "5yr_pic"   = "results/5yr_pic",
  "bcg"  = "results/BCG",
  "cord_pic"  = "results/cord_pic",
  "covid"   = "results/covid",
  "erp"  = "results/ERP",
  "lupus"  = "results/lupus",
  "sepsis"   = "results/sepsis",
  "tnbc"  = "results/TNBC",
  "cz_influenza" = "results/cz_influenza",
  "cz_hpap_t1d_islets" = "results/cz_hpap_t1d_islets",
  "cz_hnscc_hpv" = "results/cz_hnscc_hpv",
  "cz_human_kidney_v1.5" = "results/cz_human_kidney_v1.5",
  "cz_cf_bronchial_biopsy" = "results/cz_cf_bronchial_biopsy",
  "SevCOVID_Azimuthl2" = "results/SevCOVID_Azimuthl2",
  "MilCOVID_Azimuthl2" = "results/MilCOVID_Azimuthl2"
)


# Initialize empty lists for the three outputs
results_preprocessed   <- list()  # will store results_to_compare_full for each dataset
results_for_correlation <- list()  # will store results_to_compare for correlation
results_for_comparison  <- list()  # will store results_to_compare for method comparison

# Loop over each dataset
for (ds in names(datasets)) {
  # --- Here you would load/process your data for the current dataset ---
  # Define file paths and directories for the current dataset
  dataset_path <- datasets[[ds]]
  pre_processing_filepath <- file.path(dataset_path, "pre_processing")
  meta_path <- "manuscript_analysis/data_for_meta_comparisons" #is this path ok?
  output_figures_filepath <- file.path(dataset_path, "figures")
  reference_filepath <- "reference_data"
  nichenet_reference_filepath <- file.path("reference_data", "nichenet")
  decipher_filepath <- file.path(dataset_path, "data")
  nichenet_filepath <- file.path(dataset_path, "nichenet/data")
  connectome_filepath <- file.path(dataset_path, "connectome/data")
  liana_filepath <- file.path(dataset_path, "liana/data")
  natmi_filepath <- file.path(dataset_path, "natmi/data")
  cytosig_filepath <- file.path(dataset_path, "cytosig/0_outputs")
  # Create meta directory if needed
  dir.create(meta_path, recursive = TRUE, showWarnings = FALSE)
  ## FIGURE 2: Load Data ----
  liana_results <- safe_load(file.path(liana_filepath, "liana_p_interaction_results.csv"), header = TRUE, row.names = 1)
  nichenet_results <- safe_load(file.path(nichenet_filepath, "nichenet_results.rds"))
  nichenet_prior_table_all_clusters <- safe_load(file.path(nichenet_filepath, "prior_table_all_clusters.rds"))
  decipher_results <- safe_load(file.path(decipher_filepath, "decipher_scores_by_cluster.rds"))
  connectome_results <- safe_load(file.path(connectome_filepath, "connectome_results.rds"))
  natmi_results_all <- safe_load(file.path(natmi_filepath, "diff/Delta_edges_lrc2p/All_edges_mean.csv"))

  ## Data Pre-processing ----
  natmi_results_pre_processed <- preProcessNATMI(natmi_results_all)
  connectome_results_pre_processed <- preProcessConnectome(connectome_results)
  liana_pre_processed <- preProcessLIANA(liana_results)
  decipher_pre_processed <- preProcessDecipher(decipher_results)
  nichenet_pre_processed <- preProcessNicheNet(nichenet_prior_table_all_clusters)  

  # ----- 1. Preprocessed Results -----
  results_to_compare_full <- list()

  results_to_compare_full <- add_if_not_null(results_to_compare_full, "NicheNet", nichenet_pre_processed)
  results_to_compare_full <- add_if_not_null(results_to_compare_full, "Decipher", decipher_pre_processed)
  results_to_compare_full <- add_if_not_null(results_to_compare_full, "LIANA+", liana_pre_processed)
  results_to_compare_full <- add_if_not_null(results_to_compare_full, "NATMI", natmi_results_pre_processed)
  results_to_compare_full <- add_if_not_null(results_to_compare_full, "Connectome", connectome_results_pre_processed)

  results_preprocessed[[ds]] <- results_to_compare_full

  # ----- 2. Results for Correlation Analysis -----
  results_to_compare_correlation <- lapply(names(results_to_compare_full), function(name) {
    prepareForCorrelation(name, results_to_compare_full[[name]])
  })
  names(results_to_compare_correlation) <- names(results_to_compare_full)
  results_for_correlation[[ds]] <- results_to_compare_correlation

  liana_results_for_comparison <- prepareLianaForCytosigComparison(liana_results)
  nichenet_results_for_comparison <- generateComparisonObjectFromNicheNet(nichenet_results)
  names(nichenet_results_for_comparison) <- names(nichenet_results)
  decipher_results_for_comparison <- lapply(decipher_results, "renameDecipherScore")
  connectome_results_for_comparison <- prepareConnectomeForCytosigComparison(connectome_results)
  natmi_results_for_comparison <- prepareNatmiForCytosigComparison(natmi_results_all)

  # ----- 3. Results for Method Comparison -----
  results_to_compare_comparison <- list()

  results_to_compare_comparison <- add_if_not_null(results_to_compare_comparison, "NicheNet", nichenet_results_for_comparison)
  results_to_compare_comparison <- add_if_not_null(results_to_compare_comparison, "Decipher", decipher_results_for_comparison)
  results_to_compare_comparison <- add_if_not_null(results_to_compare_comparison, "LIANA+", liana_results_for_comparison)
  results_to_compare_comparison <- add_if_not_null(results_to_compare_comparison, "NATMI", natmi_results_for_comparison)
  results_to_compare_comparison <- add_if_not_null(results_to_compare_comparison, "Connectome", connectome_results_for_comparison)

  results_for_comparison[[ds]] <- results_to_compare_comparison
}

# ==== libraries ====
library(ggplot2)
library(dplyr)     # Or library(data.table) if you prefer
library(purrr)     # For map functions
library(stringr)   # For string manipulation if needed
library(patchwork) 
library(ggridges)  # For density ridgeline plots
library(tidyr)
library(Seurat)
library(ggrepel)
library(tibble)
library(scales) 
library(reshape2)   # For reshaping matrix to long format
library(gridExtra)  # For arranging multiple heatmaps in a grid
library(pROC)
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)          # for pseudo_log_trans
library(ggrepel)
library(ggnewscale)
library(ggbeeswarm) 
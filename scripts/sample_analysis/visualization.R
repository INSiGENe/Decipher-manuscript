#Hi mate, so just doing a vignette on how to visualize the data you produced
#running this in root of the Decipher folder
#first we load the necessary librares
renv::restore()
library(devtools)
load_all()
library(data.table)

set.seed(123)

# Define datasets to benchmark/visualize
# Here it might only be one dataset (depends if you did one or multiple comparisons). 
# You can all it whatever you want but I tend to use the same folder_name

dataset_path <- "results/sample_analysis"
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- file.path("reference_data")
output_data_filepath <- file.path(dataset_path,"data")
figures_folder <- "figures_25_09_25"
#meta_path <- "manuscript_analysis/data_for_meta_comparisons" #is this path ok?
output_figures_filepath <- file.path(dataset_path, "figures")

# Create meta directory if needed
#dir.create(meta_path, recursive = TRUE, showWarnings = FALSE)
## FIGURE 2: Load Data ----
# Load and pre-process the data
decipher_results <- safe_load(file.path(output_data_filepath, "decipher_scores_by_cluster.rds"))
decipher_pre_processed <- preProcessDecipher(decipher_results)

results_preprocessed <- list()
results_preprocessed <- add_if_not_null(results_preprocessed, "Decipher", decipher_pre_processed)

# ----- 2. Results for Correlation Analysis -----
results_for_correlation <- lapply(names(results_preprocessed), function(name) {
  prepareForCorrelation(name, results_preprocessed[[name]])
})
names(results_for_correlation) <- names(results_preprocessed)

decipher_results_for_comparison <- lapply(decipher_results, "renameDecipherScore")
results_for_comparison <- list()
results_for_comparison <- add_if_not_null(results_for_comparison, "Decipher", decipher_results_for_comparison)

#ok great, that should've loaded the results so we can actually visualize them, now there's several graphs you might be interested in...
#first a volcano plot of TF activity by cell-type
selected_cluster <- names(decipher_results)[1]
differentialTFActivityVolcanoPlotByCluster(
  output_data_filepath,
  selected_cluster,
  figures_folder,
  p_threshold = 0.01,
  fc_threhsold = 2.883,
  output_filename = "diff_regulong_act_volcano_plot"
)

#Then a Decipher plot
plotDecipherPrioritizedMap(dataset_path,top_n=6,priority_receiver_cells = selected_cluster,dataset_name="sample_dataset", width=21,height=11)
#note that this doesn't currently write to the figures folder
#next a visualization of TFs tied to upstream LR pairs
plotLRTFHeatmap(
  output_data_filepath,
  selected_cluster,
  figures_folder,
  p_threshold = 0.01,
  output_name = "lr_tf_heatmap",
  n_interactions = 10,
  min_abs_decipher_score = 0.4,
  min_delta_pagoda = 2
)

#sorted TF heatmap

# Define the specific cell types (receiver cells) you want to analyze
selected_clusters <- names(decipher_results)
top_n_regulons <- 10
clusters_per_group_in_output <- 1 # Adjust as needed (e.g., 1, 2, 3)

# Dynamically load data based on conditions defined above
file_path <- file.path(output_data_filepath, "regulon_deltas_by_cluster.rds")
regulon_deltas_by_cluster <- load_regulon_data(file_path, selected_clusters)

# Calculate Global Color Scale Limit
absolute_max <- find_absolute_max_simple(regulon_deltas_by_cluster)
if(is.na(absolute_max) || absolute_max <= 0) {
    warning("Could not determine a valid absolute max deltaPagoda. Setting scale limit to 1.")
    absolute_max <- 1
}

top_n = 15

for (cell_type in selected_clusters){
  p <- plotClusterSortedTFHeatmap(
    selected_cluster = cell_type,
    top_n = top_n,
    deltas_by_cluster = regulon_deltas_by_cluster,
    global_max = absolute_max   # optional: fix color scale across multiple plots
  )

  file_name <- paste0("sorted_tf_deltas_",cell_type, "_", top_n, ".png")
  ggsave(file.path(figures_folder,file_name), p, width = 8, height = 3)
  file_name_csv <- paste0("sorted_tf_deltas_",cell_type, "_", top_n, ".csv")
  write.csv(p$data, file.path(figures_folder,file_name_csv), row.names = TRUE) 
}


#networks with pubmed prioritization
selected_clusters <- names(decipher_results)

n_pubmed <- 40
output_data_filepath

plot_pubmed_tg_heatmaps(
  output_data_filepath,
  selected_clusters[1],
  n_pubmed
)

##############
#network plots
##############

decipher_scores <- readRDS(file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
decipher_scores_by_regulon_and_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
feature_statistics <- readRDS(file.path(output_data_filepath,"feature_statistics.rds"))


# 4. Specify Parameters and Generate Networks

# Define target receiver clusters and sender cell types
target_clusters <- selected_clusters[1]
sender_cts <- selected_clusters

#calculate global stats
# GLOBAL SCALING VALUES
# TF deltaPagoda
global_deltaPagoda_max <- max(sapply(regulon_deltas_by_cluster[target_clusters], function(x) max(x$deltaPagoda, na.rm = TRUE)))

# Receptor→TF imp.perm * sign(spearman.cor)
global_receptor_tf_col_max <- max(
  sapply(c(decipher_scores_by_regulon_and_cluster[target_clusters]), function(cluster_df) {
    if (!is.null(cluster_df)) {
      df <- cluster_df %>% mutate(col = imp.perm * sign(spearman.cor))
      max(abs(df$col), na.rm = TRUE)
    } else {
      0
    }
  })
)

# Sender→Ligand frac.normalized.counts

global_sender_ligand_max <- max(
  feature_statistics %>% filter(cluster %in% target_clusters) %>% pull(sum.counts) / feature_statistics %>% filter(cluster %in% target_clusters) %>% pull(n.cell),
  na.rm = TRUE
)
#global_sender_ligand_max <- 1


# Ligand→Receptor decipher_score
global_decipher_score_max <- max(
  sapply(decipher_scores[target_clusters], function(x) max(abs(x$decipher_score), na.rm = TRUE)))

# Generate network plots for Severe condition
for (cl in target_clusters) {

  generate_network_plot("sample_analysis", cl,
                      decipher_scores, 
                      decipher_scores_by_regulon_and_cluster,
                      regulon_deltas_by_cluster, 
                      feature_statistics,
                      sender_cts, 
                      figures_folder,
                      top_interactions = NULL,
                      global_deltaPagoda_max = global_deltaPagoda_max,
                      global_receptor_tf_col_max = global_receptor_tf_col_max,
                      global_sender_ligand_max = global_sender_ligand_max/1.2,
                      global_decipher_score_max = global_decipher_score_max/1.5,
                      n_top_regulons = 10)

                        }

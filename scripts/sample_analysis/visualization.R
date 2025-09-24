#Hi mate, so just doing a vignette on how to visualize the data you produced
#running this in root of the Decipher folder
#first we load the necessary librares
library(devtools)
load_all()
library(data.table)

set.seed(123)

# Define datasets to benchmark/visualize
# Here it might only be one dataset (depends if you did one or multiple comparisons). 
# You can all it whatever you want but I tend to use the same folder_name

dataset_path <- "results/5yr_pic"
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- file.path("reference_data")
output_data_filepath <- file.path(dataset_path,"data")
figures_folder <- "figures_04_08_2025"
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
results_for_comparison <- add_if_not_null(results_to_compare_comparison, "Decipher", decipher_results_for_comparison)

#ok great, that should've loaded the results so we can actually visualize them, now there's several graphs you might be interested in...
#first a volcano plot of TF activity by cell-type
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
selected_clusters <- c("CD4")
top_n_regulons <- 10
clusters_per_group_in_output <- 1 # Adjust as needed (e.g., 1, 2, 3)

# Dynamically load data based on conditions defined above
file_path <- file.path(output_data_filepath, "regulon_deltas_by_cluster.rds")
regulon_deltas_by_cluster <- load_regulon_data(file_path, selected_clusters)

# Calculate Global Color Scale Limit
absolute_max <- find_absolute_max(regulon_deltas_by_cluster)

if(is.na(absolute_max) || absolute_max <= 0) {
    warning("Could not determine a valid absolute max deltaPagoda. Setting scale limit to 1.")
    absolute_max <- 1
}


absolute_max <- find_absolute_max(regulon_deltas_by_cluster)

for (cell_type in selected_clusters){
  for(sort_by in c("case","control")){
    top_n = 20
    p <- plotConditionSortedClusterSpecificTFHeatmaps(
      cell_type   = cell_type,
      sort_by     = sort_by,
      top_n       = top_n,
      deltas_list = regulon_deltas_by_cluster,
      global_max = absolute_max
    )

    file_name <- paste0("figure_4f_",cell_type,"_sorted_by_", sort_by, "_", top_n, ".png")
    ggsave(file.path(figures_folder,file_name), p, width = 8, height = 3)
    file_name_csv <- paste0("figure_4f_",cell_type,"_sorted_by_", sort_by, "_", top_n, ".csv")
    write.csv(p$data, file.path(figures_folder,file_name_csv), row.names = TRUE) 
  }
}
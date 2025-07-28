# ==== libraries ====
library(ggplot2)
library(dplyr)    
library(purrr)     # For map functions
library(stringr)   
library(patchwork) 
library(ggridges)  
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

install.packages("ggnewscale") 
install.packages("UpSetR")
install.packages("ggbeeswarm")

library(ggnewscale)
library(UpSetR)  
library(ggbeeswarm) 



##########################
## FIGURE 2a
##########################

set.seed(1)
figures_folder <- "figures_14_06_2025"
supp_figures_folder <- "figures_14_06_2025/supp"
dir.create(figures_folder,recursive = TRUE)
dir.create(supp_figures_folder,recursive = TRUE)

# ==== clean up results from load_all_results.r ====
all_scores_list <- imap(results_preprocessed, ~{
  # For each dataset, iterate through its methods
  dataset_name <- .y
  methods_list <- .x

  # Map over the methods within the current dataset
  imap_dfr(methods_list, ~{
    method_name <- .y
    method_data <- .x

    if (!is.null(method_data) && nrow(method_data) > 0 && "prioritization_score" %in% names(method_data)) {
      # Select the score and add dataset/method identifiers
      method_data %>%
        select(prioritization_score) %>%
        mutate(
          method = method_name,
          dataset = dataset_name
        )
    } else {
      # Return an empty tibble/dataframe if data is missing or invalid
      tibble(prioritization_score = numeric(0), method = character(0), dataset = character(0))
    }
  })
})

# Combine the list of dataframes into one single dataframe
combined_scores_df <- bind_rows(all_scores_list)

# Define a consistent color scheme (adjust colors as needed)
method_colors <- c(
  "Decipher" = "#1f77b4",
  "NicheNet" = "#ff7f0e",
  "LIANA+" = "#2ca02c",
  "Connectome" = "#d62728",
  "NATMI" = "#9467bd"
)

# Optional: Remove datasets/methods with zero rows if any slipped through
combined_scores_df <- combined_scores_df %>% filter(!is.na(prioritization_score))

desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# Convert the 'method' column to a factor with the specified levels.
combined_scores_df <- combined_scores_df %>%
  mutate(method = factor(method, levels = rev(desired_method_order)))

# ==== Box-plot number of interactions ====
# Collect interaction counts into a data frame
interaction_counts <- do.call(rbind, lapply(names(results_preprocessed), function(ds) {
  df <- results_preprocessed[[ds]]
  data.frame(
    dataset = ds,
    method = names(df),
    interaction_count = sapply(df, nrow),
    stringsAsFactors = FALSE
  )
}))

# Ensure correct order and factor levels
interaction_counts$method <- factor(interaction_counts$method, levels = rev(desired_method_order))

# Define the breaks you want to keep (excluding the smallest and largest ones)
custom_breaks <- c(1e3, 1e4, 1e5)

# Plot
p <- ggplot(interaction_counts, aes(y = method, x = interaction_count, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.2, size = 2.8, alpha = 0.8) +  # changed from width to height since x/y are flipped
  scale_fill_manual(values = method_colors) +
  scale_x_log10(position = "bottom",breaks = custom_breaks) +  # place x-axis labels at the top
  labs(
    x = "Number of Interactions",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    # X-axis at the top, horizontal text
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 17, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    # Y-axis on left, horizontal (perpendicular to axis)
    axis.text.y = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    text = element_text(face = "bold")
    
  )

ggsave(file.path(figures_folder,"boxplot_number_of_interactions_by_method.png"), plot = p, width = 5, height = 4, dpi = 300)

write.csv(
  interaction_counts,
  file = file.path(figures_folder, "boxplot_number_of_interactions_by_method.csv"),
  row.names = TRUE
)

##########################
## FIGURE 2b
##########################

# Generate the Violin Plot
plot_violin_scores_by_method <- ggplot(
  combined_scores_df,
  aes(y = method,                     # flip: categorical axis on Y
      x = prioritization_score,       # numeric axis on X
      fill = method)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7,
    scale = "width",
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) +

  # same colour palette / order
  scale_fill_manual(values = method_colors, breaks = rev(desired_method_order)) +

  # numeric axis at the top, pseudo-log like in your box plot
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.1, base = 10),
    breaks = c(-10, -1, 0, 1, 10, 100),
    position = "bottom",
    name = "Prioritization Score"
  ) +

  # categorical axis order, but hide ticks and title
  scale_y_discrete(limits = rev(desired_method_order), name = NULL) +

  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    # X-axis (now across the top)
    axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 17, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    # Suppress Y-axis ticks/labels/title
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(face = "bold")
  ) +
  labs(fill = "Method")


# Save the Plot
output_filename_violin <- "fig_2b_violin_scores_by_method.png"
ggsave(file.path(figures_folder,output_filename_violin), plot = plot_violin_scores_by_method, width = 3.5, height = 4, dpi = 300)

write.csv(
  combined_scores_df,
  file = file.path(figures_folder, "fig_2b_violin_scores_by_method.csv"),
  row.names = TRUE
)


# FIGURE 2C

# FIGURE 2D


##########################
## FIGURE 2e
##########################
#load cytosig data
cytosig_significance   <- list() 
#cells_per_cluster <- list()  # New list to store cell counts per cluster

for (ds in names(datasets)) {
  dataset_path <- datasets[[ds]]
  pre_processing_filepath <- file.path(dataset_path, "pre_processing")
  cytosig_filepath <- file.path(dataset_path, "cytosig/0_outputs")
  reference_filepath <- "reference_data"
  #Cytosig results
  z_score_folder <- file.path(cytosig_filepath,"z_score/")
  p_value_folder <- file.path(cytosig_filepath,"p_value/")
  z_score_files <- list.files(z_score_folder)
  p_value_files <- list.files(p_value_folder)
  
  #seurat_object_oi <- readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))

  # Record number of cells per cluster
  #cluster_counts <- table(seurat_object_oi@meta.data$cluster)
  #cells_per_cluster[[ds]] <- as.data.frame(cluster_counts)

  # Load mapping table
  mapping_table <- read.csv(file.path(reference_filepath,"cytosig_mapping_table_ligands_genes.csv"),header=TRUE)

  # Process Cytosig significance
  cytosig_significance[[ds]] <- summarizeZScores(z_score_files, z_score_folder, mapping_table)
}

# Filter datasets with non-empty heatmaps
valid <- list()
all_ligands <- c()
all_cytosig_results <- list()
for (dataset in names(cytosig_significance)) {
  dataset_path <- datasets[[dataset]]
  decipher_filepath <- file.path(dataset_path, "data")
  df <- cytosig_significance[[dataset]]
  seurat_object_oi <- readRDS(file.path(decipher_filepath,"pseudobulk_seurat.rds"))
  # Keep only genes that exist in Seurat object
  valid_genes <- rownames(seurat_object_oi)
  valid_clusters <- unique(seurat_object_oi$cluster)
  df_unique <- df %>%
    filter(gene %in% valid_genes) %>%
    select(-gene) %>%
    distinct()  # Removes rows with identical ligand + all values

  # Keep only columns (i.e., cell types) that match valid clusters
  df_unique <- df_unique %>%
    select(ligand, all_of(valid_clusters))

  mat <- df_unique %>%
    tibble::column_to_rownames(var = "ligand") %>%
    as.matrix()
    
  row_mask <- apply(abs(mat), 1, max) > 2
  col_mask <- apply(abs(mat), 2, max) > 2
  mat_filt <- mat[row_mask, col_mask, drop = FALSE]
  if (nrow(mat_filt) > 0 && ncol(mat_filt) > 0) {
    valid[[dataset]] <- mat_filt
    all_cytosig_results[[dataset]] <- mat
    all_ligands <- union(all_ligands, rownames(mat_filt))  # keep building union

  }
}

for (dataset in names(valid)) {
  mat <- valid[[dataset]]

  # Create a completed matrix with all ligands as rows
  missing_ligands <- setdiff(all_ligands, rownames(mat))
  if (length(missing_ligands) > 0) {
    padding <- matrix(0, nrow = length(missing_ligands), ncol = ncol(mat),
                      dimnames = list(missing_ligands, colnames(mat)))
    mat <- rbind(mat, padding)
  }

  # Sort by ligand name to ensure row alignment across plots
  mat <- mat[sort(rownames(mat)), , drop = FALSE]
  valid[[dataset]] <- mat
}


#==== ROC curves ====
#results_for_comparison <- results_for_comparison[setdiff(names(results_for_comparison), "lupus")]
predictions_and_responses_all <- list()
auc_scores_by_datset <- list()
for (ds in names(datasets)) {
    dataset_path <- datasets[[ds]]
    pre_processing_filepath <- file.path(dataset_path, "pre_processing")
    meta_path <- "manuscript_analysis/data_for_meta_comparisons"
    output_figures_filepath <-  supp_figures_folder
    reference_filepath <- "reference_data"
    decipher_filepath <- file.path(dataset_path, "data")
    seurat_object_oi <- readRDS(file.path(decipher_filepath,"pseudobulk_seurat.rds"))
    
    L.set <- getForrestLRDatabase(file.path(reference_filepath,"connectomedb_forrest_lrc2p.csv"))
    L.set <- L.set %>%
    mutate(interaction = paste(ligand,receptor,sep="-"),
            lr = interaction) %>%
        unique()

    predictions_and_responses <- getPredictionsResponsesForMethods(
        results_for_comparison[[ds]],
        cytosig_significance[[ds]],
        L.set = L.set,
        seurat_object_oi,
        output_figures_filepath
    )

    predictions_and_responses_all[[ds]] <- predictions_and_responses


    all_predictions_across_methods <- predictions_and_responses$predictions
    all_responses_across_methods <- predictions_and_responses$responses

    AUC_scores <- plotROCAndExtractAUC(all_predictions_across_methods,all_responses_across_methods,output_figures_filepath,dataset_name = ds)
    saveRDS(AUC_scores,file.path(meta_path,paste(ds,"auc_scores.rds")))
    auc_scores_by_datset[[ds]] <- AUC_scores

    plotCytosigSignificanceMatrix(cytosig_significance[[ds]],output_figures_filepath)
}

#ggbeeswarm
results_df <- map_dfr(names(auc_scores_by_datset), function(dataset) {
  map_dfr(names(auc_scores_by_datset[[dataset]]), function(threshold) {
    map_dfr(names(auc_scores_by_datset[[dataset]][[threshold]]), function(method){
      vals <- auc_scores_by_datset[[dataset]][[threshold]][[method]]
          
          if (is.null(vals$auc) || is.null(vals$n_true)) return(NULL)

          data.frame(
            dataset = dataset,
            method = method,
            threshold = as.numeric(threshold),
            value = vals$auc,
            n_true = vals$n_true,
            stringsAsFactors = FALSE
          )
    })
  })
})


results_df <- results_df %>% filter(threshold == 2)  %>%
  mutate(flag = ifelse(n_true < 10, "*", ""))

# Reorder methods by median if desired
results_df$method <- reorder(results_df$method, results_df$value, FUN = median)

# Filter flagged points n_true lt 10
flagged_points <- filter(results_df, flag == "*")

# Plot
results_df$flagged <- results_df$flag == "*"

#color line by variance
# Calculate variance (or std dev) for each method
dataset_var <- results_df %>%
  group_by(dataset) %>%
  summarize(variance = var(value, na.rm = TRUE), .groups = "drop")

# Decide thresholds — top 3 lowest and 1-2 highest variance
dataset_var <- dataset_var %>%
  mutate(var_rank = rank(variance, ties.method = "first"))

# Color rules
metdataset_varhod_var <- dataset_var %>%
  mutate(line_color = case_when(
    var_rank <= 3 ~ "green",     # Most consistent
    var_rank >= (n() - 1) ~ "red",  # Most inconsistent
    TRUE ~ "gray"
  ))
results_df <- results_df %>%
  left_join(metdataset_varhod_var %>% select(dataset, line_color), by = "dataset")

# --- Prepare Data: Order Methods ---
print("Ordering methods on X-axis...")
available_methods_plot <- unique(results_df$method)
valid_order_plot <- intersect(desired_method_order, available_methods_plot)

# Filter data if some methods are not in the desired order list
if(length(valid_order_plot) < length(available_methods_plot)){
    warning("Some methods found in 'results_df$method' are not in 'desired_method_order'. Plotting only those included in the defined order and filtering data.")
    results_df <- results_df %>% filter(method %in% valid_order_plot)
}
# Stop if no data remains after filtering
if(nrow(results_df) == 0) {
    stop("No data remains after filtering 'results_df' based on 'desired_method_order'. Check method names.")
}

# Convert 'method' column to an ordered factor
results_df$method <- factor(results_df$method, levels = valid_order_plot)

# --- Define colors and desired order ---
color_k_value <- "#E69F00" # Gold/Orange
color_spearman <- "#56B4E9" # Light Blue

# --- Dynamically Create Line Color Mapping ---
# Identify the unique values in the 'line_color' column to map them
line_color_values_in_data <- sort(unique(results_df$line_color)) # Sort for consistent mapping order
if (length(line_color_values_in_data) == 0) {
    stop("No values found in the 'line_color' column.")
} else if (length(line_color_values_in_data) > 2) {
    warning(paste("Found more than 2 unique values in 'line_color':",
                  paste(line_color_values_in_data, collapse=", "),
                  "- Only the first two will be mapped to Gold/Blue."))
    # Map only the first two values found
    line_color_map <- setNames(c(color_k_value, color_spearman), line_color_values_in_data[1:2])
} else if (length(line_color_values_in_data) == 1) {
     warning(paste("Found only 1 unique value in 'line_color':", line_color_values_in_data[1],
                  "- Mapping it to Gold. All lines will have the same color."))
     line_color_map <- setNames(color_k_value, line_color_values_in_data[1])
} else {
   # Exactly two values found, map them
   line_color_map <- setNames(c(color_k_value, color_spearman), line_color_values_in_data)
}
print("Mapping for line colors:")
print(line_color_map)


# --- Generate the Plot with updated colors and order ---
print("Generating plot...")



####################
# FIGURE 2e
####################
p_updated <- ggplot(results_df, aes(x = method, y = value)) +

  # 1. Lines - Mapped to 'line_color', using first color scale
  #geom_line(aes(group = dataset, color = line_color), size = 1, alpha = 0.6) +
  geom_line(aes(group = dataset, color = "lightgray"), size = 1, alpha = 0.6) +
  #scale_color_manual(
  #  # name = "Line Group", # Optional legend name
  #  values = line_color_map, # Use the map created above
  #  guide = "none" # Hide legend for lines
  #) +

  # *** Introduce a new scale for color ***
  #new_scale_color() +

  # 2. Boxplot - Neutral colors (Plot after lines, before points?)
  geom_boxplot(outlier.shape = NA, width = 0.25, alpha = 0.4, color = "black", fill = "lightgray") +

  # 3. Beeswarm Points - Mapped to 'method', using second color scale
  geom_beeswarm(
    aes(color = method),   # Color points by method
    size = 4.0,            # Adjust point size if needed
    cex = 5,               # Use 'size', not 'cex' for ggplot point size control
    priority = "density",
    groupOnX = TRUE,        # Ensure beeswarm groups correctly on discrete x-axis
    alpha = 0.5
  ) +
  # Second color scale specifically for the beeswarm points
  scale_color_manual(
    # name = "Method",     # Optional legend name
    values = method_colors, # Use the method colors for points
    guide = "none"         # Hide legend for points
  ) +

  # Median line (same as before)
  stat_summary(fun = median, geom = "segment",
               aes(xend = after_stat(x), yend = after_stat(y)),
               size = 2.5, color = "black") +

  # Horizontal line at 0.5 (same as before)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +

  # Ensure X-axis order respects the factor levels
  scale_x_discrete(limits = levels(results_df$method), # Use levels from the factor
                   drop = FALSE) + # Prevent dropping unused levels

  # Labels and Theme
  labs(y = "AUROC target prediction", x = NULL) + # x = NULL removes x-axis title
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=17,face="bold"),
        panel.grid = element_blank()) # Rotate labels

# --- Save the Plot ---
output_filename_updated <- file.path(figures_folder,"variance_w_boxplot_beeswarm_auc_plot_updated.png")
ggsave(output_filename_updated, plot = p_updated, width = 4, height = 8, dpi = 300) # Adjusted width slightly
#THIS IS THE ONE!

write.csv(
  results_df,
  file.path(figures_folder, "variance_w_boxplot_beeswarm_raw_data.csv"),
  row.names = TRUE
)

####################
# FIGURE 2f
####################
# Visualise the results of 100 runs of Decipher without setting a seed

# Load libraries
library(dplyr)
library(purrr)
library(ggplot2)

# Initialize empty lists to store the results
all_runs_dat_1 <- list()
all_runs_dat_2 <- list()

N <- 100
# Loop through all 100 runs
for (i in 1:N) {
  # Load the data
  dat_1 <- readRDS(paste0("sample_analysis/validity/data/Run_", i, "_decipher_scores_by_cluster.rds"))
  dat_2 <- readRDS(paste0("sample_analysis/validity/data/Run_", i, "_decipher_scores_by_regulon_and_cluster.rds"))

  # Add the run number and cell type columns, and combine the results
  combined_tibble_1 <- bind_rows(
    lapply(names(dat_1), function(cell_type) {
      dat_1[[cell_type]] %>%
        mutate(Cells = cell_type,
               run_number = i)
    })
  )

  combined_tibble_2 <- bind_rows(
    lapply(names(dat_2), function(cell_type) {
      dat_2[[cell_type]] %>%
        mutate(Cells = cell_type,
               run_number = i)
    })
  )

  # Store the results in the lists
  all_runs_dat_1[[i]] <- combined_tibble_1
  all_runs_dat_2[[i]] <- combined_tibble_2
}

# Combine all runs into single tibbles
final_combined_tibble_1 <- bind_rows(all_runs_dat_1)
final_combined_tibble_2 <- bind_rows(all_runs_dat_2)

#===================================================
# Step 2: Identify Top 10 Interactions Across Runs
#===================================================

# Identify top 10 interactions for each run and cell type
top_10_interactions <- final_combined_tibble_1 %>%
  group_by(Cells, run_number) %>%
  arrange(desc(decipher_score)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Count how often each interaction appears in the top 10
interaction_consistency <- top_10_interactions %>%
  group_by(Cells, interaction) %>%
  summarise(count = n(), .groups = "drop")

# Calculate the proportion of runs each interaction was in the top 10
interaction_consistency <- interaction_consistency %>%
  mutate(proportion = count / 100)

#===================================================
# Step 3: Visualize the Consistency of Interactions
#===================================================
# Plot the proportion of runs where each interaction was in the top 10
p <- ggplot(interaction_consistency, aes(x = interaction, y = proportion, fill = Cells)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Proportion of Runs Where Interactions are in the Top 10",
       x = "Interaction",
       y = "Proportion of Runs (Top 10)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Calculate total count per interaction across all cell types
interaction_totals <- interaction_consistency %>%
  group_by(interaction) %>%
  summarise(total_count = sum(count), .groups = "drop")

# Join the total_count back into your interaction_consistency table
interaction_consistency <- interaction_consistency %>%
  left_join(interaction_totals, by = "interaction")

# Base R reordering 
interaction_consistency$interaction <- factor(
  interaction_consistency$interaction,
  levels = interaction_totals$interaction[order(interaction_totals$total_count,decreasing = FALSE)]
)

top20_interactions <- interaction_consistency %>%
  arrange(desc(proportion)) %>%
  slice_head(n = 30)

top20_mapped <- top20_interactions %>%
  mutate(Cells = case_when(
    grepl("^B_", Cells)                     ~ "B",
    grepl("^CD14_plus_Monocytes", Cells)    ~ "CD14+ Mono",
    grepl("^CD4_T", Cells)                  ~ "CD4 T",
    grepl("^CD8", Cells)                    ~ "CD8 T",
    TRUE                                    ~ Cells  # fallback if unmatched
  ))
  
p <- ggplot(top20_mapped, aes(x = Cells, y = interaction, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", name = "Count") +  # Viridis scale
  scale_y_discrete(position = "right") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 17,face="bold"),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 13), 
    panel.grid = element_blank(),
    legend.position = "bottom"     # move legend to bottom
  )

ggsave(file.path(figures_folder,"proportions_runs_top_10_flipped.png"), p, width = 4.2, height = 8)

# 1. Rename the Cells
interaction_consistency <- interaction_consistency %>%
  mutate(Cells = recode(Cells,
                        "CD8_T_cells" = "CD8 T",
                        "CD4_T_cells" = "CD4 T",
                        "CD14_plus_Monocytes" = "CD14+ Mono",
                        "B_cells" = "B"))

# 2. Plot
p <- ggplot(interaction_consistency, aes(x = interaction, y = Cells, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,   # remove y-axis label here
    fill = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 10),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# Save
ggsave(file.path(figures_folder,"proportions_runs_top_10_flipped_cleaned.png"), p, width = 10, height = 10)




#===================================================
# Step 4: Track Decipher Scores Over Runs
#===================================================
# Filter for interactions that are consistently in the top 10 across runs
consistent_interactions <- interaction_consistency %>%
  #filter(proportion == 1) %>%
  select(Cells, interaction)

# Join with the original data to track the decipher score changes
score_trends <- final_combined_tibble_1 %>%
  semi_join(consistent_interactions, by = c("Cells", "interaction")) %>%
  group_by(Cells, interaction, run_number) %>%
  summarise(mean_decipher_score = mean(decipher_score), .groups = "drop")

# Plot the trends
ggplot(score_trends, aes(x = run_number, y = mean_decipher_score, color = interaction, group = interaction)) +
  geom_line() +
  facet_wrap(~Cells, scales = "free_y") +
  theme_minimal() +
  labs(title = "Decipher Score Trends for Consistent Top 10 Interactions",
       x = "Run Number",
       y = "Mean Decipher Score")



####################
# FIGURE 2g
####################
#TODO: fix the caption letter here
plotDecipherPrioritizedMap("results/lupus",top_n=5,dataset_name="lupus", abs_decipher_plot_limit = 20,width=21,height=9)

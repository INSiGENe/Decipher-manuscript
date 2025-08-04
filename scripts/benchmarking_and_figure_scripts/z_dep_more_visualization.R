
# --- Load or Ensure 'combined_intersection_df' is available ---
# This dataframe should have columns: Intersection_Name, Count, Dataset, Degree
# It MUST contain data for all 31 intersection types before median calculation.

# Option 1: If it's already in your environment from the previous step (before filtering)
if (!exists("combined_intersection_df")) {
  # Option 2: Load it from the RDS file created earlier
  rds_file_unfiltered <- "figures/overlap_all_intersections_data.rds" # Make sure this is the UNFILTERED one
  if (file.exists(rds_file_unfiltered)) {
    print(paste("Loading unfiltered data from:", rds_file_unfiltered))
    combined_intersection_df <- readRDS(rds_file_unfiltered)
    # Ensure Degree column exists if it wasn't saved before
    if (!"Degree" %in% names(combined_intersection_df)) {
       combined_intersection_df <- combined_intersection_df %>%
        mutate(Degree = sapply(strsplit(as.character(Intersection_Name), " & "), length))
    }
  } else {
    stop("Required dataframe 'combined_intersection_df' not found. ",
         "Please ensure it's loaded or calculated (containing all 31 intersection types).")
  }
} else {
    print("Using existing 'combined_intersection_df'. Ensure it was not filtered.")
    # Ensure Degree column exists if needed
    if (!"Degree" %in% names(combined_intersection_df)) {
       combined_intersection_df <- combined_intersection_df %>%
        mutate(Degree = sapply(strsplit(as.character(Intersection_Name), " & "), length))
    }
}

# --- 1. Calculate Median Counts for ALL intersection types ---
print("Calculating median counts for UpSetR input...")
median_counts_all <- combined_intersection_df %>%
  group_by(Intersection_Name) %>%
  summarise(median_count = median(Count, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(median_count >= 0) # Keep zeros, remove NAs

# --- 2. Prepare data for UpSetR's fromExpression ---
# UpSetR expects intersection names with single '&' separators
# and the value associated with each name is the size (median count here).

# Check we have 31 intersections - if not, UpSetR might behave unexpectedly
# or the calculation function might need debugging
if(nrow(median_counts_all) != 31) {
   warning(paste("Expected 31 intersection types, but found", nrow(median_counts_all),
                 "after calculating medians. Proceeding, but results might be incomplete.",
                 "Ensure 'calculate_all_intersections' ran correctly for all datasets."))
}

# Create the named vector
upset_input_vector <- median_counts_all$median_count
# Convert names from "A & B" format to "A&B" format
names(upset_input_vector) <- gsub(" & ", "&", median_counts_all$Intersection_Name)

# --- 3. Generate the UpSet plot ---
print("Generating UpSet plot based on median intersection sizes...")

# Use png() device for saving UpSet plots
png(file.path(figures_folder,"overlap_upset_median_plot_ordered.png"), width = 10, height = 6, units = "in", res = 300)

#V2 plot
upset(
  fromExpression(upset_input_vector),  # Use the named vector of median counts

  # --- Modifications ---
  sets = rev(desired_method_order),         # 1. Set the ORDER of sets on the left
  nsets = length(rev(desired_method_order)),# Ensure nsets matches the provided sets
  queries = list(
      list(query = intersects, params = list("Decipher"), color = method_colors["Decipher"], active = T),
      list(query = intersects, params = list("NicheNet"), color = method_colors["NicheNet"], active = T),
      list(query = intersects, params = list("LIANA+"), color = method_colors["LIANA+"], active = T),
      list(query = intersects, params = list("Connectome"), color = method_colors["Connectome"], active = T),
      list(query = intersects, params = list("NATMI"), color = method_colors["NATMI"], active = T),
      list(query = intersects, params = list("Decipher", "LIANA+"), color = method_colors["Decipher"], active = T),
      list(query = intersects, params = list("Connectome", "Decipher"), color = method_colors["Decipher"], active = T),
      list(query = intersects, params = list("NicheNet", "Decipher"), color = method_colors["Decipher"], active = T),
      list(query = intersects, params = list("NicheNet", "Decipher","LIANA+"), color = method_colors["Decipher"], active = T),
      list(query = intersects, params = list("Decipher", "Connectome","NATMI"), color = method_colors["Decipher"], active = T)
   ),            # 2. Apply the coloring rules

  # --- Keep other parameters ---
  order.by = "freq",       # Order intersection bars by frequency (median counts)
  keep.order = TRUE,            
  decreasing = TRUE,                   # Show highest bars first
  mainbar.y.label = "Median Intersection Size", # Y-axis label for the main bar plot
  #sets.x.label = "Total Interactions in Top 100", # X-axis label for the set size plot
  point.size = 2.8,                    # Size of points in the matrix
  line.size = 1,                       # Size of lines in the matrix
  mb.ratio = c(0.6, 0.4),              # Ratio main bar height to matrix height
  text.scale = c(intersection_size=1.5, # Adjust text sizes
                 tick_labels=1.5,
                 set_size=1.5,
                 main_bar_text=1.5,
                 sets_names=1.5
                ),
  set_size.show = FALSE
  # min_size = 1, # Optional: uncomment if you only want combinations with median >= 1
  # show.numbers = FALSE, # Optional: uncomment to hide numbers above bars
)

# Close the PNG device
dev.off()

##########################
## end FIGURE 2? remove
##########################





# Create the final dual-axis boxplot
p_final <- ggplot(plot_data, aes(x = Method2)) +

  # Boxplots (same as before)
  geom_boxplot(aes(y = k_value, fill = "Search Space (k-value)"),
               color = "black", outlier.shape = NA, alpha = 0.7) +
  geom_boxplot(aes(y = Spearman * max(plot_data$k_value, na.rm = TRUE), fill = "Spearman Correlation"),
               color = "black", outlier.shape = NA, alpha = 0.7) +

  # --- 4. & 5. Add Horizontal Lines ---
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 100, color = color_k_value, linetype = "dashed", size = 0.6) + # Using the k-value color (gold/orange)

  # Manual color scale for fill
  scale_fill_manual(
      name = NULL,
      values = c("Search Space (k-value)" = color_k_value, "Spearman Correlation" = color_spearman)
      ) +

  # Dual Y-axis setup
  scale_y_continuous(
    name = "Search Space (k-value)",
    expand = expansion(mult = c(0.05, 0.05)), # Keep some expansion
    sec.axis = sec_axis(
        ~ . / max(plot_data$k_value, na.rm = TRUE),
        name = "Spearman Correlation"
        )
  ) +

  # X-axis - order now controlled by factor levels from plot_data$Method2
  scale_x_discrete(
      name = "Method Compared To",
      drop = TRUE # Prevent dropping levels if a method is missing in a specific facet
      ) +

  # Labels
  labs(title = "Comparison of Search Space & Spearman Correlation Across Method Pairs") +

  # Faceting by Method1 (order controlled by factor levels from plot_data$Method1)
  facet_grid(. ~ Method1, scales = "free_x", space = "free_x") +

  # Theme modifications
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", size = 0.4, linetype = "dotted"),
    panel.grid.minor.x = element_blank(),

    # --- 1. Increase Facet Separation ---
    panel.spacing.x = unit(1.5, "lines"), # Increased space BETWEEN facets

    strip.text.x = element_text(face = "bold", size=rel(1.1)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),

    # Color Y Axes (same as before)
    axis.title.y.left = element_text(color = color_k_value, face = "bold", size=rel(1.0)),
    axis.text.y.left = element_text(color = color_k_value, face="bold"),
    axis.ticks.y.left = element_line(color = color_k_value),
    axis.line.y.left = element_line(color = color_k_value),
    axis.title.y.right = element_text(color = color_spearman, face = "bold", size=rel(1.0)),
    axis.text.y.right = element_text(color = color_spearman, face="bold"),
    axis.ticks.y.right = element_line(color = color_spearman),
    axis.line.y.right = element_line(color = color_spearman),

    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save the plot (adjust dimensions if needed)
ggsave(file.path(figures_folder,"fig_2d.png"), plot = p_final, width = 8, height = 7.5, dpi = 300)

#let's not do complete
plot_data <- combined_df %>%
  filter(Method1 != Method2) %>%
  filter(Method1 %in% desired_method_order & Method2 %in% desired_method_order) %>%
  rowwise() %>%
  mutate(
    Method_A = min(c(Method1, Method2)),
    Method_B = max(c(Method1, Method2))
  ) %>%
  ungroup() %>%
  mutate(Method_Pair = paste(Method_A, Method_B, sep = "_"))
# Compute median Spearman per unique pair
pair_order <- plot_data %>%
  group_by(Method_Pair) %>%
  summarise(median_spearman = median(Spearman, na.rm = TRUE)) %>%
  arrange(desc(median_spearman)) %>%
  pull(Method_Pair)

  # Convert to factors WITH the desired levels
plot_data$Method1 <- factor(plot_data$Method1, levels = desired_method_order)
plot_data$Method2 <- factor(plot_data$Method2, levels = desired_method_order)
# Convert Method_Pair to a factor ordered by median Spearman
plot_data$Method_Pair <- factor(plot_data$Method_Pair, levels = pair_order)

# Create the final dual-axis boxplot

p_final <- ggplot(plot_data, aes(x = Method_Pair)) +

  geom_boxplot(aes(y = k_value, fill = "Search Space (k-value)"),
               color = "black", outlier.shape = NA, alpha = 0.7) +
  geom_boxplot(aes(y = Spearman * max(plot_data$k_value, na.rm = TRUE),
                   fill = "Spearman Correlation"),
               color = "black", outlier.shape = NA, alpha = 0.7) +

  # --- 4. & 5. Add Horizontal Lines ---
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 100, color = color_k_value, linetype = "dashed", size = 0.6) + # Using the k-value color (gold/orange)

  # Manual color scale for fill
  scale_fill_manual(
      name = NULL,
      values = c("Search Space (k-value)" = color_k_value, "Spearman Correlation" = color_spearman)
      ) +

  # Dual Y-axis setup
  scale_y_continuous(
    name = "Search Space (k-value)",
    expand = expansion(mult = c(0.05, 0.05)), # Keep some expansion
    sec.axis = sec_axis(
        ~ . / max(plot_data$k_value, na.rm = TRUE),
        name = "Spearman Correlation"
        )
  ) +

  # X-axis - order now controlled by factor levels from plot_data$Method2
  scale_x_discrete(name = NULL) +
  labs(title = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +

  # Theme modifications
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", size = 0.4, linetype = "dotted"),
    panel.grid.minor.x = element_blank(),

    # --- 1. Increase Facet Separation ---
    panel.spacing.x = unit(1.5, "lines"), # Increased space BETWEEN facets

    strip.text.x = element_text(face = "bold", size=rel(1.1)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),

    # Color Y Axes (same as before)
    axis.title.y.left = element_text(color = color_k_value, face = "bold", size=rel(1.0)),
    axis.text.y.left = element_text(color = color_k_value, face="bold"),
    axis.ticks.y.left = element_line(color = color_k_value),
    axis.line.y.left = element_line(color = color_k_value),
    axis.title.y.right = element_text(color = color_spearman, face = "bold", size=rel(1.0)),
    axis.text.y.right = element_text(color = color_spearman, face="bold"),
    axis.ticks.y.right = element_line(color = color_spearman),
    axis.line.y.right = element_line(color = color_spearman),

    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save the plot (adjust dimensions if needed)
ggsave(file.path(figures_folder,"combined_k_spearman_boxplot_final_ordered.png"), plot = p_final, width = 10, height = 7.5, dpi = 300)
#TODO: get data behind plot out
#where is the real matrix




#figure 2e removed
# this one had to with consistency and lines to mark that

p <- ggplot(results_df, aes(x = method, y = value)) +
  geom_line(aes(group = dataset, color = line_color), size = 1, alpha = 0.6) +
  
  geom_boxplot(outlier.shape = NA, width = 0.25, alpha = 0.4, color = "black", fill = "lightgray") +
  
  geom_beeswarm(
    aes(color = interaction(method, flagged)),
    size = 3.5, 
    cex = 5,
    priority = "density", 
    groupOnX = TRUE
  ) +
  
  scale_color_manual(
    values = c("green" = "green", "red" = "red", "gray" = "gray"),
    guide = "none"
  ) +
  
  stat_summary(fun = median, geom = "segment", 
               aes(xend = after_stat(x), yend = after_stat(y)), 
               size = 2.5, color = "black") +
  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  
  labs(y = "AUROC target prediction", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(figures_folder,"fig_2e.png"), plot = p, width = 4, height = 6, dpi = 300)



#################
#Remove this one
#################
#another type
p <- ggplot(results_df, aes(x = method, y = value)) +
  geom_line(aes(group = dataset), color = "gray", size = 1, alpha = 0.3) +
  
  # Use geom_beeswarm instead of geom_point
  geom_beeswarm(
    aes(color = interaction(method, flagged)),
    size = 3.5, 
    cex = 5,
    priority = "density", 
    groupOnX = TRUE
  ) +
  
  scale_color_manual(values = create_flag_color_scale(unique(results_df$method))) +
  
  stat_summary(fun = median, geom = "segment", 
               aes(xend = after_stat(x), yend = after_stat(y)), 
               size = 2.5, color = "black") +
  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  
  labs(y = "AUROC target prediction", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path(figures_folder,"beeswarm_auc_plot.png"), plot = p, width = 4, height = 6, dpi = 300)




######################
# start remove this 123
#####################
# ==== Decipher heatmap ====
#  # Not strictly required if using base R alternatives below

# --- Ensure 'results_for_comparison' list exists ---
if (!exists("results_for_comparison") || length(results_for_comparison) == 0) {
    stop("The 'results_for_comparison' list is empty or does not exist. Please run the data loading script.")
}
# Ensure the list of datasets is available
if (!exists("datasets")) {
    stop("The 'datasets' list mapping names to paths is required.")
}

# --- 1. Aggregate all Decipher scores across datasets/clusters ---
print("Aggregating Decipher scores...")
all_decipher_long <- list()
valid_datasets <- names(results_for_comparison) # Use datasets actually loaded

for (ds in valid_datasets) {
    if (!is.null(results_for_comparison[[ds]][["Decipher"]]) && length(results_for_comparison[[ds]][["Decipher"]]) > 0) {
        decipher_data_list <- results_for_comparison[[ds]][["Decipher"]]
        
        # Iterate through clusters within the dataset
        for (cluster_name in names(decipher_data_list)) {
            cluster_df <- decipher_data_list[[cluster_name]]
            
            # Basic check for expected columns
            if (!is.null(cluster_df) && nrow(cluster_df) > 0 && all(c("interaction", "score") %in% names(cluster_df))) {
                 # Ensure interaction is character
                cluster_df$interaction <- as.character(cluster_df$interaction)
                
                # Create a unique identifier for the column
                dataset_cluster_id <- paste(ds, cluster_name, sep = "_")
                
                # Select relevant columns and add identifiers
                temp_df <- data.frame(
                    dataset = ds,
                    cluster = cluster_name,
                    dataset_cluster = dataset_cluster_id,
                    interaction = cluster_df$interaction,
                    score = cluster_df$score
                )
                all_decipher_long[[length(all_decipher_long) + 1]] <- temp_df
            } else {
                 warning(paste("Skipping cluster", cluster_name, "in dataset", ds, "- data missing or invalid format."))
            }
        }
    } else {
         warning(paste("No Decipher results found or list is empty for dataset:", ds))
    }
}

# Combine all data frames into one
if (length(all_decipher_long) > 0) {
    combined_decipher_df <- do.call(rbind, all_decipher_long)
} else {
    stop("No valid Decipher data could be aggregated from 'results_for_comparison'.")
}

# --- 2. Identify Top 20 Interactions Overall ---
print("Identifying top 20 interactions...")
# Calculate the maximum score for each interaction across all dataset/clusters
# Using aggregate (base R)
interaction_max_scores <- aggregate(score ~ interaction, data = combined_decipher_df, FUN = max, na.rm = TRUE)

# Using dplyr (if available and preferred)
# 
# interaction_max_scores <- combined_decipher_df %>%
#   group_by(interaction) %>%
#   summarise(max_score = max(score, na.rm = TRUE), .groups = 'drop') %>%
#   arrange(desc(max_score))

# Sort by score descending and get top 20
interaction_max_scores <- interaction_max_scores[order(interaction_max_scores$score, decreasing = TRUE), ]
top_20_interactions <- head(interaction_max_scores$interaction, 20)

if (length(top_20_interactions) == 0) {
    stop("Could not determine top 20 interactions. Check aggregated data.")
}
print(paste("Top 20 interactions identified:", paste(top_20_interactions, collapse=", ")))

# --- 3. Prepare Data for Plotting ---
print("Preparing data for heatmap...")

# Filter the combined data frame for only the top 20 interactions
heatmap_data_long <- combined_decipher_df[combined_decipher_df$interaction %in% top_20_interactions, ]

# Create a complete grid of all top 20 interactions and all dataset_cluster combinations
# This ensures all cells are present for plotting, filling missing ones with 0
all_combinations <- expand.grid(
    interaction = top_20_interactions,
    dataset_cluster = unique(combined_decipher_df$dataset_cluster),
    stringsAsFactors = FALSE
)

# Merge with the actual scores, preserving dataset info
dataset_mapping <- unique(combined_decipher_df[, c("dataset_cluster", "dataset")])
heatmap_data_full <- merge(all_combinations, dataset_mapping, by = "dataset_cluster", all.x = TRUE)
heatmap_data_full <- merge(heatmap_data_full, heatmap_data_long[, c("interaction", "dataset_cluster", "score")],
                           by = c("interaction", "dataset_cluster"), all.x = TRUE)

# Replace NA scores (missing interactions for a specific dataset/cluster) with 0
heatmap_data_full$score[is.na(heatmap_data_full$score)] <- 0

# --- 4. Create Factor Levels for Ordering Axes ---
# Order interactions (Y-axis): Keep the top 20 order or sort alphabetically? Let's keep top 20 order.
y_axis_order <- top_20_interactions # Or sort(top_20_interactions) for alphabetical

# Order dataset_cluster (X-axis): Group by dataset first
x_axis_order <- unique(heatmap_data_full[order(heatmap_data_full$dataset), "dataset_cluster"])

# Apply factor levels
heatmap_data_full$interaction <- factor(heatmap_data_full$interaction, levels = rev(y_axis_order)) # Reverse for heatmap display
heatmap_data_full$dataset_cluster <- factor(heatmap_data_full$dataset_cluster, levels = x_axis_order)
heatmap_data_full$dataset <- factor(heatmap_data_full$dataset, levels = unique(heatmap_data_full[order(heatmap_data_full$dataset), "dataset"])) # Ensure dataset factor is ordered too


# --- 5. Create Plots (Dataset Strip and Heatmap) ---
print("Generating plot components...")

# Dataset strip data
dataset_strip_data <- unique(heatmap_data_full[, c("dataset_cluster", "dataset")])

# A. Dataset Strip Plot (Top Bar)
strip_plot <- ggplot(dataset_strip_data, aes(x = dataset_cluster, y = "Dataset", fill = dataset)) +
  geom_tile(width = 1) + # Color tiles fully
  # scale_fill_viridis_d(name = "Dataset") + # Or use another discrete scale if needed
  labs(fill = "Dataset") +
  theme_void() + # Minimal theme
  theme(
    legend.position = "bottom", # Legend below
    axis.text.x = element_blank(), # No text on strip plot x-axis
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 1, l = 10, unit = "pt") # Small margin, maybe adjust left margin if y-axis labels overlap later
  )

# B. Main Heatmap Plot
heatmap_plot <- ggplot(heatmap_data_full, aes(x = dataset_cluster, y = interaction, fill = score)) +
  geom_tile(color = "grey95", size = 0.1) + # Add faint lines between tiles
  # Choose a suitable color scale for scores (likely non-negative)
  scale_fill_viridis_c(name = "Decipher Score") + # Viridis 'plasma' option  option = "plasma"
  scale_fill_gradient2(
  low = "steelblue",
  mid = "white",
  high = "firebrick",  # "steelred" doesn't exist, but this is a good strong red
  midpoint = 0,        # Set this to the value you want as the center
  name = "Decipher Score"
  ) + # Alternative simple gradient
  labs(x = NULL, y = NULL) + # No axis titles
  theme_minimal(base_size = 10) + # Base theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.8)), # Smaller, rotated x labels
    axis.text.y = element_text(size = rel(0.9)), # Adjust y label size
    panel.grid = element_blank(), # No grid lines
    legend.position = "bottom", # Legend below
    legend.key.width = unit(1.5, "cm"), # Wider color bar key
    legend.key.height = unit(0.4, "cm") # Thinner color bar key
  )

# --- 6. Combine Plots using Patchwork ---
print("Combining plots...")
combined_decipher_heatmap <- strip_plot / heatmap_plot +
  plot_layout(heights = c(0.5, 10), # Adjust height ratio (more space for heatmap)
              guides = "collect") & # Collect legends at bottom
  theme(legend.position = "bottom",
        legend.box.margin = margin(t = 10) # Add margin above the collected legend
       )

# Print the combined plot
print(combined_decipher_heatmap)

# --- 7. Save Plot ---
print("Saving combined heatmap...")
output_filename <- file.path(figures_folder,"decipher_top20_heatmap.png")
dir.create("figures", showWarnings = FALSE) # Create directory if it doesn't exist

save_result <- tryCatch({
    # Adjust width/height as needed for readability
    ggsave(output_filename, plot = combined_decipher_heatmap, width = 16, height = 8, dpi = 300)
    TRUE
}, error = function(e) {
    warning("Failed to save the Decipher heatmap: ", e$message, call. = FALSE)
    FALSE
})

if (save_result) {
    print(paste("Decipher heatmap saved to", output_filename))
} else {
    print("Decipher heatmap was generated but could not be saved automatically.")
}

######################
# end remove this 123
#####################



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
output_filename_updated <- file.path(figures_folder,"figure_2e.png")
ggsave(output_filename_updated, plot = p_updated, width = 4, height = 8, dpi = 300) # Adjusted width slightly

write.csv(
  results_df,
  file = file.path(figures_folder, "figure_2e.csv"),
  row.names = TRUE
)
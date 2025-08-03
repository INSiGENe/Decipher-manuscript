# ---- Libraries ----
library(devtools)
load_all()

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(scales)
library(RColorBrewer)
library(tibble)
library(patchwork)
install.packages("ggnewscale")
library(ggnewscale)
library(purrr)


# ---- reproducible seed ----
set.seed(123)
figures_folder <- "figures_03_08_2025"

# cell-type funnel ----

#params
logfc_threshold <- log2(1.2)

dataset_paths <- list(
  SevCOVID = "results/SevCOVID_Azimuthl2",
  MilCOVID = "results/MilCOVID_Azimuthl2"
)

comparisons <- list(
  SevCOVID = list(case = "Severe", control = "Healthy", up_col = "Severe_Up", down_col = "Severe_Down"),
  MilCOVID = list(case = "Moderate", control = "Healthy", up_col = "Mild_Up", down_col = "Mild_Down")
)

#run
output <- run_all_comparisons(dataset_paths, comparisons, logfc_threshold)
df_sev <- output$results$SevCOVID
df_mild <- output$results$MilCOVID
seurat_objs <- output$seurat_objects  # Retain these for downstream analysis

#celltype_props <- calculate_celltype_proportions(seurat_objs)
celltype_props_normalized <- calculate_condition_normalized_proportions(seurat_objs)

celltype_props_normalized <- celltype_props_normalized %>%
  filter(!(severity_group == "Healthy" & Dataset == "MilCOVID")) %>%
  select(-Dataset) %>%
  group_by(cluster) %>%
  mutate(Proportion = NormalizedCount / sum(NormalizedCount) * 100) %>%
  ungroup() %>%
  mutate(severity_group = ifelse(severity_group == "Moderate", "Mild", severity_group)) %>%
  mutate(severity_group = factor(severity_group, levels = c("Healthy", "Mild", "Severe"))) 


#used to be called "final_res"
deg_summary_by_cluster <- combine_deg_counts(df_sev, df_mild)

write.csv(deg_summary_by_cluster, file.path(figures_folder,"final_deg_table.csv"))

# Reshape
abundance_degs <- deg_summary_by_cluster %>%
  select(cluster, Severe_Up, Severe_Down, Mild_Up, Mild_Down) %>%
  pivot_longer(cols = -cluster, names_to = "Group", values_to = "Count") %>%
  mutate(
    Condition = case_when(
      grepl("Severe", Group) ~ "Severe",
      grepl("Mild", Group) ~ "Mild"
    ),
    Direction = case_when(
      grepl("Up", Group) ~ "Up",
      grepl("Down", Group) ~ "Down"
    ),
    ColorGroup = paste(Condition, Direction, sep = "_")
  ) %>%
  group_by(cluster) %>%
  mutate(Percent = 100 * Count / sum(Count)) %>%
  ungroup()

# Order clusters by total DEGs if you want
abundance_degs$cluster <- factor(
  abundance_degs$cluster,
  levels = deg_summary_by_cluster %>% arrange(desc(Total)) %>% pull(cluster)
)

# Set cluster order based on total Moderate (Mild) proportion
cluster_order <- abundance_degs %>%
  filter(Condition == "Mild") %>%
  group_by(cluster) %>%
  summarise(TotalMild = sum(Percent)) %>%
  arrange(TotalMild) %>%
  pull(cluster)

# Apply ordering to the factor
abundance_degs <- abundance_degs %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

# Custom colors
custom_colors <- c(
  Severe_Up = "#B2182B",     # dark red
  Severe_Down = "#FDBBA0",   # light red
  Mild_Up = "#2166AC",       # dark blue
  Mild_Down = "#B2CDE3"      # light blue
)

# Ensure ordering matches shared_cluster_order
cluster_abundance <- celltype_props_normalized %>%
  group_by(cluster) %>%
  summarise(TotalNormalized = sum(NormalizedCount), .groups = "drop")
# Join abundance info
abundance_degs <- abundance_degs %>%
  left_join(cluster_abundance, by = "cluster") %>%
  mutate(
    Count = ifelse(TotalNormalized/3 < 0.01, 0, Count),  # set Count to 0 if under 1%,
    Percent = ifelse(TotalNormalized/3 < 0.01, 0, Percent),
    Count = ifelse(cluster == "Plasmablast", 0, Count)
  )

half_total_df <- abundance_degs %>%
  group_by(cluster) %>%
  summarise(half_total = sum(Count) / 2)

abundance_degs <- abundance_degs %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

shared_cluster_order <- levels(abundance_degs$cluster)

celltype_props_normalized <- celltype_props_normalized %>%
  mutate(cluster = factor(cluster, levels = shared_cluster_order))


# Plot
degs_Severe_moderate_plot <- ggplot(abundance_degs, aes(x = cluster, y = Count, fill = ColorGroup)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_point(data = half_total_df,
           aes(x = cluster, y = half_total),
           shape = 3,     # shape 3 = plus/cross
           size = 2,      # adjust size if needed
           color = "gray40",
           inherit.aes = FALSE)+
  #geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +  # <- Add this line
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(
    y = "Number of DEGs",
    x = "cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"  
    )+
  guides(fill = guide_legend(ncol = 1))

# Plot
celltype_props_normalized_plot <- ggplot(celltype_props_normalized, aes(x = Proportion, y = cluster, fill = severity_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Healthy" = "#A6CEE3",   # Light Blue
      "Mild" = "#FFFACD",      # Light Yellow (was Moderate)
      "Severe" = "#FFD700"     # Dark Yellow
    )
  ) +
  labs(x = "Proportion (%)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"  # 
  ) +
  guides(fill = guide_legend(ncol = 1))  # 

# Ensure ordering matches shared_cluster_order
text_data <- celltype_props_normalized %>%
  group_by(cluster) %>%
  summarise(TotalNormalized = sum(NormalizedCount), .groups = "drop") %>%
  mutate(TotalNormalized = 100*TotalNormalized/3) %>%
  mutate(
    cluster = factor(cluster, levels = shared_cluster_order),
    Label = sprintf("%.2f", TotalNormalized)
  )


text_plot <- ggplot(text_data, aes(y = cluster, x = 1, label = Label)) +
  geom_text(size = 4.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "Total cells (%)", y = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

#####################
#PCA
#####################
#parameters
conditions <- c(
  mild =  "MilCOVID_Azimuthl2", # Folder name for moderate data
  severe   =  "SevCOVID_Azimuthl2"    # Folder name for severe data
)

# Define the specific cell types (receiver cells) you want to analyze
selected_receiver_cells <- c( "Eryth", "CD16_Mono", "HSPC",    "CD4_TCM",  "Plasmablast",    "B_intermediate", "B_naive","CD8_Naive","NK","CD8_TEM","pDC","cDC2","Platelet","CD14_Mono","CD4_CTL" )

# Dynamically load data based on conditions defined above
regulon_deltas_list <- lapply(names(conditions), function(cond_name) {
  folder_name <- conditions[[cond_name]]
  file_path <- file.path("results", folder_name, "data/regulon_deltas_by_cluster.rds")
  cat("Loading data for:", cond_name, "from:", file_path, "\n")
  load_regulon_data(file_path, selected_receiver_cells)
})

# Set names of the main list (e.g., regulon_deltas_list$moderate, regulon_deltas_list$severe)
names(regulon_deltas_list) <- names(conditions)

# ==== Plot PCA of regulons ====
long_moderate <- get_long_deltas(regulon_deltas_list$mild, "Mild")
long_severe   <- get_long_deltas(regulon_deltas_list$severe, "Severe")

combined_long <- bind_rows(long_moderate, long_severe)
combined_long <- combined_long %>% rename("cluster"="Cluster")
# 2. Pivot to wide format: rows = cluster+Condition, cols = TFs
wide_deltas <- combined_long %>%
  select(cluster, Condition, name, deltaPagoda) %>%
  pivot_wider(names_from = name, values_from = deltaPagoda, values_fill = 0) %>%
  unite(cluster_Condition, cluster, Condition, sep = "_")

# 3. Run PCA
deltas_mat <- wide_deltas %>% column_to_rownames("cluster_Condition") %>% as.matrix()
#filter for informative features and run PCA
min_presence <- 0.8 * nrow(deltas_mat)
deltas_mat_filtered <- deltas_mat[, colSums(deltas_mat != 0) >= min_presence]
pca_res <- prcomp(deltas_mat_filtered, scale. = TRUE)

# 4. Extract scores and metadata
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$cluster_Condition <- rownames(pca_df)

cluster_condition_map <- combined_long %>%
  mutate(cluster_Condition = paste(cluster, Condition, sep = "_")) %>%
  distinct(cluster_Condition, cluster, Condition)

pca_df <- pca_df %>%
  left_join(cluster_condition_map, by = "cluster_Condition")

moderate_clusters <- pca_df %>% filter(Condition == "Mild") %>% pull(cluster)
severe_clusters   <- pca_df %>% filter(Condition == "Severe") %>% pull(cluster)

# clusters that appear only in Moderate
unique_moderate_clusters <- setdiff(moderate_clusters, severe_clusters)

pca_df <- pca_df %>%
  mutate(label_text = ifelse(
    Condition == "Severe" | cluster %in% unique_moderate_clusters,
    cluster,
    NA_character_
  ))

# 5. Plot
 # for pseudo log transformation

# Define custom colors
condition_colors <- c("Mild" = "#91C8F6", "Severe" = "#E4B731")  # light blue and gold

# Match Moderate and Severe by cluster
cluster_pairs <- intersect(moderate_clusters, severe_clusters)

segment_df <- map_dfr(cluster_pairs, function(clust) {
  pt1 <- pca_df %>% filter(cluster == clust, Condition == "Mild")
  pt2 <- pca_df %>% filter(cluster == clust, Condition == "Severe")
  
  dist <- sqrt((pt1$PC1 - pt2$PC1)^2 + (pt1$PC2 - pt2$PC2)^2)
  
  tibble(
    cluster = clust,
    x = pt1$PC1, y = pt1$PC2,
    xend = pt2$PC1, yend = pt2$PC2,
    dist = dist
  )
})

# Get loadings from PCA result
# Get variance explained
pca_var <- pca_res$sdev^2
pca_var_explained <- pca_var / sum(pca_var)

# View percentage of variance explained by PC1 and PC2
pca_var_explained[1:2]

#ok so let's do the loadings as part of the plot above
loadings <- as.data.frame(pca_res$rotation[, 1:2])  # Just PC1 and PC2 for now
loadings$TF <- rownames(loadings)

# Arrange by absolute contribution to PC1
top_PC1 <- loadings %>% arrange(desc(abs(PC1))) %>% head(10)

# Arrange by absolute contribution to PC2
top_PC2 <- loadings %>% arrange(desc(abs(PC2))) %>% head(10)


# ==== test 2 =====
# Define custom colors
condition_colors <- c("Mild" = "#91C8F6", "Severe" = "#E4B731")  # light blue and gold

# Calculate segments for connecting points
cluster_pairs <- intersect(
  pca_df %>% filter(Condition == "Mild") %>% pull(cluster),
  pca_df %>% filter(Condition == "Severe") %>% pull(cluster)
)

segment_df <- purrr::map_dfr(cluster_pairs, function(clust) {
  pt1 <- pca_df %>% filter(cluster == clust, Condition == "Mild")
  pt2 <- pca_df %>% filter(cluster == clust, Condition == "Severe")
  dist <- sqrt((pt1$PC1 - pt2$PC1)^2)
  tibble(
    cluster = clust,
    x = pt1$PC1, xend = pt2$PC1,
    y = pt1$cluster, yend = pt2$cluster,
    dist = dist
  )
})

pca_df <- pca_df %>%
  mutate(cluster = factor(cluster, levels = cleanSymbols(shared_cluster_order)))


# Plot
pca_embedding_plot <- ggplot(pca_df, aes(x = PC1, y = cluster)) +
  # Dashed line by distance
  geom_segment(data = segment_df,
               aes(x = x, xend = xend, y = y, yend = yend, color = dist),
               linetype = "dashed", size = 1) +
  scale_color_gradient(
    name = "Distance",
    low = "#DCE6F2",
    high = "#08306B",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                         barwidth = 4, barheight = 0.4)
  ) +
  ggnewscale::new_scale_color() +
  scale_y_discrete(limits = cleanSymbols(shared_cluster_order)) +

  # Plot points
  geom_point(aes(color = Condition), size = 4) +
  scale_color_manual(values = condition_colors, name = "Condition",
                     guide = guide_legend(override.aes = list(size = 4),ncol = 1,title.position = "top",
    title.hjust = 0.5)) +

  # Axis transformation and theme tweaks
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    #axis.text = element_text(face = "bold"),
    #axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  ) +
  labs(x = "PC1", y = "cluster")

# ok last one
# 1. Prepare diff_df using cluster_condition_map
diff_df <- deltas_mat_filtered %>%
  as.data.frame() %>%
  rownames_to_column("cluster_Condition") %>%
  left_join(cluster_condition_map, by = "cluster_Condition") %>%
  pivot_longer(
    cols = -c(cluster_Condition, cluster, Condition),
    names_to = "TF",
    values_to = "Value"
  ) %>% select(!cluster_Condition) %>%
  pivot_wider(
    names_from   = Condition,
    values_from  = Value,
    values_fn    = mean
  ) %>%
  mutate(
    Diff = Severe - Mild
  ) %>%
  filter(TF %in% top_PC1$TF) %>%
  mutate(
    TF      = factor(TF, levels = top_PC1$TF),
    cluster = factor(cluster, levels = cleanSymbols(shared_cluster_order))
  )

# Create all combinations of TF and cluster
full_grid <- expand_grid(
  TF      = levels(diff_df$TF),
  cluster = levels(diff_df$cluster)
)

# Join your existing data to the full grid
diff_df_complete <- full_grid %>%
  left_join(diff_df, by = c("TF", "cluster"))


# Explicitly relevel cluster with your desired order
diff_df_complete$cluster <- factor(diff_df_complete$cluster, levels = cleanSymbols(shared_cluster_order))


# Format PC1 variance as percentage
pc1_var_percent <- round(pca_var_explained[1] * 100, 1)

# Now plot
heatmap_plot <- ggplot(diff_df_complete, aes(x = TF, y = cluster, fill = Diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    name     = "Δ TF Activity",
    na.value = "grey90",  # Color for missing values
    guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 6,
    barheight = 0.4
  )
  ) +
  labs(x = paste0("PC1 loadings")) +  # Add dynamic x-axis title
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.y      = element_blank(),
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(2, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10)
  )

# 3. Combine with existing plots

final_plot <- (
     text_plot
  +  celltype_props_normalized_plot
  +  degs_Severe_moderate_plot
  +  pca_embedding_plot
  +  heatmap_plot
) +
  plot_layout(
    ncol   = 5,
    widths = c(0.12, 0.30, 0.30, 0.20, 0.20)
  ) &
  theme(legend.position = "bottom")

# 4. Save the composite figure
ggsave(
  file.path(figures_folder, "figure_4a_e.png"),
  final_plot,
  width  = 11,
  height = 8,
  dpi    = 300
)


write.csv(text_data,
          file.path(figures_folder, "figure_4a.csv"),
          row.names = FALSE)

write.csv(celltype_props_normalized,
          file.path(figures_folder, "figure_4b.csv"),
          row.names = FALSE)

write.csv(abundance_degs,
          file.path(figures_folder, "figure_4c.csv"),
          row.names = FALSE)

write.csv(pca_df,
          file.path(figures_folder, "figure_4d.csv"),
          row.names = FALSE)

write.csv(diff_df_complete,
          file.path(figures_folder, "figure_4e.csv"),
          row.names = FALSE)


# TF activity deltas
#didn't change moderate to mild here

# ==== TF activity deltas (Severe vs Mild) ====

# Define the specific cell types (receiver cells) you want to analyze
selected_receiver_cells <- c( "Eryth", "CD16_Mono", "HSPC",    "CD4_TCM",  "Plasmablast",    "B_intermediate", "B_naive","CD8_Naive","NK","CD8_TEM","pDC","cDC2","Platelet","CD14_Mono","CD4_CTL" )

# Number of top regulons (TFs) to display in each heatmap
top_n_regulons <- 10

# Number of cell types to combine per output PNG file
clusters_per_group_in_output <- 1 # Adjust as needed (e.g., 1, 2, 3)

# Define condition names and the subfolders where their data resides
# The names ('moderate', 'severe') will be used throughout the script
conditions <- c(
  mild =  "MilCOVID_Azimuthl2", # Folder name for moderate data
  severe   =  "SevCOVID_Azimuthl2"    # Folder name for severe data
)

# Dynamically load data based on conditions defined above
regulon_deltas_list <- lapply(names(conditions), function(cond_name) {
  folder_name <- conditions[[cond_name]]
  file_path <- file.path("results", folder_name, "data/regulon_deltas_by_cluster.rds")
  cat("Loading data for:", cond_name, "from:", file_path, "\n")
  load_regulon_data(file_path, selected_receiver_cells)
})

# Set names of the main list (e.g., regulon_deltas_list$moderate, regulon_deltas_list$severe)
names(regulon_deltas_list) <- names(conditions)

# Refine selected_receiver_cells to only those present in the loaded data (at least in the first condition)
# This prevents errors if a specified cell type wasn't found in the files.
initial_cell_count <- length(selected_receiver_cells)
if (length(regulon_deltas_list) > 0 && !is.null(regulon_deltas_list[[1]])) {
    selected_receiver_cells <- intersect(selected_receiver_cells, names(regulon_deltas_list[[1]]))
} else {
    warning("Could not validate selected_receiver_cells against loaded data. Proceeding with the original list.")
}
final_cell_count <- length(selected_receiver_cells)
cat("Initial selected cell types:", initial_cell_count, "\n")
cat("Validated selected cell types found in data:", final_cell_count, "\n")
if(final_cell_count == 0) {
    stop("No valid selected cell types found in the loaded data. Please check 'selected_receiver_cells' and data files.")
}

# Calculate Global Color Scale Limit
absolute_max <- find_absolute_max(regulon_deltas_list)
cat("Global absolute max deltaPagoda for scaling:", absolute_max, "\n")
# Ensure absolute_max is not zero or negative, set a minimum limit if needed
if(is.na(absolute_max) || absolute_max <= 0) {
    warning("Could not determine a valid absolute max deltaPagoda. Setting scale limit to 1.")
    absolute_max <- 1
}

# Execute plot generation
#generated_plots <- generate_sorted_plots(selected_receiver_cells, regulon_deltas_list, conditions, top_n_regulons, absolute_max)

# Create combined plots with titles (one combined plot per cell type)
#celltype_combined_plots <- create_combined_plots_per_celltype(generated_plots, selected_receiver_cells)

# Function to save combined plots in groups
#save_grouped_plots(
#    combined_plots = celltype_combined_plots,
#    clusters_per_group = clusters_per_group_in_output,
#    output_dir_base = ".", # Use the base path
#    output_folder_name = figures_folder
#)

single_TF_heatmap <- function(cell_type,
                              sort_by    = c("mild","severe"),
                              top_n      = 20,
                              deltas_list,
                              global_max) {
  sort_by   <- match.arg(sort_by)
  sort_cond <- if (sort_by == "mild") "Mild" else "Severe"

  # helper: extract a named vector from either a numeric vector or your 3-col df
  vec_from <- function(x) {
    if (is.data.frame(x)) {
      if (all(c("deltaPagoda","name") %in% names(x))) {
        v <- x$deltaPagoda
        names(v) <- x$name
        return(v)
      } else {
        stop("Data.frame must have columns 'deltaPagoda' and 'name'.")
      }
    } else if (is.numeric(x) && !is.null(names(x))) {
      return(x)
    } else {
      stop("Each element must be a named numeric vector or a data.frame with 'deltaPagoda'+'name'.")
    }
  }

  # pull out moderate & severe vectors
  mod_vec <- vec_from(deltas_list$mild[[cell_type]])
  sev_vec <- vec_from(deltas_list$severe  [[cell_type]])

  # turn into tibbles
  mod_df <- enframe(mod_vec, name="TF", value="Mild")
  sev_df <- enframe(sev_vec, name="TF", value="Severe")

  # join & pivot
  df_long <- full_join(mod_df, sev_df, by="TF") %>%
    pivot_longer(c(Mild,Severe),
                 names_to="Condition",
                 values_to="Delta")

  # pick top N by |Delta| in sort_cond
  top_tfs <- df_long %>%
    filter(Condition == sort_cond, !is.na(Delta)) %>%
    arrange(desc(abs(Delta))) %>%
    slice_head(n = top_n) %>%
    arrange(Delta) %>%
    pull(TF)

  df_top <- df_long %>%
    filter(TF %in% top_tfs) %>%
    mutate(
      TF        = factor(TF, levels = top_tfs),
      Condition = factor(Condition, levels = c("Mild","Severe"))
    )

  # symmetric color limits
  lim <- max(abs(df_top$Delta), na.rm = TRUE)

  # plot
  ggplot(df_top, aes(TF, Condition, fill = Delta)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      limits   = c(-global_max, global_max),
      na.value = "grey80"
    ) +
    labs(title = paste0(cell_type, " (sorted by ",sort_by,")"), x = "TF", y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title       = element_text(hjust = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid       = element_blank(),
      legend.position = "bottom"
    )
}


absolute_max <- find_absolute_max(regulon_deltas_list)

for (cell_type in c("CD16_Mono","CD14_Mono")){
  for(sort_by in c("severe","mild")){
    top_n = 20
    p <- single_TF_heatmap(
      cell_type   = cell_type,
      sort_by     = sort_by,
      top_n       = top_n,
      deltas_list = regulon_deltas_list,
      global_max = absolute_max
    )

    file_name <- paste0("figure_4f_",cell_type,"_sorted_by_", sort_by, "_", top_n, ".png")
    ggsave(file.path(figures_folder,file_name), p, width = 8, height = 3)
    file_name_csv <- paste0("figure_4f_",cell_type,"_sorted_by_", sort_by, "_", top_n, ".csv")
    write.csv(p$data, file.path(figures_folder,file_name_csv), row.names = TRUE) 

  }
}
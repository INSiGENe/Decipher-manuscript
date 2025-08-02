
library(igraph)
library(dplyr)
install.packages("rentrez")
library(rentrez)
library(data.table)
library(ggplot2)

#functions ----

# Initialize an empty data frame for edges and combined log2FC data
initialize_data_frames <- function() {
  list(
    edges = data.frame(),
    combined_data = data.frame()
  )
}

# Extract target genes and average log fold change (avg_log2FC) for each regulon
extract_regulon_data <- function(regulons, capped_regulons_all_clusters, significant_regulon_markers_by_cluster,selected_ct) {
  data <- initialize_data_frames()

  for (gene in regulons) {
    regulon_data <- significant_regulon_markers_by_cluster[[selected_ct]][[gene]]
    em_data <- capped_regulons_all_clusters[[selected_ct]] %>% filter(source == gene)

    data$edges <- rbind(data$edges, data.frame(from = gene, to = em_data$target))
    data$combined_data <- rbind(data$combined_data, regulon_data[, c("regulon", "tg_gene", "avg_log2FC")])
  }

  data
}

# Create color mapping for log2FC values
generate_log2fc_colors <- function(log2fc_values, num_colors = 100) {
  valid_log2fc <- log2fc_values[log2fc_values != -999 & !is.na(log2fc_values)]
  max_abs_log2fc <- max(abs(valid_log2fc), na.rm = TRUE)
  breaks <- seq(-max_abs_log2fc, max_abs_log2fc, length.out = num_colors + 1)
  color_palette <- colorRampPalette(c("cornflowerblue", "white", "coral1"))(num_colors)

  colors <- cut(log2fc_values, breaks = breaks, labels = color_palette, include.lowest = TRUE)
  colors <- as.character(colors)
  colors[log2fc_values == -999 | is.na(log2fc_values)] <- "white"

  colors
}

# Set vertex attributes (color, size, label size) for the graph
set_vertex_attributes <- function(g, log2fc_colors, regulons) {
  vertex_colors <- rep("white", vcount(g))
  vertex_colors[V(g)$name %in% names(log2fc_colors)] <- log2fc_colors[V(g)$name %in% names(log2fc_colors)]
  vertex_colors[V(g)$name %in% regulons] <- "darkgoldenrod1"

  vertex_size <- ifelse(V(g)$name %in% regulons, 10, 5)
  vertex_label_cex <- ifelse(V(g)$name %in% regulons, 1.0, 0.6)

  list(colors = vertex_colors, sizes = vertex_size, label_cex = vertex_label_cex)
}

# Main function to create and plot the graph
plot_regulon_graph <- function(data, regulons) {
  g <- graph_from_data_frame(data$edges, directed = FALSE)
  log2fc_values <- setNames(data$combined_data$avg_log2FC, data$combined_data$tg_gene)
  log2fc_colors <- generate_log2fc_colors(log2fc_values)

  vertex_attrs <- set_vertex_attributes(g, log2fc_colors, regulons)
}

get_pubmed_count <- function(gene) {
    search_result <- entrez_search(db = "pubmed", term = gene)
    return(search_result$count)
  }


#parameters ----

# Read this https://www.nature.com/articles/srep16923
set.seed(123)
figures_folder <- "figures_01_08_2025"
comparison_name <- "SevCOVID_Azimuthl2"
data_path <- file.path("results",comparison_name,"data")
capped_regulons_all_clusters <- readRDS(file.path(data_path,"capped_regulons_all_clusters.rds"))
regulon_deltas_by_cluster <- readRDS(file.path(data_path,"regulon_deltas_by_cluster.rds"))
significant_regulon_markers_by_cluster <- readRDS(file.path(data_path,"significant_regulon_markers_by_cluster.rds"))
n_pubmed <- 40
#scale_y_discrete(position = "right") #do this for severe and comment out for moderate

for(selected_ct in c("CD14_Mono","CD16_Mono")){
  # Run the refactored code
  regulons <- regulon_deltas_by_cluster[[selected_ct]] %>%
    filter(class == "real") %>%
    arrange(desc(deltaPagoda)) %>%
    dplyr::slice_head(n=10) %>%
    pull(name)
  data <- extract_regulon_data(regulons, capped_regulons_all_clusters, significant_regulon_markers_by_cluster, selected_ct)
  plot_regulon_graph(data, regulons)

  # Create a binary matrix
  binary_matrix <- data$edges %>%
    mutate(value = 1) %>%  # Assign 1 for each edge
    tidyr::pivot_wider(names_from = to, values_from = value, values_fill = 0) %>%
    tibble::column_to_rownames("from")


  # Filter columns where more than one 'from' talks to a 'to'
  binary_matrix <- binary_matrix[, colSums(binary_matrix) > 1]

  #see frequency of genes mentioned in pubmed
  # Define genes
  genes <- colnames(binary_matrix)  # Replace with your genes of interest

  # # Apply the function to each gene and store the results in a data frame
  PubMed_Count = sapply(genes, get_pubmed_count)
  gene_counts <- data.frame(
    Gene = genes,
    PubMed_Count = PubMed_Count
  )

  top_genes <- gene_counts %>%
    filter(Gene %in% genes) %>%
    arrange(desc(PubMed_Count)) %>%
    top_n(n = n_pubmed, wt = PubMed_Count) %>%
    pull(Gene)

  binary_matrix <- binary_matrix[,top_genes]

  # Convert binary_matrix to matrix format
  binary_matrix_mat <- as.matrix(binary_matrix)
  # Define the color palette with custom logic for red and white
  for(i in 1:nrow(binary_matrix_mat)){
    for(j in 1:ncol(binary_matrix_mat)){
      if(binary_matrix_mat[i,j] == 1 & colnames(binary_matrix_mat)[j] %in% top_genes)
        binary_matrix_mat[i,j] <- 1.2
    }
  }
  my_palette <- c("white", "black", "red")
  file_name <- paste(selected_ct,"tg_heatmaps.png")
  file_dir <- file.path("figures")

  # Perform hierarchical clustering on rows and columns
  row_dendrogram <- hclust(dist(binary_matrix_mat))   # Hierarchical clustering on rows
  col_dendrogram <- hclust(dist(t(binary_matrix_mat))) # Hierarchical clustering on columns

  # Order the matrix based on clustering
  binary_matrix_ordered <- binary_matrix_mat[row_dendrogram$order, col_dendrogram$order]

  # Convert the ordered matrix to a data frame for ggplot2
  binary_matrix_df <- as.data.frame(binary_matrix_ordered)
  binary_matrix_df$row <- factor(rownames(binary_matrix_ordered), levels = rownames(binary_matrix_ordered)) # Set ordered row levels

  # Melt the data for ggplot2
  binary_matrix_melted <- melt(binary_matrix_df, id.vars = "row")

  # Create the heatmap plot
  p <- ggplot(binary_matrix_melted, aes(x = variable, y = row, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = my_palette,
                         breaks = c(-0.1, 0.9, 1.1, 1.2),
                         limits = c(-0.1, 1.2),
                         guide = "none") +
    labs(title = NULL) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      # Conditionally set y-axis labels
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    ) +
    scale_y_discrete(position = "right") #do this for severe and comment out for moderate
    print(comparison_name)
    print(selected_ct)
    #scale_y_discrete(labels = function(x) ifelse(x %in% top_genes, x, ""))  # Conditional row labels
  filename <- paste0(comparison_name,"_",selected_ct,"_tgs_heatmap_pubmed.png")
  ggsave(file.path(figures_folder,filename),p,width = 12,height = 5.5,units="cm")
}


#network plots

library(igraph)
library(scales)
library(dplyr)
library(tidyr)


# 1. Define Helper Functions

# Function to reformat cell type names for plotting
formatCellTypeNamesForPlotting <- function(cluster_names) {
  formatted <- gsub("_minus_", "-", cluster_names)
  formatted <- gsub("_plus_", "+", formatted)
  formatted <- gsub("_", " ", formatted)
  return(formatted)
}

# Function to shorten long strings (first 3 & last 3 characters)
get_first_and_last_three <- function(x) {
  if (nchar(x) <= 8) {
    return(x)  # If the string is too short, return as is
  } else {
    return(paste0(substr(x, 1, 4), ".", substr(x, nchar(x) - 3, nchar(x))))
  }
}

# Function 1: Check Data Inputs
check_cluster_exists <- function(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, regulon_deltas_by_cluster) {
  if (!cluster_name %in% names(regulon_deltas_by_cluster)) {
    stop(paste("Cluster", cluster_name, "not found in regulon_deltas_by_cluster"))
  }
  if (!cluster_name %in% names(decipher_scores_by_regulon_and_cluster)) {
    stop(paste("Cluster", cluster_name, "not found in decipher_scores_by_regulon_and_cluster"))
  }
  if (!cluster_name %in% names(decipher_scores)) {
    stop(paste("Cluster", cluster_name, "not found in decipher_scores"))
  }
}

# Function 2: Extract Top Regulons
get_top_regulons <- function(cluster_name, regulon_deltas_by_cluster, n = 10) {
  this_regulon_deltas <- regulon_deltas_by_cluster[[cluster_name]]
  top_regulons_df <- this_regulon_deltas %>%
    filter(class == "real") %>%
    arrange(desc(deltaPagoda)) %>%
    head(n)
  top_regulons <- top_regulons_df$name
  list(top_regulons_df = top_regulons_df, top_regulons = top_regulons)
}

# Function 3: Extract Top Interactions
get_top_interactions <- function(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, top_regulons,top_interactions = NULL) {
  this_decipher_scores_by_regulon <- decipher_scores_by_regulon_and_cluster[[cluster_name]]
  this_decipher_scores <- decipher_scores[[cluster_name]]
  if(is.null(top_interactions)){
      top_20_interactions <- this_decipher_scores %>%
      arrange(desc(decipher_score)) %>%
      slice_head(n = 20) %>%
      pull(interaction)
  } else {
    top_20_interactions <- top_interactions
  }

  top_interactions <- this_decipher_scores_by_regulon %>%
    filter(regulon %in% top_regulons) %>%
    filter(interaction %in% top_20_interactions) %>%
    group_by(regulon) %>%
    slice_max(order_by = imp.perm, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    drop_na(interaction)
  top_interactions
}

# Function 4: Build Sender → Ligand Edges
build_sender_ligand_edges <- function(top_interactions, feature_statistics, sender_cts,global_sender_ligand_max) {
  normalized_fs <- feature_statistics %>%
    mutate(normalized.counts = sum.counts / n.cell) %>%
    group_by(condition, feature) %>%
    mutate(total.normalized.counts = sum(normalized.counts)) %>%
    ungroup() %>%
    mutate(frac.normalized.counts = normalized.counts / total.normalized.counts)
  
  sender_ligand <- normalized_fs %>%
    filter(condition == 'case', cluster %in% sender_cts, feature %in% top_interactions$ligand) %>%
    group_by(feature) %>%
    slice_max(order_by = frac.normalized.counts, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(total_counts = frac.normalized.counts) %>%
    select(ligand = feature, sender_cluster = cluster, total_counts) %>%
    #mutate(weight = rescale(total_counts, to = c(1, 5)))
    mutate(weight = 5 * (total_counts / global_sender_ligand_max))

  # Format sender cell type names
  sender_ligand$sender_cluster <- formatCellTypeNamesForPlotting(sender_ligand$sender_cluster)
  sender_ligand$sender_cluster <- sapply(sender_ligand$sender_cluster, get_first_and_last_three)
  
  # Build edge table
  edges_sender_ligand <- sender_ligand %>%
    select(from = sender_cluster, to = ligand, weight) %>%
    mutate(edge_type = "Sender_Ligand", colour = weight)
  
  list(sender_ligand = sender_ligand, edges_sender_ligand = edges_sender_ligand)
}

# Function 5: Build Ligand → Receptor Edges
build_ligand_receptor_edges <- function(this_decipher_scores, sender_ligand, top_interactions,global_decipher_score_max) {
  ligand_receptor_edges <- this_decipher_scores %>%
    filter(ligand %in% sender_ligand$ligand, receptor %in% top_interactions$receptor) %>%
    select(ligand, receptor, decipher_score) %>%
    rename(from = ligand, to = receptor, weight = decipher_score) %>%
    mutate(edge_type = "Ligand_Receptor", colour = weight) %>%
    mutate(weight = 5 * abs(weight) / global_decipher_score_max)
    #mutate(weight = rescale(abs(weight), to = c(1, 5)))
  ligand_receptor_edges
}

# Function 6: Build Receptor → TF (Regulon) Edges
build_receptor_tf_edges <- function(top_interactions,global_receptor_tf_col_max) {
  receptor_tf_edges <- top_interactions %>%
    select(receptor, regulon, imp.perm, spearman.cor) %>%
    distinct() %>%
    mutate(weight = imp.perm,
           edge_type = "Receptor_TF",
           colour = imp.perm * sign(spearman.cor)) %>%
    rename(from = receptor, to = regulon) %>%
    #mutate(weight = rescale(weight, to = c(1, 5)))
    #mutate(weight = 5 * imp.perm / max(imp.perm, na.rm = TRUE),
    mutate(weight = 7.5 * imp.perm / global_receptor_tf_col_max,
       colour = imp.perm * sign(spearman.cor))  # scale colour later
  receptor_tf_edges
}

# Function 7: Combine Edges and Create Nodes Data Frame
build_graph_components <- function(edges_sender_ligand, ligand_receptor_edges, receptor_tf_edges, top_interactions, sender_ligand, top_regulons_df,global_deltaPagoda_max) {
  
  all_edges <- bind_rows(
    edges_sender_ligand %>% select(from, to, weight, edge_type, colour),
    ligand_receptor_edges %>% select(from, to, weight, edge_type, colour),
    receptor_tf_edges %>% select(from, to, weight, edge_type, colour)
  )
  
  # Nodes from different sources
  nodes_sender <- sender_ligand$sender_cluster
  nodes_ligand <- sender_ligand$ligand
  nodes_receptor <- top_interactions$receptor
  nodes_tf <- top_interactions$regulon
  
  nodes <- data.frame(
    name = unique(c(nodes_sender, nodes_ligand, nodes_receptor, nodes_tf)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(layer = case_when(
      name %in% nodes_sender ~ "Sender Cell Type",
      name %in% nodes_ligand ~ "Ligand",
      name %in% nodes_receptor ~ "Receptor",
      name %in% nodes_tf ~ "TF",
      TRUE ~ "Other"
    ))
  
  # Merge deltaPagoda information for TF nodes if available
  nodes <- nodes %>%
    left_join(top_regulons_df %>% select(name, deltaPagoda), by = "name") %>%
    mutate(color = case_when(
      layer == "TF" ~ col_numeric(palette = c("white", "tomato"),
                                   #domain = c(0, max(top_regulons_df$deltaPagoda, na.rm = TRUE)))(deltaPagoda),
                                   domain = c(0, global_deltaPagoda_max))(deltaPagoda),
      layer == "Sender Cell Type" ~ "cadetblue1",
      layer == "Ligand" ~ "darkolivegreen2",
      layer == "Receptor" ~ "darkorange",
      TRUE ~ "grey"
    ))
  
  list(all_edges = all_edges, nodes = nodes)
}

# Function 8: Create and Configure Graph Object
create_graph <- function(all_edges, nodes) {
  g <- graph_from_data_frame(d = all_edges, vertices = nodes, directed = TRUE)
  V(g)$color <- nodes$color[match(V(g)$name, nodes$name)]
  V(g)$size <- 15
  g
}

# Function 9: Assign Edge Colors
assign_edge_colors <- function(g, all_edges,global_receptor_tf_col_max) {
  # For Receptor_TF edges: define a gradient
  edges_rt <- all_edges %>% filter(edge_type == "Receptor_TF")
  max_val <- max(abs(edges_rt$colour), na.rm = TRUE)
  max_value <- max_val + 0.1 * max_val
  gradient_func <- col_numeric(palette = c("blue", "white", "tomato"),
                               domain = c(-max_value, max_value))
  gradient_func <- col_numeric(palette = c("blue", "white", "tomato"),
                             domain = c(-global_receptor_tf_col_max, global_receptor_tf_col_max))

  
  # For Sender_Ligand edges:
  edges_sl <- all_edges %>% filter(edge_type == "Sender_Ligand")
  max_val_sl <- max(abs(edges_sl$colour), na.rm = TRUE)
  max_value_sl <- max_val_sl + 0.1 * max_val_sl
  gradient_func_sl <- col_numeric(palette = c("white", "grey"),
                                  domain = c(0, max_value_sl))
  
  # Assign colors based on edge type
  all_edges <- all_edges %>%
    mutate(
      edge_color = case_when(
        edge_type == "Ligand_Receptor" ~ case_when(
          colour > 0 ~ "tomato",
          colour < 0 ~ "blue",
          TRUE ~ "white"
        ),
        edge_type == "Receptor_TF" ~ gradient_func(colour),
        TRUE ~ gradient_func_sl(colour)
      )
    )
  
  E(g)$color <- all_edges$edge_color
  E(g)$width <- abs(all_edges$weight)
  g
}

# Function 10: Plot the Graph and Save
plot_graph <- function(g, output_file, cluster_name, condition_label) {
  #layout <- layout_with_sugiyama(g)$layout
  #layout <- layout_as_tree(g, root = V(g)[V(g)$layer == "Sender Cell Type"], circular = FALSE)
  # Custom layout: even horizontal spacing per layer
  node_layers <- V(g)$layer
  unique_layers <- c("Sender Cell Type", "Ligand", "Receptor", "TF")

  layout <- matrix(NA, nrow = length(V(g)), ncol = 2)

  for (i in seq_along(unique_layers)) {
    layer <- unique_layers[i]
    nodes_in_layer <- which(node_layers == layer)
    n <- length(nodes_in_layer)
    if (n > 0) {
      x_pad <- 4  # stretch width of layout
      #x_positions <- seq(from = 0, to = x_pad, length.out = n)
      #x_positions <- seq(from = 0, to = n * 0.3, length.out = n)
      x_positions <- seq(from = 0, to = 1, length.out = n + 2)[2:(n + 1)]
      y_position <- length(unique_layers) - i
      layout[nodes_in_layer, 1] <- x_positions
      layout[nodes_in_layer, 2] <- y_position
    }
  }

  png(output_file, width = 10, height = 6.7, units = "in", res = 400)
  plot(g,
       layout = layout,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.9,       # Font size
       vertex.label.font = 2,        # Font weight (2 = bold)
       vertex.label.color = "black",
       main = paste(cluster_name, " (", condition_label, ")", sep = " "),
       edge.arrow.size = 0.5,
       vertex.frame.color = NA,
       asp = 0.5)
  dev.off()
  message("Saved network plot to: ", output_file)
}

# Main Function: Generate Network Plot
generate_network_plot <- function(condition_label, cluster_name, 
                                  decipher_scores, decipher_scores_by_regulon_and_cluster, 
                                  regulon_deltas_by_cluster, feature_statistics,
                                  sender_cts, output_dir,top_interactions = NULL,
                                  global_deltaPagoda_max = global_deltaPagoda_max,
                                  global_receptor_tf_col_max = global_receptor_tf_col_max,
                                  global_sender_ligand_max = global_sender_ligand_max,
                                  global_decipher_score_max = global_decipher_score_max,
                                  n_top_regulons = 10) {
  # 1. Check inputs
  check_cluster_exists(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, regulon_deltas_by_cluster)
  
  # 2. Get top regulons
  top_regulons_out <- get_top_regulons(cluster_name, regulon_deltas_by_cluster, n = n_top_regulons)
  top_regulons_df <- top_regulons_out$top_regulons_df
  top_regulons <- top_regulons_out$top_regulons
  
  # 3. Get top interactions
  top_interactions <- get_top_interactions(cluster_name, decipher_scores, decipher_scores_by_regulon_and_cluster, top_regulons,top_interactions = top_interactions)
  
  # 4. Build Sender → Ligand edges
  sl_out <- build_sender_ligand_edges(top_interactions, feature_statistics, sender_cts,global_sender_ligand_max)
  sender_ligand <- sl_out$sender_ligand
  edges_sender_ligand <- sl_out$edges_sender_ligand
  
  # 5. Build Ligand → Receptor edges
  this_decipher_scores <- decipher_scores[[cluster_name]]
  ligand_receptor_edges <- build_ligand_receptor_edges(this_decipher_scores, sender_ligand, top_interactions,global_decipher_score_max)
  
  # 6. Build Receptor → TF edges
  receptor_tf_edges <- build_receptor_tf_edges(top_interactions,global_receptor_tf_col_max)
  
  # 7. Combine edges and build nodes
  graph_components <- build_graph_components(edges_sender_ligand, ligand_receptor_edges, receptor_tf_edges, top_interactions, sender_ligand, top_regulons_df,global_deltaPagoda_max)
  all_edges <- graph_components$all_edges
  nodes <- graph_components$nodes
  
  # 8. Create graph
  g <- create_graph(all_edges, nodes)
  
  # 9. Assign edge colors and widths
  g <- assign_edge_colors(g, all_edges,global_receptor_tf_col_max)
  
  # 10. Plot and save the network
  output_file <- file.path(output_dir, paste0(condition_label, "_", cluster_name, "_network_map.png"))
  plot_graph(g, output_file, cluster_name, condition_label)
}


# 3. Load Data for Severe and Moderate Conditions
# Severe data
decipher_scores_severe <- readRDS("results/SevCOVID_Azimuthl2/data/decipher_scores_by_cluster.rds")
decipher_scores_by_regulon_and_cluster_severe <- readRDS("results/SevCOVID_Azimuthl2/data/decipher_scores_by_regulon_and_cluster.rds")
regulon_deltas_by_cluster_severe <- readRDS("results/SevCOVID_Azimuthl2/data/regulon_deltas_by_cluster.rds")
feature_statistics_severe <- readRDS("results/SevCOVID_Azimuthl2/data/feature_statistics.rds")

# Mild data
decipher_scores_moderate <- readRDS("results/MilCOVID_Azimuthl2/data/decipher_scores_by_cluster.rds")
decipher_scores_by_regulon_and_cluster_moderate <- readRDS("results/MilCOVID_Azimuthl2/data/decipher_scores_by_regulon_and_cluster.rds")
regulon_deltas_by_cluster_moderate <- readRDS("results/MilCOVID_Azimuthl2/data/regulon_deltas_by_cluster.rds")
feature_statistics_moderate <- readRDS("results/MilCOVID_Azimuthl2/data/feature_statistics.rds")

# 4. Specify Parameters and Generate Networks

# Define target receiver clusters and sender cell types
target_clusters <- c("CD14_Mono", "CD16_Mono")
sender_cts <- c("Eryth", "NK", "cDC2", "CD16_Mono", "CD14_Mono", "CD8_TEM","Platelet","pDC")
output_dir <- "figures_01_08_2025"  # adjust if needed

#calculate global stats
# GLOBAL SCALING VALUES
# TF deltaPagoda
global_deltaPagoda_max <- max(
  sapply(regulon_deltas_by_cluster_severe[target_clusters], function(x) max(x$deltaPagoda, na.rm = TRUE)),
  sapply(regulon_deltas_by_cluster_moderate[target_clusters], function(x) max(x$deltaPagoda, na.rm = TRUE))
)

# Receptor→TF imp.perm * sign(spearman.cor)
global_receptor_tf_col_max <- max(
  sapply(c(decipher_scores_by_regulon_and_cluster_severe[target_clusters], decipher_scores_by_regulon_and_cluster_moderate[target_clusters]), function(cluster_df) {
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
  feature_statistics_severe %>% filter(cluster %in% target_clusters) %>% pull(sum.counts) / feature_statistics_severe %>% filter(cluster %in% target_clusters) %>% pull(n.cell),
  feature_statistics_moderate  %>% filter(cluster %in% target_clusters) %>% pull(sum.counts) / feature_statistics_moderate  %>% filter(cluster %in% target_clusters) %>% pull(n.cell),
  na.rm = TRUE
)
global_sender_ligand_max <- 1


# Ligand→Receptor decipher_score
global_decipher_score_max <- max(
  sapply(decipher_scores_severe[target_clusters], function(x) max(abs(x$decipher_score), na.rm = TRUE)),
  sapply(decipher_scores_moderate[target_clusters], function(x) max(abs(x$decipher_score), na.rm = TRUE))
)

# Generate network plots for Severe condition
for (cl in target_clusters) {

  generate_network_plot("SevCOVID_Azimuthl2", cl,
                      decipher_scores_severe, 
                      decipher_scores_by_regulon_and_cluster_severe,
                      regulon_deltas_by_cluster_severe, 
                      feature_statistics_severe,
                      sender_cts, 
                      output_dir,
                      top_interactions = NULL,
                      global_deltaPagoda_max = global_deltaPagoda_max,
                      global_receptor_tf_col_max = global_receptor_tf_col_max,
                      global_sender_ligand_max = global_sender_ligand_max/1.2,
                      global_decipher_score_max = global_decipher_score_max/1.5,
                      n_top_regulons = 10)

                        }


# Generate network plots for Moderate condition
for (cl in target_clusters) {

  generate_network_plot("MilCOVID_Azimuthl2", cl,
                      decipher_scores_moderate, 
                      decipher_scores_by_regulon_and_cluster_moderate,
                      regulon_deltas_by_cluster_moderate, 
                      feature_statistics_moderate,
                      sender_cts, 
                      output_dir,
                      top_interactions = NULL,
                      global_deltaPagoda_max = global_deltaPagoda_max,
                      global_receptor_tf_col_max = global_receptor_tf_col_max,
                      global_sender_ligand_max = global_sender_ligand_max/1.2,
                      global_decipher_score_max = global_decipher_score_max/1.5,
                      n_top_regulons = 10)

                        }

                        #top_interactions = c("CD38-PECAM1","ENAM-CD63","SLAMF7-SLAMF7", "SIGLEC1-SPN", "CCL3L1-CCR1"))

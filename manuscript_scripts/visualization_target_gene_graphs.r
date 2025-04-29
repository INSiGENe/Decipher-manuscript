library(dplyr)
library(igraph)
library(scales)
library(dplyr)
library(igraph)
library(scales)

#-----------------------------------------------
# Helper: Parse Ligand and Receptor from Interaction
#-----------------------------------------------
parse_interaction_pair <- function(interaction_pair) {
  parts <- strsplit(interaction_pair, "-")[[1]]
  list(ligand = parts[1], receptor = parts[2])
}

#-----------------------------------------------
# Helper: Get Top Regulons for the Interaction
#-----------------------------------------------
get_top_regulons_for_interaction <- function(cluster_name, interaction_pair, decipher_scores_by_regulon_and_cluster, regulon_deltas_by_cluster, n_regulons) {
  regs_df <- decipher_scores_by_regulon_and_cluster[[cluster_name]] %>%
    filter(interaction == interaction_pair) %>%
    arrange(desc(imp.perm)) %>%
    slice_head(n = n_regulons) %>%
    select(regulon, imp.perm)

  regs_info <- regulon_deltas_by_cluster[[cluster_name]] %>%
    filter(name %in% regs_df$regulon) %>%
    rename(regulon = name) %>%
    inner_join(regs_df, by = "regulon") %>%
    select(regulon, imp.perm, deltaPagoda)

  return(regs_info)
}

#-----------------------------------------------
# Helper: Get Top Target Genes for Each Regulon
#-----------------------------------------------
get_top_targets_per_regulon <- function(cluster_name, grns_by_cluster, regs_info, n_targets) {
  grn_df <- grns_by_cluster[[cluster_name]] %>%
    filter(source %in% regs_info$regulon) %>%
    group_by(source) %>%
    slice_max(order_by = coef_abs, n = n_targets, with_ties = FALSE) %>%
    ungroup() %>%
    rename(regulon = source, target_gene = target, weight = coef_abs)

  return(grn_df)
}

#-----------------------------------------------
# Helper: Assemble Edges
#-----------------------------------------------
assemble_edges <- function(ligand, receptor, regs_info, grn_df) {
  edges_lr <- tibble(
    from      = ligand,
    to        = receptor,
    weight    = 1,
    edge_type = "ligand_receptor",
    colour    = NA_real_  # <- NA numeric to match types
  )
  
  edges_rt <- regs_info %>%
    transmute(
      from      = receptor,
      to        = regulon,
      weight    = abs(deltaPagoda),
      edge_type = "receptor_tf",
      colour    = deltaPagoda  # numeric
    )
  
  edges_tg <- grn_df %>%
    transmute(
      from      = regulon,
      to        = target_gene,
      weight    = weight,
      edge_type = "tf_target",
      colour    = weight  # numeric
    )
  
  bind_rows(edges_lr, edges_rt, edges_tg)
}


#-----------------------------------------------
# Helper: Build Node Data Frame
#-----------------------------------------------
build_nodes_df <- function(all_edges, ligand, receptor, regs_info, grn_df) {
  nodes <- unique(c(all_edges$from, all_edges$to))
  
  max_dp   <- max(abs(regs_info$deltaPagoda), na.rm = TRUE)
  max_coef <- max(abs(grn_df$weight), na.rm = TRUE)

  color_tf_palette <- col_numeric(
    palette = c("blue", "white", "tomato"), 
    domain = c(-max(abs(regs_info$deltaPagoda), na.rm = TRUE), 
                max(abs(regs_info$deltaPagoda), na.rm = TRUE))
  )
  color_target_palette <- col_numeric(c("white", "purple"), domain = c(0, max_coef * 1.1))

  nodes_df <- tibble(name = nodes) %>%
    mutate(
      layer = case_when(
        name == ligand             ~ "Ligand",
        name == receptor           ~ "Receptor",
        name %in% regs_info$regulon ~ "TF",
        TRUE                       ~ "Target"
      )
    ) %>%
    left_join(regs_info %>% rename(name = regulon), by = "name") %>%
    left_join(grn_df    %>% rename(name = target_gene), by = "name") %>%
    distinct(name, .keep_all = TRUE) %>%
    mutate(
      size = case_when(
        layer == "Ligand"   ~ 15,
        layer == "Receptor" ~ 15,
        layer == "TF"       ~ rescale(deltaPagoda, c(6, 12), c(-max_dp, max_dp)),
        layer == "Target"   ~ rescale(weight, c(4, 8), c(-max_coef, max_coef)),
        TRUE                ~ 5
      ),
      # 🔥 force safe color assignment
      color = case_when(
        layer == "Ligand"   ~ "skyblue",
        layer == "Receptor" ~ "lightgreen",
        layer == "TF"       ~ ifelse(!is.na(deltaPagoda), color_tf_palette(deltaPagoda), "grey"),
        layer == "Target"   ~ ifelse(!is.na(weight),      color_target_palette(weight), "grey"),
        TRUE                ~ "grey"
      )
    )

  nodes_df <- nodes_df %>%
    mutate(size = pmax(size, 2))   # minimum size 2

}




#-----------------------------------------------
# Helper: Create Graph
#-----------------------------------------------
create_interaction_graph <- function(all_edges, nodes_df, ligand) {
  g <- graph_from_data_frame(all_edges, vertices = nodes_df, directed = TRUE)
  V(g)$color <- nodes_df$color[match(V(g)$name, nodes_df$name)]
  V(g)$size  <- nodes_df$size[match(V(g)$name, nodes_df$name)]
  V(g)$layer <- nodes_df$layer[match(V(g)$name, nodes_df$name)]  # Add layer for label sizing
  
  # Edge color mapping
  color_tf_edges     <- col_numeric(c("blue", "white", "tomato"), domain = c(-20, 20))
  color_target_edges <- col_numeric(c("blue", "white", "tomato"), domain = c(-1, 1))

  # Correct edge color handling:
  E(g)$color <- sapply(seq_len(ecount(g)), function(i) {
    e <- all_edges[i, ]
    if (e$edge_type == "receptor_tf" && !is.na(e$colour)) {
      return(color_tf_edges(as.numeric(e$colour)))
    } else if (e$edge_type == "tf_target" && !is.na(e$colour)) {
      return(color_target_edges(as.numeric(e$colour)))
    } else {
      return("grey")  # << fallback!
    }
  })

  E(g)$width <- rescale(all_edges$weight, to = c(1, 5))
  
  layout <- layout_as_tree(g, root = which(V(g)$name == ligand), circular = TRUE)
  
  list(g = g, layout = layout)
}



#-----------------------------------------------
# Helper: Plot Graph
#-----------------------------------------------
plot_interaction_graph <- function(g, layout, interaction_pair, cluster_name, dataset_name,output_file = NULL) {
  label_size <- ifelse(V(g)$layer == "Target", 0.6, 0.8)  # smaller text for targets
  
  if (!is.null(output_file)) {
    png(output_file, width = 15, height = 15, units = "cm", res = 400)
  }
  plot(
    g,
    layout          = layout,
    vertex.label.cex= label_size,
    edge.arrow.size = 0.4,
    main            = paste0(interaction_pair, " → downstream in ", cluster_name, " (",dataset_name,")")
  )
  if (!is.null(output_file)) dev.off()
}


#=============================================================
# Main Function: generate_interaction_tree (clean version)
#=============================================================
generate_interaction_tree <- function(cluster_name,
                                      interaction_pair,
                                      dataset_name,
                                      decipher_scores_by_regulon_and_cluster,
                                      regulon_deltas_by_cluster,
                                      grns_by_cluster,
                                      n_regulons = 10,
                                      n_targets  = 5,
                                      output_file = NULL) {
  # Step 0 - Check if interaction exists
  if (!interaction_pair %in% decipher_scores_by_regulon_and_cluster[[cluster_name]]$interaction) {
    message(paste0("Skipping: interaction ", interaction_pair, " not found in cluster ", cluster_name))
    return(NULL)
  }
  
  parsed <- parse_interaction_pair(interaction_pair)
  ligand <- parsed$ligand
  receptor <- parsed$receptor
  
  # Step 1 - Get top regulons
  regs_info <- get_top_regulons_for_interaction(cluster_name, interaction_pair, decipher_scores_by_regulon_and_cluster, regulon_deltas_by_cluster, n_regulons)
  
  # Check if regulons found
  if (nrow(regs_info) == 0) {
    message(paste0("Skipping: no regulons downstream of interaction ", interaction_pair, " in cluster ", cluster_name))
    return(NULL)
  }
  
  # Step 2 - Get GRN
  grn_df <- get_top_targets_per_regulon(cluster_name, grns_by_cluster, regs_info, n_targets)
  
  # ⚡ NEW: check if target genes found
  if (nrow(grn_df) == 0) {
    message(paste0("Skipping: no target genes for regulons downstream of interaction ", interaction_pair, " in cluster ", cluster_name))
    return(NULL)
  }
  
  # Step 3 onward
  all_edges <- assemble_edges(ligand, receptor, regs_info, grn_df)
  nodes_df <- build_nodes_df(all_edges, ligand, receptor, regs_info, grn_df)
  graph_list <- create_interaction_graph(all_edges, nodes_df, ligand)
  plot_interaction_graph(graph_list$g, graph_list$layout, interaction_pair, cluster_name, dataset_name,output_file)
}



dataset <- "MilCOVID_Azimuthl2"
dataset_path <- file.path("results",dataset,"data")

# load your GRNs and Decipher results
decipher_by_reg  <- readRDS(file.path(dataset_path,"decipher_scores_by_regulon_and_cluster.rds"))
delta_by_cluster <- readRDS(file.path(dataset_path,"regulon_deltas_by_cluster.rds"))
grns_by_cluster  <- readRDS(file.path(dataset_path,"regulon_grns_by_cluster.rds"))

# generate and save a tree for NK cells, APP–CD74 interaction
selected_cluster <- "cDC2"
selected_interaction <- "C1QB-CD33"
file_name <- file.path("figures",paste0(dataset,"_",selected_cluster,"_",selected_interaction,"_tree.png"))
generate_interaction_tree(
  cluster_name   = selected_cluster,
  interaction_pair = selected_interaction,
  dataset_name = dataset,
  decipher_scores_by_regulon_and_cluster = decipher_by_reg,
  regulon_deltas_by_cluster            = delta_by_cluster,
  grns_by_cluster                     = grns_by_cluster,
  n_regulons  = 10,
  n_targets   = 5,
  output_file = file_name
)

library(devtools)
load_all()
set.seed(123)
library(data.table)

#functions
safe_load <- function(filepath, ...) {
  if (!file.exists(filepath)) {
    message("Skipping missing file: ", filepath)
    return(NULL)
  }

  ext <- tools::file_ext(filepath)

  tryCatch({
    switch(ext,
      "rds" = readRDS(filepath),
      "csv" = read.csv(filepath, ...),
      stop(paste("Unsupported file type:", ext))
    )
  }, error = function(e) {
    message("Error loading file: ", filepath, " — ", e$message)
    return(NULL)
  })
}

add_if_not_null <- function(lst, name, value) {
  if (!is.null(value)) {
    if (is.data.frame(value)) {
      lst[[name]] <- value %>% select(sender, receiver, interaction, prioritization_score, scaled_score)
    } else {
      lst[[name]] <- value  # If it's not a data frame, just store it as-is
    }
  }
  lst
}


prepareForCorrelation <- function(name, df) {
  if (is.null(df)) return(NULL)
  
  if (name == "Decipher") {
    #different for Decipher since scores are already aggregated across sender cell types
    df %>%
      select(interaction, receiver, prioritization_score) %>%
      arrange(prioritization_score)
  } else {
    prepareDataForCorrelationAnalysis(df)
  }
}



# Define your datasets (you can add more as needed)
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
  meta_path <- "manuscript_analysis/data_for_meta_comparisons"
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

install.packages("ggnewscale") 
install.packages("UpSetR")
install.packages("ggbeeswarm")

library(ggnewscale)
library(UpSetR)  
library(ggbeeswarm) 

# --- Assume your previous code has run and populated 'results_preprocessed' ---

# ==== functions ====
# Assumed function to load regulon data (replace with your actual loading mechanism)
# It should return a list where each element corresponds to a cell type,
# and contains a dataframe with 'name' (regulon) and 'deltaPagoda' columns.
load_regulon_data <- function(file_path, cell_types) {
  tryCatch({
      data_list <- readRDS(file_path)
      # Ensure it's filtered for selected cell types if the file contains more
      data_list <- data_list[intersect(names(data_list), cell_types)]
      # Add basic validation
      if(!is.list(data_list)) stop("Loaded data is not a list.")
      if(length(data_list) > 0) {
          first_el <- data_list[[1]]
          if(!is.data.frame(first_el) || !all(c("name", "deltaPagoda") %in% colnames(first_el))) {
              stop("Dataframe structure is incorrect. Needs 'name' and 'deltaPagoda' columns.")
          }
      }
      return(data_list)
  }, error = function(e) {
      warning(paste("Error loading or validating file:", file_path, "-", e$message))
      # Return an empty list or handle appropriately
      return(list())
  })
}

get_deltaPagoda <- function(identity, regulon_list, regulon) {
  identity_char <- as.character(identity)
  if (is.null(regulon_list) || !identity_char %in% names(regulon_list)) {
    # warning(paste("Identity", identity_char, "not found in provided regulon list. Returning NA."))
    return(NA)
  }
  df <- regulon_list[[identity_char]]
  if (is.null(df) || !is.data.frame(df) || !all(c("name", "deltaPagoda") %in% colnames(df))) {
    # warning(paste("Required columns ('name', 'deltaPagoda') missing or data invalid for identity", identity_char, ". Returning NA."))
    return(NA)
  }
  if (regulon %in% df$name) {
    # Handle potential multiple matches (shouldn't happen with unique names)
    return(df$deltaPagoda[df$name == regulon][1])
  } else {
    # warning(paste("Regulon", regulon, "not found for identity", identity_char, ". Returning NA."))
    return(NA)
  }
}

# Helper function to find the absolute max deltaPagoda, handling NULL/NA/empty lists
find_absolute_max <- function(deltas_list_of_lists) {
  max_val <- -Inf # Initialize with negative infinity
  for (cond_list in deltas_list_of_lists) {
      if (!is.null(cond_list) && length(cond_list) > 0) {
          cond_max <- sapply(cond_list, function(df) {
              if (!is.null(df) && is.data.frame(df) && "deltaPagoda" %in% colnames(df) && nrow(df) > 0) {
                  current_max <- max(abs(df$deltaPagoda), na.rm = TRUE)
                  # Handle case where all deltaPagoda are NA or df is empty after NA removal
                  if (is.infinite(current_max)) {
                    return(-Inf) # Return -Inf if no valid values
                  } else {
                    return(current_max)
                  }
              } else {
                  return(-Inf) # Return -Inf for invalid/empty dataframes
              }
          })
          # Filter out -Inf before taking the max for the condition
          valid_cond_max <- cond_max[is.finite(cond_max)]
          if (length(valid_cond_max) > 0) {
             max_val <- max(max_val, max(valid_cond_max, na.rm = TRUE), na.rm = TRUE)
          }
      }
  }
   # If max_val is still -Inf (no valid data found), return a default like 1 or NA
  return(ifelse(is.finite(max_val), max_val, 1))
}

# Function to generate the combined data needed for plotting
generate_heatmap_data_for_celltype <- function(selected_ct, regulon_deltas_list, conditions) {
  heatmap_data <- data.frame()

  # Gather all unique regulons for this cell type across all conditions
  all_regulons <- character(0)
  for (cond_name in names(conditions)) {
    if (!is.null(regulon_deltas_list[[cond_name]]) && !is.null(regulon_deltas_list[[cond_name]][[selected_ct]])) {
      current_regulons <- regulon_deltas_list[[cond_name]][[selected_ct]]$name
      all_regulons <- union(all_regulons, current_regulons)
    } else {
        warning("No data found for cell type '", selected_ct, "' in condition '", cond_name, "'")
    }
  }
  # Remove potential "sample*" regulons if they exist
  all_regulons <- unique(all_regulons[!grepl("sample", all_regulons)])

  if(length(all_regulons) == 0) {
      warning("No valid regulons found for cell type: ", selected_ct)
      return(data.frame()) # Return empty frame
  }

  # Populate the data frame with deltaPagoda values for all regulons and conditions
  for (cond_name in names(conditions)) {
    regulon_deltas <- regulon_deltas_list[[cond_name]]
    for (regulon in all_regulons) {
      mean_delta <- get_deltaPagoda(selected_ct, regulon_deltas, regulon)
      heatmap_data <- rbind(heatmap_data, data.frame(
        TF = regulon,
        Comparison = cond_name,
        DeltaPagoda = mean_delta,
        stringsAsFactors = FALSE # Important for factor handling later
      ))
    }
  }
  return(heatmap_data)
}

# Generate plots (two per cell type: one sorted by moderate, one by severe)
# Helper function: Get top TFs for a given condition
get_top_tfs <- function(data, condition, top_n) {
  data %>%
    filter(Comparison == condition, !is.na(DeltaPagoda)) %>%
    arrange(desc(abs(DeltaPagoda))) %>%
    slice_head(n = top_n) %>%
    arrange(DeltaPagoda) %>%
    pull(TF)
}

# Helper function: Generate ggplot for given TFs
generate_heatmap_plot <- function(data, tfs, condition_names, absolute_max, condition_label, top_n) {
  if (length(tfs) == 0) return(NULL)

  plot_data <- data %>%
    filter(TF %in% tfs) %>%
    mutate(
      TF = factor(TF, levels = tfs),
      Comparison = factor(Comparison, levels = condition_names)
    )

  ggplot(plot_data, aes(x = Comparison, y = TF, fill = DeltaPagoda)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      na.value = "grey80", name = "TF Activity\nDelta",
      limits = c(-absolute_max, absolute_max)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = rel(1.1)),
      axis.text.y = element_text(size = rel(0.9)),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste("Top", top_n, "Regulons (Sorted by", condition_label, ")"))
}

# Main function
generate_sorted_plots <- function(selected_receiver_cells, regulon_deltas_list, conditions, top_n, absolute_max) {
  plots <- list()
  condition_names <- names(conditions)

  for (selected_ct in selected_receiver_cells) {
    cat("Generating plots for cell type:", selected_ct, "\n")
    heatmap_data_full <- generate_heatmap_data_for_celltype(selected_ct, regulon_deltas_list, conditions)

    if (nrow(heatmap_data_full) == 0) {
      warning("Skipping plots for ", selected_ct, " due to lack of data.")
      next
    }

    plots[[selected_ct]] <- list()

    for (cond in c("moderate", "severe")) {
      top_tfs <- get_top_tfs(heatmap_data_full, cond, top_n)

      plot <- generate_heatmap_plot(
        data = heatmap_data_full,
        tfs = top_tfs,
        condition_names = condition_names,
        absolute_max = absolute_max,
        condition_label = cond,
        top_n = top_n
      )

      if (is.null(plot)) {
        warning("No regulons passed filtering for '", cond, "_sorted' plot in ", selected_ct)
      }

      plots[[selected_ct]][[paste0(cond, "_sorted")]] <- plot
    }
  }

  return(plots)
}

create_combined_plots_per_celltype <- function(plots, selected_receiver_cells) {
  combined_plots <- list()
  for (selected_ct in selected_receiver_cells) {
    # Check if both plots exist for this cell type
    plot_moderate <- plots[[selected_ct]][["moderate_sorted"]]
    plotsevere <- plots[[selected_ct]][["severe_sorted"]]

    # Create placeholder plots if one or both are missing
    placeholder_plot <- ggplot() + theme_void() + ggtitle("Data Not Available") + theme(plot.title = element_text(hjust = 0.5))
    if (is.null(plot_moderate)) plot_moderate <- placeholder_plot
    if (is.null(plotsevere)) plotsevere <- placeholder_plot

    # Create title plot
    title <- ggplot() + theme_void() + ggtitle(selected_ct) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

    # Combine the two heatmaps side-by-side
    heatmaps <- wrap_plots(plot_moderate, plotsevere, ncol = 2)

    # Stack title above heatmaps
    combined_plots[[selected_ct]] <- wrap_plots(title, heatmaps, ncol = 1, heights = c(0.1, 1)) # Adjust height ratio for title
  }
  return(combined_plots)
}

save_grouped_plots <- function(combined_plots, clusters_per_group, output_dir_base, output_folder_name) {
  # Construct the full output path
  output_dir <- file.path(output_dir_base, output_folder_name) # Match structure from file 1
  if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
  }

  cluster_names <- names(combined_plots)
  if (length(cluster_names) == 0) {
      warning("No plots available to save.")
      return()
  }

  for (i in seq(1, length(cluster_names), by = clusters_per_group)) {
    # Select clusters for the current group
    selected_clusters_group <- cluster_names[i:min(i + clusters_per_group - 1, length(cluster_names))]

    # Combine the plots for the selected group (e.g., arrange in 2 columns)
    # Adjust ncol based on how many plots per page you want vertically vs horizontally
    plots_to_save_list <- combined_plots[selected_clusters_group]

    # Check if list is empty before plotting
    if(length(plots_to_save_list) > 0){
        to_plot <- wrap_plots(plots_to_save_list, ncol = 1) + # Stack cell types vertically by default
            plot_layout(guides = "collect") & # Collect legends
            theme(
                legend.position = "bottom",
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14, face = "bold"),
                legend.key.size = unit(1.2, "cm"), # Adjust key size
                # Adjust key height/width if needed, size often covers it
                # legend.key.height = unit(1.0, "cm"),
                # legend.key.width = unit(1.0, "cm")
            )

        # Define filename
        group_num <- (i - 1) %/% clusters_per_group + 1
        file_name <- paste0("Combined_Sorted_Heatmaps_Group_", group_num, ".png")
        file_path <- file.path(output_dir, file_name)

        # Calculate dynamic height (adjust base height and multiplier as needed)
        plot_height <- 7 * length(selected_clusters_group) # Base height * number of plots stacked
        plot_width <- 12 # Fixed width (adjust if needed)

        cat("Saving group", group_num, "to:", file_path, "\n")
        # Save the plot
        ggsave(
          filename = file_path,
          plot = to_plot,
          width = plot_width,
          height = plot_height,
          dpi = 300,
          limitsize = FALSE # Important for potentially large plots
        )
    } else {
        cat("Skipping save for group", (i - 1) %/% clusters_per_group + 1, "as no plots were generated for the selected clusters.\n")
    }
  }
}

# 1. Combine both moderate and severe into a long format
get_long_deltas <- function(regulon_deltas, condition) {
  bind_rows(lapply(names(regulon_deltas), function(ct) {
    df <- regulon_deltas[[ct]] %>%
      filter(class == "real")  # Keep only real TFs
    df$Cluster <- ct
    df$Condition <- condition
    return(df)
  }))
}

# Define a pseudo-log transformation that handles negatives
pseudo_log_trans <- function(base = 10) {
  trans_new(
    name = paste0("pseudo_log", base),
    transform = function(x) sign(x) * log1p(abs(x)) / log(base),
    inverse = function(x) sign(x) * (base^abs(x) - 1),
    domain = c(-Inf, Inf)
  )
}

# function that extract the top 'n' interaction identifiers from the input object (which is results_for_correlation[[ds]][[method_name]]).
getSet <- function(method_result, n) {
  # Basic Input Checks
  if (is.null(method_result)) {
    # warning("Input to getSet is NULL. Returning empty set.")
    return(character(0))
  }
  # Check if it's dataframe-like and has rows
  if (!is.data.frame(method_result) && !is(method_result, "DataFrame") && !is.matrix(method_result)) {
     # If it's some other structure, maybe it's already a vector of interactions?
     if(is.character(method_result)) {
         return(head(method_result, n))
     }
     warning("Input to getSet is not a recognized dataframe-like structure or character vector. Returning empty set.")
     return(character(0))
  }
   if (nrow(method_result) == 0) {
    # warning("Input to getSet has 0 rows. Returning empty set.")
    return(character(0))
  }

  # Check for 'interaction' column - Adapt column name if necessary!
  interaction_col <- "interaction" # <--- CHANGE THIS if your column name is different
  if (!interaction_col %in% colnames(method_result)) {
     warning(paste("Column '", interaction_col, "' not found in data for a method. Returning empty set.", sep=""))
     return(character(0))
  }

  # Extract top N interactions
  # Ensure the column is character type
  score_col <- "prioritization_score"
  sorted_indices <- order(method_result[[score_col]], decreasing = TRUE, na.last = TRUE)
  sorted_method_result <- method_result[sorted_indices, , drop = FALSE] # Use drop=FALSE for safety
  interactions <- as.character(sorted_method_result[[interaction_col]])
  return(head(interactions, n))
}

# --- Helper function to calculate sizes of ALL specific intersections ---
# Input: list_input = named list of character vectors (sets) for 5 methods
# Output: dataframe with Intersection_Name and Count for all 31 non-empty intersections
calculate_all_intersections <- function(list_input) {
  method_names <- names(list_input)
  expected_methods <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")
  n_methods <- length(expected_methods)

  # Ensure all 5 methods are potentially present, use empty set if NULL or missing
  sets <- list()
  for(m in expected_methods) {
    sets[[m]] <- if (!is.null(list_input[[m]]) && length(list_input[[m]]) > 0) list_input[[m]] else character(0)
  }

  # List to store results: Intersection Name -> Count
  all_intersections_info <- list()

  # Loop through all possible intersection sizes (k = 1 to 5)
  for (k in 1:n_methods) {
    # Generate all combinations (subsets) of method names of size k
    combinations_k <- combn(expected_methods, k, simplify = FALSE)

    # Process each combination (subset S)
    for (subset_S in combinations_k) {
      # Sort subset for consistent naming
      subset_S_sorted <- sort(subset_S)
      # Create the intersection name
      intersection_name <- paste(subset_S_sorted, collapse = " & ")

      # Identify methods NOT in the current subset (subset NotS)
      subset_NotS <- setdiff(expected_methods, subset_S)

      # --- Calculate the count for elements ONLY in subset_S ---

      # 1. Find the intersection of all sets IN subset_S
      # Need to handle case where subset_S has only 1 element
      if (length(subset_S) == 1) {
        intersect_S <- sets[[subset_S[[1]]]]
      } else {
        intersect_S <- Reduce(intersect, sets[subset_S])
      }

      # If the intersection is already empty, the specific count is 0
      if (length(intersect_S) == 0) {
        count <- 0
      } else {
        # 2. Find the union of all sets NOT in subset_S (if any)
        union_NotS <- character(0) # Initialize as empty set
        if (length(subset_NotS) > 0) {
          union_NotS <- Reduce(union, sets[subset_NotS])
        }

        # 3. Find elements in intersect_S that are NOT in union_NotS
        specific_intersect_S_elements <- setdiff(intersect_S, union_NotS)
        count <- length(specific_intersect_S_elements)
      }

      # Store result
      all_intersections_info[[intersection_name]] <- count
    }
  }

  # Convert list to dataframe
  results_df <- data.frame(
    Intersection_Name = names(all_intersections_info),
    Count = unlist(all_intersections_info)
  )
  rownames(results_df) <- NULL # Clean up row names
  return(results_df)
}

# Function: Create a custom color scale for flagged points.
# For each method, assign:
#   - A darker color when flagged (flagged == TRUE)
#   - The base color when not flagged (flagged == FALSE)
create_flag_color_scale <- function(methods) {
  # Get a base color for each method using a hue palette.
  #base_cols <- scales::hue_pal()(length(methods))
  base_cols <- method_colors
  #names(base_cols) <- methods
  names(base_cols) <- names(method_colors)
  
  # Helper to darken a given color.
  darker <- function(col) {
    adjustcolor(col, red.f = 0.6, green.f = 0.6, blue.f = 0.6)
  }
  
  # For each method, return a vector with:
  #   - A darker color (flagged)
  #   - The base color (not flagged)
  color_values <- unlist(lapply(methods, function(m) {
    c(darker(base_cols[m]), base_cols[m])
  }))
  
  # Name the colors to match the values produced by interaction(method, flagged).
  # That is, names like "MethodName.TRUE" and "MethodName.FALSE".
  names(color_values) <- unlist(lapply(methods, function(m) {
    c(paste0(m, ".TRUE"), paste0(m, ".FALSE"))
  }))
  
  return(color_values)
}

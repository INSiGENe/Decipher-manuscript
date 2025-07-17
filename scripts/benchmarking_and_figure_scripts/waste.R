

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
    ggtitle(paste("Top", top_n, "TFs (Sorted by ", condition_label, ")"))
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

    for (cond in c("mild", "severe")) {
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
    plot_mild <- plots[[selected_ct]][["mild_sorted"]]
    plotsevere <- plots[[selected_ct]][["severe_sorted"]]

    # Create placeholder plots if one or both are missing
    placeholder_plot <- ggplot() + theme_void() + ggtitle("Data Not Available") + theme(plot.title = element_text(hjust = 0.5))
    if (is.null(plot_mild)) plot_mild <- placeholder_plot
    if (is.null(plotsevere)) plotsevere <- placeholder_plot

    # Create title plot
    title <- ggplot() + theme_void() + ggtitle(selected_ct) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

    # Combine the two heatmaps side-by-side
    heatmaps <- wrap_plots(plot_mild, plotsevere, ncol = 2)

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
  # Get unique interactions only
  unique_interactions <- unique(as.character(sorted_method_result[[interaction_col]]))
  return(head(unique_interactions, n))

  #return(head(interactions, n))
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

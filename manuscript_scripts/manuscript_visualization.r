library(ggplot2)
library(dplyr)     # Or library(data.table) if you prefer
library(purrr)     # For map functions
library(stringr)   # For string manipulation if needed
library(patchwork) # For combining plots (optional)
library(ggridges)  # For density ridgeline plots
library(UpSetR) # Make sure UpSetR is installed install.packages("UpSetR")
library(tidyr)
# --- Assume your previous code has run and populated 'results_preprocessed' ---

# Combine results into a single dataframe for plotting
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
  # Add more if needed
)

# Optional: Remove datasets/methods with zero rows if any slipped through
combined_scores_df <- combined_scores_df %>% filter(!is.na(prioritization_score))


desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# Convert the 'method' column to a factor with the specified levels.
# We use rev() because ggplot plots factors bottom-up on the y-axis.
# To have "Decipher" at the top visually, it needs to be the last level.
combined_scores_df <- combined_scores_df %>%
  mutate(method = factor(method, levels = rev(desired_method_order)))







############################
###### Box-plot of number of interactions reported on by each method ##
############################


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
interaction_counts$method <- factor(interaction_counts$method, levels = desired_method_order)

# Plot
p <- ggplot(interaction_counts, aes(x = method, y = interaction_count, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = method_colors) +
  scale_y_log10() +
  labs(
    x = "Method",
    y = "Number of Interactions (log10)"
  ) +
  theme_minimal(base_size = 14) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 1,size = 14), # Rotate labels for readability
    axis.text.y = element_text(size=14)
  )


ggsave("figures/boxplot_number_of_interactions_by_method.png", plot = p, width = 8, height = 4, dpi = 300)






############################
###### violin plot ##
############################
# Load necessary library
library(ggplot2)
library(scales) # For pseudo_log_trans

# --- Check if variables exist ---
if (!exists("combined_scores_df")) stop("DataFrame 'combined_scores_df' not found.")
if (!exists("method_colors")) stop("Color vector 'method_colors' not found.")
if (!exists("desired_method_order")) stop("Order vector 'desired_method_order' not found.")

# --- Generate the Violin Plot ---
plot_violin_scores_by_method <- ggplot(combined_scores_df,
                                        # Swap x and y aesthetics
                                        aes(x = method,
                                            y = prioritization_score,
                                            fill = method)) + # Grouping defaults to x aesthetic

  # Use geom_violin instead of geom_density_ridges
  geom_violin(trim = FALSE,    # Keep tails
              alpha = 0.7,     # Adjust transparency
              scale = "width", # Makes violins comparable across methods
              # Optional: Show median and quartiles
              draw_quantiles = c(0.25, 0.5, 0.75)
              ) +

  # Use consistent colors; ensure names match factor levels
  scale_fill_manual(values = method_colors, name = "Method", breaks = desired_method_order) +

  # Apply pseudo-log scale to the Y axis
  scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
      breaks = c(-10, -1, 0, 1, 10, 100), # Adjust breaks if needed for score range
      name = "Interaction Prioritization Score (Pseudo-Log Scale)" # Set y-axis label here
  ) +

  # Ensure method order on the X axis
  scale_x_discrete(limits = desired_method_order, name = "Method") + # Set x-axis label here

  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", # Keep legend off
    axis.text.x = element_text(angle = 0, vjust = 1,size = 14), # Rotate labels for readability
    axis.text.y = element_text(size=14)
    # Add other theme adjustments if needed
  ) +
  labs(
    # x and y labels are set within the scale functions now
    fill = "Method"
  )


# --- Save the Plot ---
# Adjust dimensions as needed for violin plot layout
output_filename_violin <- "figures/violin_scores_by_method.png"
dir.create("figures", showWarnings = FALSE)

save_result_violin <- tryCatch({
    ggsave(output_filename_violin, plot = plot_violin_scores_by_method, width = 8, height = 4, dpi = 300)
    TRUE
}, error = function(e) {
    warning("Failed to save the violin plot: ", e$message, call. = FALSE)
    FALSE
})

if (save_result_violin) {
    print(paste("Violin plot saved to", output_filename_violin))
} else {
    print("Violin plot was generated but could not be saved automatically.")
}













#############################
######### overlap ###########
#############################
# function that extract the top 'n' interaction identifiers
# from the input object (which is results_for_correlation[[ds]][[method_name]]).

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

# --- Main Calculation Loop ---
n_top <- 100
method_names_expected <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")
all_intersection_data <- list()

# Ensure 'results_for_correlation' list is populated
if (!exists("results_for_correlation") || length(results_for_correlation) == 0) {
    stop("The 'results_for_correlation' list was not found or is empty.")
}

print(paste("Processing", length(results_for_correlation), "datasets for full overlap analysis..."))
dataset_names_to_process <- names(results_for_correlation)

for (ds in dataset_names_to_process) {
    print(paste("Calculating overlaps for:", ds))
    dataset_results <- results_for_correlation[[ds]]

    if (is.null(dataset_results) || !is.list(dataset_results) || length(dataset_results) == 0) {
        warning(paste("Skipping dataset", ds, "- invalid entry."))
        next
    }

    list_input_current_ds <- list()
    available_methods_in_ds <- names(dataset_results)
    for (method_name in method_names_expected) {
        if (method_name %in% available_methods_in_ds) {
             list_input_current_ds[[method_name]] <- tryCatch({
                  getSet(dataset_results[[method_name]], n_top)
             }, error = function(e) {
                  warning(paste("Error in getSet for", method_name, ds, ":", e$message)); return(character(0))
             })
        } else {
             list_input_current_ds[[method_name]] <- character(0)
        }
    }

    # *** Use the new function here ***
    dataset_summary <- tryCatch({
         calculate_all_intersections(list_input_current_ds) # Call the new function
    }, error = function(e) {
         warning(paste("Error calculating all intersections for", ds, ":", e$message)); return(NULL)
    })

    if (!is.null(dataset_summary)) {
        dataset_summary$Dataset <- ds
        all_intersection_data[[ds]] <- dataset_summary
    } else {
         warning(paste("No intersection summary generated for dataset:", ds))
    }
}

# --- Combine Results and Plot ---
if (length(all_intersection_data) > 0) {
    combined_intersection_df <- bind_rows(all_intersection_data)

    # --- Define Logical Order for Plotting (by Degree, then Name) ---
    # Calculate degree (number of methods in the intersection)
    combined_intersection_df <- combined_intersection_df %>%
      mutate(Degree = sapply(strsplit(Intersection_Name, " & "), length))

    # Order the dataframe to get factor levels easily
    ordered_df <- combined_intersection_df %>%
      arrange(desc(Degree), Intersection_Name) # Order by Degree (high to low), then alphabetically

    # Get the unique names in the desired order
    intersection_order_all <- unique(ordered_df$Intersection_Name)

    # Apply the factor levels
    combined_intersection_df <- combined_intersection_df %>%
        mutate(Intersection_Name = factor(Intersection_Name, levels = intersection_order_all)) %>%
        filter(!is.na(Intersection_Name)) # Clean up just in case

    if(nrow(combined_intersection_df) == 0) {
        stop("No valid intersection data remained after ordering.")
    }

    # Save the combined data used for the plot
    saveRDS(combined_intersection_df, file = "figures/overlap_all_intersections_data.rds")
    print("Full overlap data saved to figures/overlap_all_intersections_data.rds")

} else {
    warning("No intersection data was successfully calculated for any dataset. Cannot create plot.")
}

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
png("figures/overlap_upset_median_plot_ordered.png", width = 10, height = 6, units = "in", res = 300)

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
      list(query = intersects, params = list("Decipher", "Connectome","LIANA+"), color = method_colors["Decipher"], active = T),
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
## Correlation ###########
##########################

library(ggplot2)
library(reshape2)   # For reshaping matrix to long format
library(gridExtra)  # For arranging multiple heatmaps in a grid
library(dplyr)

# Initialize list to store Spearman matrices
spearman_matrices <- list()

# Loop through each dataset
for (dataset in names(results_for_correlation)) {
  
  # Compute correlation & search space
  interaction_results_correlation_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(results_for_correlation[[dataset]])
  
  # Extract Spearman correlation matrix
  spearman_matrices[[dataset]] <- interaction_results_correlation_search_space$spearman
}

# Initialize an empty list for spearman matrix data
spearman_long_list <- list()

# Loop through datasets to extract method pairs
for (dataset in names(spearman_matrices)) {
  
  # Convert Spearman matrix to long format
  spearman_long <- melt(spearman_matrices[[dataset]])
  
  # Add dataset name for tracking
  spearman_long$Dataset <- dataset
  
  # Store in the list
  spearman_long_list[[dataset]] <- spearman_long
}

# Combine all datasets into one long data frame
combined_spearman_df <- do.call(rbind, spearman_long_list)

# Convert factors to character before processing
combined_spearman_df <- combined_spearman_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Method1 = sub("-.*", "", paste(Var1, Var2, sep = "-")),  # Extract first method
         Method2 = sub(".*-", "", paste(Var1, Var2, sep = "-")),  # Extract second method
         Method_Pair = paste(Var1, Var2, sep = "-"))

#search space -----
# Initialize list to store k-matrices
k_matrices <- list()

# Loop through each dataset
for (dataset in names(results_for_correlation)) {
  
  # Compute correlation & search space
  interaction_results_correlation_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(results_for_correlation[[dataset]])
  
  # Extract k-matrix
  k_matrices[[dataset]] <- interaction_results_correlation_search_space$k_matrix
}

# Initialize an empty list for k-matrix data
k_matrix_long_list <- list()

# Loop through datasets to extract method pairs
for (dataset in names(k_matrices)) {
  
  # Convert k-matrix to long format
  k_long <- melt(k_matrices[[dataset]])
  
  # Add dataset name for tracking
  k_long$Dataset <- dataset
  
  # Store in the list
  k_matrix_long_list[[dataset]] <- k_long
}

# Combine all datasets into one long data frame
combined_k_df <- do.call(rbind, k_matrix_long_list)

# Convert factors to character before sorting
combined_k_df <- combined_k_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Method_Pair = paste(Var1, Var2, sep = "-"))  # Keep both orders

# Order alphabetically first by Var1, then by Var2 (but keep both orders)
combined_k_df <- combined_k_df %>% arrange(Var1, Var2)

# Split method pairs into two separate columns for hierarchical labeling
combined_k_df <- combined_k_df %>%
  mutate(Method1 = sub("-.*", "", Method_Pair),  # Extract first method
         Method2 = sub(".*-", "", Method_Pair))  # Extract second method

# Initialize lists for both matrices
k_matrix_long_list <- list()
spearman_long_list <- list()

# Loop through datasets to extract method pairs for both k-matrix and Spearman
for (dataset in names(k_matrices)) {
  
  # Process k-matrix
  k_long <- melt(k_matrices[[dataset]], varnames = c("Var1", "Var2"), value.name = "k_value") %>%
    mutate(Dataset = dataset)
  
  # Process Spearman matrix
  spearman_long <- melt(spearman_matrices[[dataset]], varnames = c("Var1", "Var2"), value.name = "Spearman") %>%
    mutate(Dataset = dataset)
  
  # Store in lists
  k_matrix_long_list[[dataset]] <- k_long
  spearman_long_list[[dataset]] <- spearman_long
}

# Combine all datasets into one long data frame for k-values and Spearman correlations
combined_k_df <- do.call(rbind, k_matrix_long_list)
combined_spearman_df <- do.call(rbind, spearman_long_list)

# Convert factors to characters before processing
combined_k_df <- combined_k_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Method1 = sub("-.*", "", paste(Var1, Var2, sep = "-")),
         Method2 = sub(".*-", "", paste(Var1, Var2, sep = "-")),
         Method_Pair = paste(Var1, Var2, sep = "-"))

combined_spearman_df <- combined_spearman_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         Method1 = sub("-.*", "", paste(Var1, Var2, sep = "-")),
         Method2 = sub(".*-", "", paste(Var1, Var2, sep = "-")),
         Method_Pair = paste(Var1, Var2, sep = "-"))

# Merge the two data frames by Dataset and Method_Pair
combined_df <- full_join(combined_k_df, combined_spearman_df, by = c("Dataset", "Method_Pair", "Var1", "Var2", "Method1", "Method2"))

# --- Load necessary libraries ---
library(ggplot2)
library(dplyr)

# --- Assume 'combined_df' is ready ---
# It should have columns like: Var1, Var2, k_value, Spearman

# --- Define colors and desired order ---
color_k_value <- "#E69F00" # Gold/Orange
color_spearman <- "#56B4E9" # Light Blue
# Make sure this exists in your environment
desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# --- 1. Filter out self-comparisons AND Prepare Method Columns ---
# Ensure Method1/Method2 exist, using Var1/Var2 if necessary
if (!all(c("Method1", "Method2") %in% names(combined_df))) {
   combined_df <- combined_df %>% mutate(Method1 = Var1, Method2 = Var2)
   warning("Recreated Method1/Method2 columns from Var1/Var2.", call. = FALSE)
}

plot_data <- combined_df %>%
  filter(Method1 != Method2) %>% # Keep only pairs where methods are different
  # Ensure methods are within the desired order list
  filter(Method1 %in% desired_method_order & Method2 %in% desired_method_order)

if(nrow(plot_data) == 0) {
  stop("Filtering removed all data. Check input 'combined_df' and filtering conditions.")
}


# --- 2. & 3. Order Facets (Method1) and X-axis (Method2) ---
# Convert to factors WITH the desired levels
plot_data$Method1 <- factor(plot_data$Method1, levels = desired_method_order)
plot_data$Method2 <- factor(plot_data$Method2, levels = desired_method_order)

# Check if conversion worked (optional)
# print(levels(plot_data$Method1))
# print(levels(plot_data$Method2))
# print(head(plot_data))


# --- Create the final dual-axis boxplot ---
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

# Print the final plot
print(p_final)

# Save the plot (adjust dimensions if needed)
ggsave("figures/combined_k_spearman_boxplot_final_ordered.png", plot = p_final, width = 16, height = 7.5, dpi = 300)









##########################
### AUC plots ####
########################
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


#library(gplots)

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


#ROC curves

library(ggplot2)
library(pROC)
library(patchwork)
library(data.table)

#results_for_comparison <- results_for_comparison[setdiff(names(results_for_comparison), "lupus")]


predictions_and_responses_all <- list()
auc_scores_by_datset <- list()
for (ds in names(datasets)) {
    dataset_path <- datasets[[ds]]
    pre_processing_filepath <- file.path(dataset_path, "pre_processing")
    meta_path <- "manuscript_analysis/data_for_meta_comparisons"
    output_figures_filepath <- file.path(dataset_path, "figures")
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


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
# If not installed yet:
# install.packages("ggbeeswarm")
library(ggbeeswarm)
library(dplyr)

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

# Filter flagged points (n_true < 10)
flagged_points <- filter(results_df, flag == "*")

# Plot
results_df$flagged <- results_df$flag == "*"

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


#color line by variance
library(dplyr)

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
ggsave("figures/variance_w_boxplot_beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)

#########################
#V2

# --- Load necessary libraries ---
library(ggplot2)
library(ggbeeswarm) # Ensure this is loaded for geom_beeswarm
library(dplyr)      # For data manipulation (filtering/mutating)
#install.packages("ggnewscale") #add this to docker
library(ggnewscale) # For multiple color scales


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
  geom_line(aes(group = dataset, color = line_color), size = 1, alpha = 0.6) +
  scale_color_manual(
    # name = "Line Group", # Optional legend name
    values = line_color_map, # Use the map created above
    guide = "none" # Hide legend for lines
  ) +

  # *** Introduce a new scale for color ***
  new_scale_color() +

  # 2. Boxplot - Neutral colors (Plot after lines, before points?)
  geom_boxplot(outlier.shape = NA, width = 0.25, alpha = 0.4, color = "black", fill = "lightgray") +

  # 3. Beeswarm Points - Mapped to 'method', using second color scale
  geom_beeswarm(
    aes(color = method),   # Color points by method
    size = 4.0,            # Adjust point size if needed
    cex = 1,               # Use 'size', not 'cex' for ggplot point size control
    priority = "density",
    groupOnX = TRUE        # Ensure beeswarm groups correctly on discrete x-axis
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
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) # Rotate labels

# Print the plot
print(p_updated)

# --- Save the Plot ---
output_filename_updated <- "figures/variance_w_boxplot_beeswarm_auc_plot_updated.png"
dir.create("figures", showWarnings = FALSE)

save_result_updated <- tryCatch({
    ggsave(output_filename_updated, plot = p_updated, width = 5, height = 6, dpi = 300) # Adjusted width slightly
    TRUE
}, error = function(e) {
    warning("Failed to save the updated plot: ", e$message, call. = FALSE)
    FALSE
})

if (save_result_updated) {
    print(paste("Updated plot saved to", output_filename_updated))
} else {
    print("Updated plot was generated but could not be saved automatically.")
}




#basic

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

ggsave("figures/beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)

















##########################
### Decipher heatmap ####
########################
# --- Load Necessary Libraries ---
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This script requires 'ggplot2'. Please install it.", call. = FALSE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("This script requires 'patchwork'. Please install it.", call. = FALSE)
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("This script requires 'tidyr'. Please install it.", call. = FALSE)
}
library(ggplot2)
library(patchwork)
library(tidyr)
# library(dplyr) # Not strictly required if using base R alternatives below

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
# library(dplyr)
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
output_filename <- "figures/decipher_top20_heatmap.png"
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













###################################
##### COVID Severe vs Mild#########
###################################
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel) # Although not used in heatmap, kept for consistency if extending
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration & Placeholders
# -----------------------------------------------------------------------------

# Base path where condition-specific result folders are located
# Adjust this path to your actual directory structure
base_comparison_path <- "Manuscript_jan_2025"
results_path <- file.path(base_comparison_path, "results")
# Subfolder name within the base path for saving results
# This will contain the final grouped PNG images
output_folder_name <- "figures" 

# Define the specific cell types (receiver cells) you want to analyze
# Example using B cell types from file 1:
selected_receiver_cells <- c( "Eryth", "CD16_Mono", "HSPC",    "CD4_TCM",  "Plasmablast",    "B_intermediate", "B_naive","CD8_Naive","NK","CD8_TEM","pDC","cDC2","Platelet","CD14_Mono","CD4_CTL" )
# Example using CD8 T cell types from file 3:
# selected_receiver_cells <- c("Naive_CD8", "CM_CD8", "EM_CD8", "Activ_CTL", "CTL_CD8", "GZMK_CD8", "ISG_CTL_CD8")
# Example using NK cell types from file 1:
# selected_receiver_cells <- c("ISG_NK","NKT","CD56_dim_NK","CD56_brt_CCL4_NK","CD56_brt_NK","Early_NK","NK","Activ_NK")

# Number of top regulons (TFs) to display in each heatmap
top_n_regulons <- 20

# Number of cell types to combine per output PNG file
clusters_per_group_in_output <- 1 # Adjust as needed (e.g., 1, 2, 3)

# Define condition names and the subfolders where their data resides
# The names ('moderate', 'severe') will be used throughout the script
conditions <- c(
  moderate =  "MilCOVID_Azimuthl2", # Folder name for moderate data
  severe   =  "SevCOVID_Azimuthl2"    # Folder name for severe data
)

# -----------------------------------------------------------------------------
# Helper Functions (Adapted from provided files)
# -----------------------------------------------------------------------------

# Assumed function to load regulon data (replace with your actual loading mechanism)
# It should return a list where each element corresponds to a cell type,
# and contains a dataframe with 'name' (regulon) and 'deltaPagoda' columns.
load_regulon_data <- function(file_path, cell_types) {
  # --- This is a PLACEHOLDER ---
  # Replace this with your actual code to load the .rds file
  # and subset/process it as needed.
  # Example:
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
  # --- End PLACEHOLDER ---
}


# Function to get deltaPagoda for a specified regulon and identity (cell type)
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


# -----------------------------------------------------------------------------
# Data Loading
# -----------------------------------------------------------------------------
# Dynamically load data based on conditions defined above
regulon_deltas_list <- lapply(names(conditions), function(cond_name) {
  folder_name <- conditions[[cond_name]]
  file_path <- file.path(results_path, folder_name, "data/regulon_deltas_by_cluster.rds")
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


# -----------------------------------------------------------------------------
# Calculate Global Color Scale Limit
# -----------------------------------------------------------------------------

absolute_max <- find_absolute_max(regulon_deltas_list)
cat("Global absolute max deltaPagoda for scaling:", absolute_max, "\n")
# Ensure absolute_max is not zero or negative, set a minimum limit if needed
if(is.na(absolute_max) || absolute_max <= 0) {
    warning("Could not determine a valid absolute max deltaPagoda. Setting scale limit to 1.")
    absolute_max <- 1
}

# -----------------------------------------------------------------------------
# Heatmap Generation
# -----------------------------------------------------------------------------

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
generate_sorted_plots <- function(selected_receiver_cells, regulon_deltas_list, conditions, top_n, absolute_max) {
  plots <- list()
  condition_names <- names(conditions)

  for (selected_ct in selected_receiver_cells) {
    cat("Generating plots for cell type:", selected_ct, "\n")
    # Get the combined data for this cell type
    heatmap_data_full <- generate_heatmap_data_for_celltype(selected_ct, regulon_deltas_list, conditions)

    if (nrow(heatmap_data_full) == 0) {
        warning("Skipping plots for ", selected_ct, " due to lack of data.")
        next # Skip to the next cell type
    }

    plots[[selected_ct]] <- list()

    # Create plot sorted by 'moderate'

      # Determine top N TFs based on 'moderate' condition
      heatmap_data_filtered_moderate <- heatmap_data_full %>%
        filter(Comparison == "Moderate" & !is.na(DeltaPagoda)) %>%
        arrange(desc(abs(DeltaPagoda))) %>% # Sort by absolute value first to get strongest effects
        #arrange(desc(DeltaPagoda)) %>% # Alternative: Sort by signed value
        slice_head(n = top_n) %>% # Take top N based on chosen sort
        arrange(DeltaPagoda) # Arrange for y-axis order (ascending is common)

      # Get the ordered list of TFs
      ordered_tfs_moderate <- heatmap_data_filtered_moderate$TF

      if(length(ordered_tfs_moderate) > 0) {
          # Filter the full data to include only these TFs
          plot_data_moderate_sorted <- heatmap_data_full %>%
            filter(TF %in% ordered_tfs_moderate) %>%
            mutate(
              TF = factor(TF, levels = ordered_tfs_moderate), # Set factor levels for y-axis order
              Comparison = factor(Comparison, levels = condition_names) # Ensure consistent x-axis order
            )

          # Generate the ggplot object
          plots[[selected_ct]][["moderate_sorted"]] <- ggplot(plot_data_moderate_sorted, aes(x = Comparison, y = TF, fill = DeltaPagoda)) +
            geom_tile(color = "white", linewidth = 0.5) + # Added line thickness
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80", name = "TF Activity\nDelta", limits = c(-absolute_max, absolute_max)) +
            theme_minimal(base_size = 14) + # Increased base size
            theme(
              axis.text.x = element_text(size = rel(1.1)), # Relative sizing
              axis.text.y = element_text(size = rel(0.9)),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              legend.position = "bottom",
              plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5) # Centered title
            ) +
            ggtitle(paste("Top", top_n, "Regulons (Sorted by moderate)"))
      } else {
          warning("No regulons passed filtering for 'moderate_sorted' plot in ", selected_ct)
          plots[[selected_ct]][["moderate_sorted"]] <- NULL # Indicate missing plot
      }
  }

  # Create plot sorted by 'severe'
  # Determine top N TFs based on 'severe' condition
  heatmap_data_filtered_severe <- heatmap_data_full %>%
    filter(Comparison == "Severe" & !is.na(DeltaPagoda)) %>%
    arrange(desc(abs(DeltaPagoda))) %>% # Sort by absolute value
    #arrange(desc(DeltaPagoda)) %>% # Alternative: sort by signed value
    slice_head(n = top_n) %>%
    arrange(DeltaPagoda) # Arrange for y-axis

  ordered_tfs_severe <- heatmap_data_filtered_severe$TF

  if(length(ordered_tfs_severe) > 0) {
      # Filter the full data to include only these TFs
      plot_data_severe_sorted <- heatmap_data_full %>%
        filter(TF %in% ordered_tfs_severe) %>%
        mutate(
          TF = factor(TF, levels = ordered_tfs_severe),
          Comparison = factor(Comparison, levels = condition_names)
        )

      plots[[selected_ct]][["severe_sorted"]] <- ggplot(plot_data_severe_sorted, aes(x = Comparison, y = TF, fill = DeltaPagoda)) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80", name = "TF Activity\nDelta", limits = c(-absolute_max, absolute_max)) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(size = rel(1.1)),
          axis.text.y = element_text(size = rel(0.9)),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5)
        ) +
        ggtitle(paste("Top", top_n, "Regulons (Sorted by severe)"))
  } else {
      warning("No regulons passed filtering for 'severe_sorted' plot in ", selected_ct)
      plots[[selected_ct]][["severe_sorted"]] <- NULL
  }
  return(plots)
}


#re-factored code
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



# Execute plot generation
generated_plots <- generate_sorted_plots(selected_receiver_cells, regulon_deltas_list, conditions, top_n_regulons, absolute_max)

# -----------------------------------------------------------------------------
# Combine and Save Plots
# -----------------------------------------------------------------------------

# Create combined plots with titles (one combined plot per cell type)
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

celltype_combined_plots <- create_combined_plots_per_celltype(generated_plots, selected_receiver_cells)

# Function to save combined plots in groups
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

# Run the function to save the plots in groups
save_grouped_plots(
    combined_plots = celltype_combined_plots,
    clusters_per_group = clusters_per_group_in_output,
    output_dir_base = base_comparison_path, # Use the base path
    output_folder_name = output_folder_name
)

cat("Script finished. Plots saved in:", file.path(base_comparison_path, , output_folder_name), "\n")

# now DecipherPlots
plotDecipherPrioritizedMap("results/SevCOVID_Azimuthl2",top_n=10,dataset_name="SevCOVID_Azimuthl2",direction = "pos",selected_receiver_cells = c("CD16_Mono"),split_by_direction = TRUE)
plotDecipherPrioritizedMap("results/MilCOVID_Azimuthl2",top_n=10,dataset_name="MilCOVID_Azimuthl2",direction = "pos",selected_receiver_cells = c("CD16_Mono"),split_by_direction=TRUE)















##############################
#### heatmap of signalling by cell-type
##############################
# --- Load Necessary Libraries ---
# Ensure these libraries are installed: install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- Ensure 'results_preprocessed' list exists ---
if (!exists("results_preprocessed") || !is.list(results_preprocessed) || length(results_preprocessed) == 0) {
    stop("The 'results_preprocessed' list is empty, not a list, or does not exist. Please ensure the initial data loading script has run successfully.")
}

# --- Define Desired Method Order ---
desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# --- 1. Aggregate Data and Calculate Percentages ---
print("Aggregating scores and calculating percentages...")
all_percentages_list <- list()

# Get the list of datasets that were actually processed (had folders and potentially data)
processed_datasets <- names(results_preprocessed)

for (ds in processed_datasets) {
    dataset_results <- results_preprocessed[[ds]]
    
    # Check if dataset results exist
    if (is.null(dataset_results) || length(dataset_results) == 0) {
        warning(paste("No results found for dataset:", ds, "- skipping."))
        next
    }
    
    available_methods <- names(dataset_results)
    
    for (method in desired_method_order) {
        
        # Check if the method exists for this dataset
        if (!method %in% available_methods) {
             # It's okay if a method is missing for a dataset, we'll handle NAs later
             # message(paste("Method '", method, "' not found for dataset:", ds))
             next
        }

        method_data <- dataset_results[[method]]
        
        # Check if data exists and has the required columns
        if (!is.null(method_data) && is.data.frame(method_data) && nrow(method_data) > 0 && all(c("receiver", "prioritization_score") %in% names(method_data))) {
            
            # Filter for positive scores only
            positive_scores_df <- method_data %>%
                filter(prioritization_score > 0)
                
            # Check if there are any positive scores
            if (nrow(positive_scores_df) > 0) {
                
                # Calculate total positive score for this method in this dataset
                total_dataset_method_score <- sum(positive_scores_df$prioritization_score, na.rm = TRUE)
                
                # Calculate sum of scores per receiver cell type and percentage
                receiver_percentages <- positive_scores_df %>%
                    group_by(receiver) %>%
                    summarise(total_receiver_score = sum(prioritization_score, na.rm = TRUE), .groups = 'drop') %>%
                    mutate(
                        dataset = ds,
                        method = method,
                        # Calculate percentage, handle division by zero
                        percentage_score = if (total_dataset_method_score > 0) {
                                             (total_receiver_score / total_dataset_method_score) * 100
                                           } else {
                                             0 # If total score is 0, percentage is 0
                                           }
                    ) %>%
                    select(dataset, method, receiver, percentage_score)
                    
                all_percentages_list[[length(all_percentages_list) + 1]] <- receiver_percentages
                
            } else {
                 message(paste("No positive scores found for Method '", method, "' in dataset:", ds))
                 # We still need to represent this method/dataset combo, potentially with 0s.
                 # This will be handled by the `complete` function later.
            }
        } else {
            warning(paste("Invalid or missing data structure for Method '", method, "' in dataset:", ds))
        }
    }
}

# Combine all percentage data frames into one
if (length(all_percentages_list) > 0) {
    combined_percentages_df <- bind_rows(all_percentages_list)
} else {
    stop("No valid percentage data could be calculated. Check input data and filtering steps.")
}

print("Filtering for Top 2 receivers per Dataset/Method...")
top2_percentages_df <- combined_percentages_df %>%
  group_by(dataset, method) %>%
  # Arrange by percentage score descending within each group
  arrange(desc(percentage_score)) %>%
  # Keep only the top 2 rows for each group
  slice_head(n = 2) %>%
  # Ungroup after slicing
  ungroup()

# Check if filtering resulted in data
if(nrow(top2_percentages_df) == 0) {
    stop("Filtering for Top 2 receivers resulted in an empty dataframe. Check initial percentages.")
}

# --- 2. Prepare Data for Plotting ---
print("Preparing data for heatmap...")

# Define a safer separator unlikely to be in names
safe_separator <- "." # Using a period instead of underscore


# Create the combined 'dataset_receiver' column using the safe separator
combined_percentages_df <- combined_percentages_df %>%
  mutate(dataset_receiver = paste(dataset, receiver, sep = safe_separator)) %>% # Changed sep="_" to sep=safe_separator
  ungroup() # Make sure it's ungrouped after calculations

# Identify all unique methods, datasets, and dataset_receiver combinations present
all_methods <- desired_method_order
# Need to re-calculate all_datasets AFTER correcting the separation later
# all_datasets <- unique(combined_percentages_df$dataset) # Incorrect placement
all_dataset_receivers <- unique(combined_percentages_df$dataset_receiver)

# Ensure all combinations of Method and Dataset_Receiver are present, filling missing with 0%
heatmap_data_full <- combined_percentages_df %>%
  complete(method = factor(desired_method_order, levels=desired_method_order),
           dataset_receiver = factor(all_dataset_receivers), # Use factor to ensure all are included
           fill = list(percentage_score = 0)) %>%
  # Re-extract dataset info after completion using the CORRECT separator
  # Note: Need to escape '.' as it's a special regex character
  tidyr::separate(dataset_receiver,
                  into = c("dataset", "receiver"),
                  sep = paste0("\\", safe_separator), # Escape the separator for regex, Changed sep="_"
                  remove = FALSE,
                  extra = "merge", # Helps handle potential edge cases if separator somehow appears
                  fill = "right") %>% # Fills with NA if separation fails, useful for debugging
  # Ensure dataset factor is ordered (e.g., alphabetically) AFTER correct separation
  # Filter out potential rows where separation might have failed (resulting in NA dataset) before creating factor
  filter(!is.na(dataset)) %>%
  mutate(dataset = factor(dataset, levels = sort(unique(dataset))))

# Check if heatmap_data_full is empty after filtering NAs
if(nrow(heatmap_data_full) == 0) {
    stop("No valid data remained after separating dataset and receiver. Check separator and input names.")
}
# --- 3. Create Factor Levels for Ordering Axes ---

# Order methods (Y-axis): Use the predefined order
y_axis_order <- desired_method_order

# Order dataset_receiver (X-axis): Group by dataset first (using the factor level created above), then alphabetically by receiver within dataset
x_axis_order <- unique(heatmap_data_full[order(heatmap_data_full$dataset, heatmap_data_full$dataset_receiver), "dataset_receiver"])

# Apply factor levels
heatmap_data_full$method <- factor(heatmap_data_full$method, levels = rev(y_axis_order)) # Reverse for typical heatmap display (top row is first element)
heatmap_data_full$dataset_receiver <- factor(heatmap_data_full$dataset_receiver, levels = x_axis_order$dataset_receiver)
# Dataset factor is already set

# --- 4. Create Plots (Dataset Strip and Heatmap) ---
print("Generating plot components...")

# Dataset strip data - needs to be based on the *final* heatmap_data_full factors
dataset_strip_data <- heatmap_data_full %>%
  select(dataset_receiver, dataset) %>%
  distinct() %>%
  # Ensure factors match the main plot for alignment
  mutate(
      # dataset_receiver factor levels should be correct from Step 3
      dataset_receiver = factor(dataset_receiver, levels = levels(heatmap_data_full$dataset_receiver)),
      # dataset factor levels should be correct from the mutate() call above
      dataset = factor(dataset, levels = levels(heatmap_data_full$dataset))
      )

# Check if dataset_strip_data has valid levels
if(length(levels(dataset_strip_data$dataset)) == 0 || length(levels(dataset_strip_data$dataset_receiver)) == 0) {
    warning("Factor levels for plotting are empty or invalid. Check data preparation steps.")
}

# A. Dataset Strip Plot (Top Bar)
# Define a color palette for datasets - using viridis discrete scale or define manually if needed
num_datasets <- length(levels(dataset_strip_data$dataset))
# dataset_colors <- viridis::viridis_pal(option = "D")(num_datasets) # Example using viridis
# names(dataset_colors) <- levels(dataset_strip_data$dataset)

strip_plot <- ggplot(dataset_strip_data, aes(x = dataset_receiver, y = "Dataset", fill = dataset)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_brewer(palette = "Set3", name = "Dataset", drop = FALSE) +
  labs(fill = "Dataset") +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 1, b = 1, l = 10, unit = "pt")
  )


# B. Main Heatmap Plot
heatmap_plot <- ggplot(heatmap_data_full, aes(x = dataset_receiver, y = method, fill = percentage_score)) +
  geom_tile(color = "grey90", size = 0.1) + # Faint lines between tiles
  # Use a sequential color scale suitable for percentages (0-100)
  scale_fill_gradient(
    name = "% Total Score\n(Positive Only)",
    low = "white",
    high = "tomato",
    limits = c(0, 100),
    oob = scales::squish
  ) +
  # Or use a simpler gradient:
  # scale_fill_gradient(name = "% Total Score\n(Positive Only)", low = "white", high = "red", limits = c(0, 100)) +
  labs(x = "Dataset - Receiving Cell Type", y = "Method") + # Add axis titles
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.7)), # Smaller, rotated x labels
    axis.text.y = element_text(size = rel(0.9)),
    panel.grid = element_blank(), # No grid lines
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.4, "cm"),
    axis.title.x = element_text(margin = margin(t = 10)), # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)) # Add space next to y-axis title
  )

# --- 5. Combine Plots using Patchwork ---
print("Combining plots...")
# Ensure the left margins align - might need manual adjustment via theme(plot.margin) if y-axis labels differ significantly in width
# This attempts to align based on plot area but sometimes needs tweaking.
combined_signal_heatmap <- strip_plot / heatmap_plot +
  plot_layout(heights = c(1, 15), # Adjust height ratio (more space for heatmap)
              guides = "collect") & # Collect legends at bottom
  theme(legend.position = "bottom",
        legend.box = "horizontal", # Arrange legends side-by-side if multiple
        legend.box.margin = margin(t = 15) # Add margin above the collected legend
       )

# Print the combined plot to the viewer
print(combined_signal_heatmap)

# --- 6. Save Plot ---
print("Saving combined heatmap...")
output_dir <- "figures" # Or specify another directory
dir.create(output_dir, showWarnings = FALSE) # Create directory if it doesn't exist
output_filename <- file.path(output_dir, "method_receiver_signal_percentage_heatmap.png")

# Adjust width/height based on the number of columns (dataset_receivers) and rows (methods)
# These are estimates, you might need to fine-tune them
heatmap_width <- max(10, length(levels(heatmap_data_full$dataset_receiver)) * 0.15 + 2) # Base width + width per column + axis label space
heatmap_height <- max(6, length(levels(heatmap_data_full$method)) * 0.3 + 3) # Base height + height per row + title/legend space

save_result <- tryCatch({
    ggsave(output_filename, plot = combined_signal_heatmap, width = heatmap_width, height = heatmap_height, dpi = 300, limitsize = FALSE)
    TRUE
}, error = function(e) {
    warning(paste("Failed to save the signal percentage heatmap:", e$message), call. = FALSE)
    FALSE
})

if (save_result) {
    print(paste("Signal percentage heatmap saved to", output_filename))
} else {
    print("Signal percentage heatmap was generated but could not be saved automatically.")
}


#############
#v2
#############
# --- Load Necessary Libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- Ensure 'results_preprocessed' list exists ---
if (!exists("results_preprocessed") || !is.list(results_preprocessed) || length(results_preprocessed) == 0) {
    stop("The 'results_preprocessed' list is empty, not a list, or does not exist. Please ensure the initial data loading script has run successfully.")
}

# --- Define Desired Method Order ---
desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# --- 1. Aggregate Data and Calculate Percentages ---
# --- (This section remains unchanged - calculate percentages for ALL receivers first) ---
print("Aggregating scores and calculating percentages...")
all_percentages_list <- list()
processed_datasets <- names(results_preprocessed)

for (ds in processed_datasets) {
    dataset_results <- results_preprocessed[[ds]]
    if (is.null(dataset_results) || length(dataset_results) == 0) {
        warning(paste("No results found for dataset:", ds, "- skipping."))
        next
    }
    available_methods <- names(dataset_results)
    for (method in desired_method_order) {
        if (!method %in% available_methods) { next }
        method_data <- dataset_results[[method]]
        if (!is.null(method_data) && is.data.frame(method_data) && nrow(method_data) > 0 && all(c("receiver", "prioritization_score") %in% names(method_data))) {
            positive_scores_df <- method_data %>% filter(prioritization_score > 0)
            if (nrow(positive_scores_df) > 0) {
                total_dataset_method_score <- sum(positive_scores_df$prioritization_score, na.rm = TRUE)
                receiver_percentages <- positive_scores_df %>%
                    # >>> Added explicit filter for NA/empty receivers here too <<<
                    filter(!is.na(receiver), nzchar(trimws(as.character(receiver)))) %>%
                    group_by(receiver) %>%
                    summarise(total_receiver_score = sum(prioritization_score, na.rm = TRUE), .groups = 'drop') %>%
                    mutate(
                        dataset = ds,
                        method = method,
                        percentage_score = if (total_dataset_method_score > 0) {
                                             (total_receiver_score / total_dataset_method_score) * 100
                                           } else { 0 }
                    ) %>%
                    select(dataset, method, receiver, percentage_score)
                # Only add if receiver_percentages is not empty after filtering NAs
                if(nrow(receiver_percentages) > 0) {
                   all_percentages_list[[length(all_percentages_list) + 1]] <- receiver_percentages
                }
            } else { message(paste("No positive scores found for Method '", method, "' in dataset:", ds)) }
        } else { warning(paste("Invalid or missing data structure for Method '", method, "' in dataset:", ds)) }
    }
}

# Combine all percentage data frames into one
if (length(all_percentages_list) > 0) {
    combined_percentages_df <- bind_rows(all_percentages_list)
} else {
    stop("No valid percentage data could be calculated (after NA filter). Check input data and filtering steps.")
}
# --- END OF STEP 1 ---


# --- ADDED STEP: Filter for Top 2 Receivers per Dataset/Method ---
print("Filtering for Top 2 receivers per Dataset/Method...")
top2_percentages_df <- combined_percentages_df %>%
  group_by(dataset, method) %>%
  # Arrange by percentage score descending within each group
  arrange(desc(percentage_score)) %>%
  # Keep only the top 2 rows for each group
  slice_head(n = 2) %>%
  # Ungroup after slicing
  ungroup()

# Check if filtering resulted in data
if(nrow(top2_percentages_df) == 0) {
    stop("Filtering for Top 2 receivers resulted in an empty dataframe. Check initial percentages.")
}
# --- End of Added Step ---


# --- 2. Prepare Data for Plotting (using the filtered data) ---
# ** IMPORTANT: Use top2_percentages_df instead of combined_percentages_df **
print("Preparing data for heatmap (using Top 2 data)...")

safe_separator <- "." # Using a period instead of underscore

# Create the combined 'dataset_receiver' column using the safe separator
# Using the filtered dataframe 'top2_percentages_df'
# Add trimming here as well for safety
heatmap_data_to_complete <- top2_percentages_df %>%
  mutate(
    dataset = trimws(as.character(dataset)),
    receiver = trimws(as.character(receiver)),
    dataset_receiver = paste(dataset, receiver, sep = safe_separator)
    ) %>%
  filter(dataset_receiver != safe_separator) %>% # Ensure paste worked
  ungroup()

# Identify unique methods, datasets, and TOP 2 dataset_receiver combinations present
# Note: desired_method_order is still used to ensure all methods are included in the grid
all_top2_dataset_receivers <- unique(heatmap_data_to_complete$dataset_receiver)

# Check if any top receivers were identified
if(length(all_top2_dataset_receivers) == 0) {
    stop("No unique dataset-receiver combinations found after filtering for top 2.")
}

# Ensure all combinations of DESIRED Methods and the SUBSET of TOP 2 Dataset_Receivers are present,
# filling missing percentage scores with 0.
heatmap_data_full <- heatmap_data_to_complete %>%
  complete(method = factor(desired_method_order, levels=desired_method_order),
           # Complete using only the dataset_receiver values that made the top 2 cut somewhere
           dataset_receiver = factor(all_top2_dataset_receivers),
           fill = list(percentage_score = 0)) %>%
  # Re-extract dataset info after completion using the CORRECT separator
  tidyr::separate(dataset_receiver,
                  into = c("dataset", "receiver"),
                  sep = paste0("\\", safe_separator),
                  remove = FALSE,
                  extra = "merge",
                  fill = "right") %>%
  # Filter any rows where separation might have failed or resulted in NAs
  # And trim again just to be absolutely sure before factor creation
  mutate(
      dataset = trimws(as.character(dataset)),
      receiver = trimws(as.character(receiver))
      ) %>%
  filter(!is.na(dataset) & dataset != "", !is.na(receiver) & receiver != "") %>%
  # Create dataset factor AFTER cleaning and separation
  mutate(dataset = factor(dataset, levels = sort(unique(dataset))))

# Check if heatmap_data_full is empty after filtering NAs
if(nrow(heatmap_data_full) == 0) {
    stop("No valid data remained after separating dataset and receiver (Top 2 filtering). Check separator and input names.")
}

# --- 3. Create Factor Levels for Ordering Axes ---
# --- (This section remains the same, but operates on the reduced 'heatmap_data_full') ---
print("Creating factor levels for Top 2 heatmap...")
y_axis_order <- desired_method_order

# Order dataset_receiver (X-axis): Group by dataset factor levels first, then alphabetically by receiver string
x_axis_order <- heatmap_data_full %>%
    arrange(dataset, receiver) %>% # Sort by dataset factor, then receiver name
    pull(dataset_receiver) %>%     # Get the sorted unique names
    unique()                       # Keep only unique values in that order

# Apply factor levels
heatmap_data_full$method <- factor(heatmap_data_full$method, levels = rev(y_axis_order))
heatmap_data_full$dataset_receiver <- factor(heatmap_data_full$dataset_receiver, levels = x_axis_order)

# Check for NAs introduced by factor assignment
if(anyNA(heatmap_data_full$dataset_receiver)) {
   warning("NAs introduced into dataset_receiver during factor assignment!")
}


# --- 4. Create Plots (Dataset Strip and Heatmap) ---
# --- (This section remains unchanged - it plots the reduced data) ---
print("Generating plot components (Top 2)...")

# Dataset strip data - needs to be based on the *final* heatmap_data_full factors
dataset_strip_data <- heatmap_data_full %>%
  select(dataset_receiver, dataset) %>%
  distinct() %>%
  mutate(
      dataset_receiver = factor(dataset_receiver, levels = levels(heatmap_data_full$dataset_receiver)),
      dataset = factor(dataset, levels = levels(heatmap_data_full$dataset))
      )

# Check levels
if(length(levels(dataset_strip_data$dataset)) == 0 || length(levels(dataset_strip_data$dataset_receiver)) == 0) {
    warning("Factor levels for plotting (Top 2) are empty or invalid.")
}

# A. Dataset Strip Plot
strip_plot <- ggplot(dataset_strip_data, aes(x = dataset_receiver, y = "Dataset", fill = dataset)) +
  geom_tile(width = 1, height = 1) +
  # Consider a different palette if Set3 doesn't have enough distinct colors for your datasets
  scale_fill_brewer(palette = "Set3", name = "Dataset", drop = FALSE) + # Or scale_fill_viridis_d(...)
  labs(fill = "Dataset") +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 1, b = 1, l = 10, unit = "pt") # Adjust left margin if needed
  )


# B. Main Heatmap Plot (using the filtered heatmap_data_full)
heatmap_plot <- ggplot(heatmap_data_full, aes(x = dataset_receiver, y = method, fill = percentage_score)) +
  geom_tile(color = "grey90", linewidth = 0.1) + # Use linewidth instead of size
  scale_fill_gradient(
    name = "% Total Score\n(Top 2 Receivers)", # Updated legend title
    low = "white",
    high = "tomato", # Or another color like "steelblue", "forestgreen"
    limits = c(0, 100),
    oob = scales::squish
  ) +
  labs(x = "Dataset - Top Receiving Cell Type", y = "Method") + # Updated x-axis label
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.7)),
    axis.text.y = element_text(size = rel(0.9)),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.4, "cm"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# --- 5. Combine Plots using Patchwork ---
# --- (This section remains unchanged) ---
print("Combining plots (Top 2)...")
combined_signal_heatmap <- strip_plot / heatmap_plot +
  plot_layout(heights = c(1, 15), guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.margin = margin(t = 15))

print(combined_signal_heatmap)

# --- 6. Save Plot ---
# --- (This section remains unchanged, but consider a new filename) ---
print("Saving combined heatmap (Top 2)...")
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)
# Suggest a new filename to avoid overwriting the full heatmap
output_filename <- file.path(output_dir, "method_receiver_signal_percentage_heatmap_TOP2.png")

heatmap_width <- max(8, length(levels(heatmap_data_full$dataset_receiver)) * 0.2 + 2) # Adjust multiplier if needed
heatmap_height <- max(6, length(levels(heatmap_data_full$method)) * 0.3 + 3)
heatmap_width <- 20
heatmap_height <- 10
save_result <- tryCatch({
    ggsave(output_filename, plot = combined_signal_heatmap, width = heatmap_width, height = heatmap_height, dpi = 300, limitsize = FALSE)
    TRUE
}, error = function(e) {
    warning(paste("Failed to save the Top 2 signal percentage heatmap:", e$message), call. = FALSE)
    FALSE
})

if (save_result) {
    print(paste("Top 2 signal percentage heatmap saved to", output_filename))
} else {
    print("Top 2 signal percentage heatmap was generated but could not be saved automatically.")
}





















################################
##### Old Code #################
################################

######### Ridge Plots ############

# --- Generate the Plot ---
plot_ridgeline_overlaid_datasets_custom <- ggplot(combined_scores_df,
                                              aes(x = prioritization_score,
                                                  y = method, # method is now an ordered factor
                                                  fill = method,
                                                  group = interaction(method, dataset))) +
  ggridges::geom_density_ridges(alpha = 0.15,
                                  scale = 1.5,
                                  rel_min_height = 0.01) +
  # Use consistent colors; ensure names in method_colors match factor levels
  # The 'name' argument sets the legend title. Ensure 'breaks' match the desired order if needed.
  scale_fill_manual(values = method_colors, name = "Method", breaks = desired_method_order) +
  scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
      breaks = c(-10, -1, 0, 1, 10, 100)
  ) +
  theme_bw(base_size = 12) +
  theme(
    # --- Modification 3: Move Legend ---
    legend.position = "none" # Set legend position to the right
  ) +
  labs(
    # --- Modification 2: Remove Subtitle ---
    title = "Overlaid Ridgeline Densities of Interaction Scores by Method",
    # subtitle = "..." # Subtitle line removed
    x = "Interaction Prioritization Score (Pseudo-Log Scale)",
    y = "Method",
    fill = "Method" # Setting legend title here also works
  ) +
  coord_cartesian(clip = "off") # Keep allowing drawing outside plot area slightly if needed

print(plot_ridgeline_overlaid_datasets_custom)

# Adjust dimensions for saving - might need less width now legend is on right
ggsave("figures/ridgeline_scores_overlaid_datasets_custom.png", plot_ridgeline_overlaid_datasets_custom, width = 10, height = 2.5) # Adjusted width/height


#### overlap box plot 
print("Generating full overlap summary plot...")
    # --- Plotting ---
    plot_overlap_full_summary <- ggplot(combined_intersection_df,
                                        # Filter out intersections with 0 counts across all datasets? Maybe not.
                                        # filter(Count > 0), # Optional: uncomment to only plot non-empty intersections
                                        aes(x = Intersection_Name, y = Count)) +
      geom_boxplot(aes(fill = factor(Degree)), outlier.shape = NA, show.legend = TRUE) + # Color by degree
      geom_jitter(width = 0.25, height = 0, alpha = 0.4, size = 1.5, shape = 16) +
      scale_fill_viridis_d(name = "Intersection Degree") + # Color scale for degree
      theme_bw(base_size = 11) + # Slightly smaller base size might be needed
      theme(
          axis.text.x = element_text(angle = 60, hjust = 1, size = 8), # Steeper angle, smaller text
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom" # Keep legend for degree coloring
      ) +
      labs(
          title = "Distribution of All Specific Interaction Overlaps (Top 100)",
          subtitle = paste("Across", length(unique(combined_intersection_df$Dataset)), "datasets"),
          x = "Method(s) Defining Intersection",
          y = "Number of L-R Pairs"
      )

    print(plot_overlap_full_summary)

    # Save the plot - likely needs to be wide
    ggsave("figures/overlap_all_intersections_boxplot.png", plot = plot_overlap_full_summary, width = 16, height = 7, dpi = 300) # Increased width
    print("Full overlap plot saved to figures/overlap_all_intersections_boxplot.png")




    ##################
    ### boxplot that aligns with upset

# Create the plot
upset(
  fromExpression(upset_input_vector),  # Use the named vector of median counts
  nsets = 5,                           # Number of methods (original sets)
  order.by = "freq",                   # Order intersection bars by frequency (the median counts)
  decreasing = TRUE,                   # Show highest bars first
  mainbar.y.label = "Median Intersection Size", # Y-axis label for the main bar plot
  sets.x.label = "Total Interactions in Top 100", # X-axis label for the set size plot
  point.size = 2.8,                    # Size of points in the matrix
  line.size = 1,                       # Size of lines in the matrix
  mb.ratio = c(0.6, 0.4),              # Ratio main bar height to matrix height
  text.scale = c(intersection_size=1.5, # Adjust text sizes
                 tick_labels=1.5,
                 set_size=1.5,
                 main_bar_text=1.5,    # Use main_bar_text for y-axis label size if needed
                 sets_names=1.5        # Use sets_names for set name size if needed
                ),
  # Optional: Show only intersections with median count > 0?
  # min_size = 1, # If you only want combinations with median >= 1
  # Keep empty intersections = TRUE/FALSE? Default usually shows non-empty.
  # Optional: Hide the set size plot numbers if confusing (shows N=100)
  # show.numbers = FALSE,
  # Optional: Color specific intersections using queries
   queries = list(
      list(query = intersects, params = list("Decipher"), color = "orange", active = T),
      list(query = intersects, params = list("Decipher", "LIANA+"), color = "orange", active = T),
      list(query = intersects, params = list("Connectome", "Decipher"), color = "orange", active = T)
   )
)

# Close the PNG device
dev.off()

print("UpSet plot saved to figures/overlap_upset_median_plot.png")




############################
##### now a box plot that aligns with the upset 
############################
# --- Load Necessary Library ---
# Ensure ggplot2 is loaded. Stop if not available.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This script requires the 'ggplot2' package. Please install it using install.packages('ggplot2')", call. = FALSE)
}
library(ggplot2)

# --- Load or Ensure 'combined_intersection_df' is available ---
# (This part remains the same - assumes data is loaded/calculated earlier)
if (!exists("combined_intersection_df")) {
  rds_file_unfiltered <- "figures/overlap_all_intersections_data.rds"
  if (file.exists(rds_file_unfiltered)) {
    print(paste("Loading unfiltered data from:", rds_file_unfiltered))
    combined_intersection_df <- readRDS(rds_file_unfiltered)
  } else {
    stop("Required dataframe 'combined_intersection_df' not found. Ensure it is calculated or loaded.")
  }
} else {
  print("Using existing 'combined_intersection_df'.")
}

# --- Ensure necessary columns exist and calculate Degree if needed (Base R) ---
if (!"Intersection_Name" %in% names(combined_intersection_df) || !"Count" %in% names(combined_intersection_df)) {
    stop("The 'combined_intersection_df' dataframe must contain 'Intersection_Name' and 'Count' columns.")
}

if (!"Degree" %in% names(combined_intersection_df)) {
    print("Calculating 'Degree' column using base R...")
    # Ensure Intersection_Name is character before splitting
    combined_intersection_df$Intersection_Name <- as.character(combined_intersection_df$Intersection_Name)
    # Calculate degree
    combined_intersection_df$Degree <- sapply(strsplit(combined_intersection_df$Intersection_Name, " & "), length)
} else {
     print("'Degree' column found.")
}
# Flag to check if we should attempt coloring by Degree
use_degree_color <- TRUE # Assume yes initially, check later if column exists after potential filtering


# --- 1. Calculate Median Counts to Determine Order (Base R) ---
print("Calculating median counts using base R aggregate...")

# Ensure Intersection_Name is character for aggregation
if(is.factor(combined_intersection_df$Intersection_Name)) {
    combined_intersection_df$Intersection_Name <- as.character(combined_intersection_df$Intersection_Name)
}

# Use aggregate to calculate median count per intersection
# Handling potential errors if data is empty or columns missing
median_counts_agg <- tryCatch({
    aggregate(Count ~ Intersection_Name,
              data = combined_intersection_df,
              FUN = median,
              na.rm = TRUE)
}, error = function(e) {
    stop("Error during median calculation with aggregate(): ", e$message)
    return(NULL) # Return NULL if aggregate fails
})

if (is.null(median_counts_agg) || nrow(median_counts_agg) == 0) {
    stop("Median count calculation resulted in no data.")
}


# --- 2. Determine Order (Base R) ---
print("Determining plot order based on median counts...")
# Order the aggregated results by median_count descending
# The column name for median count in 'median_counts_agg' is 'Count'
order_indices <- order(median_counts_agg$Count, decreasing = TRUE)
median_counts_ordered <- median_counts_agg[order_indices, ]

# Get the specific order of intersection names
plot_order <- median_counts_ordered$Intersection_Name


# --- Optional: Filter data to match UpSetR default (median > 0) ---
# Comment out the next 4 lines if you want boxplots for *all* intersections, including those with median 0.
 print("Optional: Filtering intersections to those with median count > 0...")
 median_counts_to_keep <- median_counts_ordered[median_counts_ordered$Count > 0, ]
 plot_order <- median_counts_to_keep$Intersection_Name # Update order if filtering
 combined_intersection_df_filtered <- combined_intersection_df[combined_intersection_df$Intersection_Name %in% plot_order, ]
 plot_data <- combined_intersection_df_filtered # Use filtered data for plotting
 print(paste("Filtered down to", length(plot_order), "intersection types."))
# --- End Optional Filter ---

# If the filtering section above is commented out, use the full data
if (!exists("plot_data")) {
    print("Using data for all intersection types (including median 0).")
    plot_data <- combined_intersection_df
}


# --- 3. Prepare Plotting Data (Base R Factor) ---
print("Setting factor levels for plot order...")
# Ensure the plotting column is a factor with the desired levels
plot_data$Intersection_Name <- factor(plot_data$Intersection_Name, levels = plot_order)

# Remove any rows where Intersection_Name became NA after factor conversion (shouldn't happen if logic is correct)
plot_data <- plot_data[!is.na(plot_data$Intersection_Name), ]

if(nrow(plot_data) == 0) {
    stop("No data left to plot after ordering and potential filtering.")
}

# Re-check if Degree column exists in the final plot_data
use_degree_color <- "Degree" %in% names(plot_data)
if(use_degree_color) print("Will attempt to color boxplots by Degree.") else print("Degree column not found in final plot data, skipping color.")


# --- 4. Create the Boxplot (ggplot2) ---
print("Generating ordered ggplot2 boxplot with no x-axis labels...")

# --- Plotting ---
plot_distribution_only <- ggplot(plot_data, aes(x = Intersection_Name, y = Count)) +
  # Add geom_boxplot - conditionally add fill aesthetic
  geom_boxplot(aes(fill = if(use_degree_color) factor(Degree) else NULL), # Use fill only if Degree exists
               outlier.shape = NA,
               show.legend = use_degree_color) + # Show legend only if coloring
  # Optional: Jitter plot
  geom_jitter(width = 0.25, height = 0, alpha = 0.4, size = 1.0, shape = 16, na.rm = TRUE) +
  # Apply the specific order to the x-axis
  scale_x_discrete(limits = plot_order, drop = FALSE) + # drop=FALSE ensures all levels appear even if data is missing
  # --- Optional: Add color scale if Degree is used ---
  # Check if viridisLite is available for nicer colors, otherwise use grey or remove fill
  {
      if (use_degree_color) {
          if (requireNamespace("viridisLite", quietly = TRUE)) {
              scale_fill_viridis_d(name = "Intersection Degree")
          } else {
              warning("Package 'viridisLite' not found. Using grey scale for fill.", call. = FALSE)
              scale_fill_grey(name = "Intersection Degree", start = 0.8, end = 0.4) # Use grey scale as fallback
          }
      }
      # If not use_degree_color, this block returns NULL and no scale is added
  } +
  # --- Key Theme Modifications ---
  theme_bw(base_size = 11) +
  theme(
    # Remove X axis text, ticks, and title
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    # Optional: Adjust other theme elements if needed
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = if(use_degree_color) "right" else "none" # Adjust legend position
  ) +
  # --- Labels ---
  labs(
    # title = "Distribution of Specific Interaction Counts", # Optional title
    y = "Number of L-R Pairs" # Keep Y axis label
  )

# It's good practice to explicitly print the ggplot object
print(plot_distribution_only)

# --- 5. Save the Plot (ggplot2) ---
# Adjust filename and dimensions as needed
# Enclose ggsave in tryCatch in case of file permission issues etc.
print("Saving the plot...")
save_result <- tryCatch({
    ggsave("figures/overlap_distribution_boxplot_no_xaxis.png", plot = plot_distribution_only, width = 10, height = 2, dpi = 300)
    TRUE # Indicate success
}, error = function(e) {
    warning("Failed to save the plot: ", e$message, call. = FALSE)
    FALSE # Indicate failure
})

if (save_result) {
    print("Ordered boxplot saved to figures/overlap_distribution_boxplot_no_xaxis.png")
} else {
    print("Plot was generated but could not be saved automatically.")
}


############ V2
# --- Load Necessary Library ---
# Ensure ggplot2 is loaded. Stop if not available.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This script requires the 'ggplot2' package. Please install it using install.packages('ggplot2')", call. = FALSE)
}
library(ggplot2)

# --- Load or Ensure 'combined_intersection_df' is available ---
# (This part remains the same - assumes data is loaded/calculated earlier)
if (!exists("combined_intersection_df")) {
  rds_file_unfiltered <- "figures/overlap_all_intersections_data.rds"
  if (file.exists(rds_file_unfiltered)) {
    print(paste("Loading unfiltered data from:", rds_file_unfiltered))
    combined_intersection_df <- readRDS(rds_file_unfiltered)
  } else {
    stop("Required dataframe 'combined_intersection_df' not found. Ensure it is calculated or loaded.")
  }
} else {
  print("Using existing 'combined_intersection_df'.")
}

# --- Ensure necessary columns exist ---
# (Removed Degree calculation as it's not used for fill anymore)
if (!"Intersection_Name" %in% names(combined_intersection_df) || !"Count" %in% names(combined_intersection_df)) {
    stop("The 'combined_intersection_df' dataframe must contain 'Intersection_Name' and 'Count' columns.")
}


# --- 1. Calculate Median Counts to Determine Order (Base R) ---
print("Calculating median counts using base R aggregate...")
# Ensure Intersection_Name is character for aggregation
if(is.factor(combined_intersection_df$Intersection_Name)) {
    combined_intersection_df$Intersection_Name <- as.character(combined_intersection_df$Intersection_Name)
}
median_counts_agg <- tryCatch({
    aggregate(Count ~ Intersection_Name, data = combined_intersection_df, FUN = median, na.rm = TRUE)
}, error = function(e) {
    stop("Error during median calculation with aggregate(): ", e$message); return(NULL)
})
if (is.null(median_counts_agg) || nrow(median_counts_agg) == 0) {
    stop("Median count calculation resulted in no data.")
}


# --- 2. Determine Order (Base R) ---
print("Determining plot order based on median counts...")
order_indices <- order(median_counts_agg$Count, decreasing = TRUE)
median_counts_ordered <- median_counts_agg[order_indices, ]
plot_order <- median_counts_ordered$Intersection_Name


# --- Optional: Filter data to match UpSetR default (median > 0) ---
 print("Optional: Filtering intersections to those with median count > 0...")
 median_counts_to_keep <- median_counts_ordered[median_counts_ordered$Count > 0, ]
 plot_order <- median_counts_to_keep$Intersection_Name # Update order if filtering
 combined_intersection_df_filtered <- combined_intersection_df[combined_intersection_df$Intersection_Name %in% plot_order, ]
 plot_data <- combined_intersection_df_filtered # Use filtered data for plotting
 print(paste("Filtered down to", length(plot_order), "intersection types."))
# --- End Optional Filter ---

# # If the filtering section above is commented out, use the full data
# if (!exists("plot_data")) {
#     print("Using data for all intersection types (including median 0).")
#     plot_data <- combined_intersection_df
# }


# --- 3. Prepare Plotting Data (Base R Factor) ---
print("Setting factor levels for plot order...")
plot_data$Intersection_Name <- factor(plot_data$Intersection_Name, levels = plot_order)
plot_data <- plot_data[!is.na(plot_data$Intersection_Name), ]
if(nrow(plot_data) == 0) { stop("No data left to plot after ordering and potential filtering.") }


# --- 4. Create the Boxplot (ggplot2) - Revised Style ---
print("Generating ordered ggplot2 boxplot with revised style...")

# --- Plotting ---
plot_distribution_revised <- ggplot(plot_data, aes(x = Intersection_Name, y = Count)) +

  # --- Jitter Points (Plot FIRST to be behind boxplot) ---
  geom_jitter(
      width = 0.25,       # Adjust jitter width as needed
      height = 0,
      alpha = 0.6,        # Adjust alpha for point visibility
      size = 1.8,         # Increased point size
      shape = 16,         # Solid circle shape
      color = "grey75",   # Light grey points
      na.rm = TRUE
   ) +

  # --- Boxplot (Plot SECOND to be on top) ---
  geom_boxplot(
      fill = "grey20",     # Dark grey fill for the box
      color = "grey60",    # Light grey for box outline and whiskers
      outlier.shape = NA,  # Hide default outliers (jitter points show distribution)
      show.legend = FALSE  # No legend needed for fixed fill/color
      # fatten = NULL # Default median line thickness
      # lwd = 0.5 # Default line width for box/whiskers
   ) +

  # Apply the specific order to the x-axis
  scale_x_discrete(limits = plot_order, drop = FALSE) + # drop=FALSE prevents dropping levels with no data

  # --- Theme Modifications ---
  theme_bw(base_size = 11) + # Start with a clean theme
  theme(
    # Remove X axis elements
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    # Remove All Grid Lines
    panel.grid.major.y = element_blank(), # No major horizontal lines
    panel.grid.minor.y = element_blank(), # No minor horizontal lines
    panel.grid.major.x = element_blank(), # No vertical lines

    # Optional: Adjust panel border if needed
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),

    # Ensure legend is off
    legend.position = "none"
  ) +
  # --- Labels ---
  labs(
    y = "Number of L-R Pairs" # Keep Y axis label
    # title = "Distribution of Interaction Counts" # Optional
  )

# Print the plot to the viewer/device
print(plot_distribution_revised)

# --- 5. Save the Plot (ggplot2) ---
print("Saving the revised plot...")
# Adjust filename and dimensions as needed (height might need adjustment)
save_result_revised <- tryCatch({
    # You might need to adjust the height to get the desired look
    ggsave("figures/overlap_distribution_boxplot_revised_style.png", plot = plot_distribution_revised, width = 10, height = 2.5, dpi = 300)
    TRUE # Indicate success
}, error = function(e) {
    warning("Failed to save the revised plot: ", e$message, call. = FALSE)
    FALSE # Indicate failure
})

if (save_result_revised) {
    print("Revised ordered boxplot saved to figures/overlap_distribution_boxplot_revised_style.png")
} else {
    print("Revised plot was generated but could not be saved automatically.")
}




####### V3
plot_distribution_final <- ggplot(plot_data, aes(x = Intersection_Name, y = Count)) +
  # Jitter Points (Behind)
  geom_jitter(width = 0.25, height = 0, alpha = 0.6, size = 1.8, shape = 16, color = "grey75", na.rm = TRUE) +
  # Boxplot (On Top)
  geom_boxplot(fill = "grey20", color = "grey60", outlier.shape = NA, show.legend = FALSE) +
  # Order X axis
  scale_x_discrete(limits = plot_order, drop = FALSE) +
  # Expand Y axis slightly? May help prevent points hitting bottom axis line if count is 0
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + # Optional

  # --- Theme Modifications for Final Look ---
  theme_classic(base_size = 12) + # Use theme_classic, slightly larger base font
  theme(
    # Remove X axis elements COMPLETELY (including line)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(), # Remove x-axis line

    # Y-axis Styling
    axis.title.y = element_text(
        size = rel(1.1),   # Make title slightly larger than base
        face = "bold",     # Bold title
        margin = margin(t = 0, r = 15, b = 0, l = 0) # Add space between title and axis line/labels (adjust r value)
        ),
    axis.text.y = element_text(
        size = rel(1.0),   # Base size for labels
        face = "bold"      # Bold labels
        ),
    axis.line.y = element_line(color = "black", size = 0.5), # Ensure y-axis line is visible
    axis.ticks.y = element_line(color = "black"), # Ensure y-axis ticks are visible

    # Remove Grid Lines (classic removes major, also remove minor)
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # Legend off
    legend.position = "none",

    # Ensure plot margins don't cause clipping (usually okay, but can adjust if needed)
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5) # Default is often 5.5 points
  ) +
  # --- Labels ---
  labs(
    y = "Number of L-R Pairs" # Y axis title
  )

  # --- 5. Save the Plot (ggplot2) ---
print("Saving the final styled plot...")
# Adjust filename and dimensions
save_result_final <- tryCatch({
    ggsave("figures/overlap_distribution_boxplot_final_style.png", plot = plot_distribution_final, width = 10, height = 2.5, dpi = 300)
    TRUE # Indicate success
}, error = function(e) {
    warning("Failed to save the final plot: ", e$message, call. = FALSE)
    FALSE # Indicate failure
})

if (save_result_final) {
    print("Final styled ordered boxplot saved to figures/overlap_distribution_boxplot_final_style.png")
} else {
    print("Final plot was generated but could not be saved automatically.")
}




############################
## spearman and search space boxplots

# Create the dual-axis boxplot
p <- ggplot(combined_df, aes(x = Method2)) +
  
  # Search Space (k-value) Boxplot - Red (Slightly Transparent)
  geom_boxplot(aes(y = k_value, fill = "Search Space (k-value)"), 
               color = "black", outlier.shape = NA, alpha = 0.6) +
  
  # Spearman Correlation Boxplot - Green (Slightly Transparent)
  geom_boxplot(aes(y = Spearman * max(combined_df$k_value, na.rm = TRUE), fill = "Spearman Correlation"), 
               color = "black", outlier.shape = NA, alpha = 0.6) +
  
  # Add horizontal gridlines
  geom_hline(yintercept = seq(0, max(combined_df$k_value, na.rm = TRUE), length.out = 6), 
             linetype = "dashed", color = "grey60", alpha = 0.5) +
  
  scale_fill_manual(values = c("Search Space (k-value)" = "#E69F00", "Spearman Correlation" = "#56B4E9")) +  # Color mapping
  
  scale_y_continuous(
    name = "Search Space (k-value)",
    sec.axis = sec_axis(~ . / max(combined_df$k_value, na.rm = TRUE), name = "Spearman Correlation")
  ) +
  
  labs(x = "Method Pair", title = "Comparison of Search Space & Spearman Correlation Across Method Pairs") +
  
  facet_grid(. ~ Method1, scales = "free_x", space = "free_x") +  # Hierarchical grouping
  
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(face = "bold"),  # Bold text for main categories
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "bottom"
  )

# Save the plot
ggsave("figures/combined_k_spearman_boxplot_7apr.png", plot = p, width = 12, height = 6, dpi = 300)





######################
#V2
# --- Load necessary libraries (ensure ggplot2, dplyr are loaded) ---
library(ggplot2)
library(dplyr)

# --- Assume 'combined_df' is already created and merged ---
# It should have columns like: Var1, Var2, Method1, Method2, k_value, Spearman

# --- Define colors (matching your scale_fill_manual) ---
color_k_value <- "#E69F00" # Gold/Orange
color_spearman <- "#56B4E9" # Light Blue

# --- 1. Filter out self-comparisons ---
# Ensure Method1 and Method2 columns are correctly assigned first
# If they are not present, recreate them from Var1, Var2
if (!all(c("Method1", "Method2") %in% names(combined_df))) {
   combined_df <- combined_df %>% mutate(Method1 = Var1, Method2 = Var2)
   warning("Had to recreate Method1/Method2 columns from Var1/Var2.")
}

plot_data <- combined_df %>%
  filter(Method1 != Method2) # Keep only pairs where methods are different

# Check if data remains after filtering
if(nrow(plot_data) == 0) {
  stop("Filtering removed all data. Check input 'combined_df' and filter condition.")
}


# --- Create the updated dual-axis boxplot ---
p_updated <- ggplot(plot_data, aes(x = Method2)) +

  # Search Space (k-value) Boxplot - Gold/Orange
  geom_boxplot(aes(y = k_value, fill = "Search Space (k-value)"),
               color = "black", outlier.shape = NA, alpha = 0.7) + # Slightly more alpha

  # Spearman Correlation Boxplot - Light Blue
  geom_boxplot(aes(y = Spearman * max(plot_data$k_value, na.rm = TRUE), fill = "Spearman Correlation"),
               color = "black", outlier.shape = NA, alpha = 0.7) + # Slightly more alpha

  # --- Removed manual geom_hline ---

  # Manual color scale for fill aesthetic
  scale_fill_manual(
      name = NULL, # Hide legend title for fill
      values = c("Search Space (k-value)" = color_k_value, "Spearman Correlation" = color_spearman)
      ) +

  # Dual Y-axis setup
  scale_y_continuous(
    # --- Primary Y Axis (k-value) ---
    name = "Search Space (k-value)",
    expand = expansion(mult = c(0.05, 0.05)), # Add a little space at limits

    # --- Secondary Y Axis (Spearman) ---
    sec.axis = sec_axis(
        ~ . / max(plot_data$k_value, na.rm = TRUE), # Rescale function
        name = "Spearman Correlation"
        )
  ) +

  # Labels
  labs(
      x = "Method Compared To", # More descriptive x-axis label
      title = "Comparison of Search Space & Spearman Correlation Across Method Pairs",
      fill = NULL # Ensure fill legend title is NULL
      ) +

  # Faceting by Method1
  facet_grid(. ~ Method1, scales = "free_x", space = "free_x") +

  # Theme modifications
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), # Ensure plot background is white

    # --- Remove Grid Lines & Add Demarcations ---
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", size = 0.4, linetype = "dotted"), # Subtle vertical lines at ticks
    panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0.75, "lines"), # Increase space BETWEEN facets (Method1 groups)

    # Strip text (Method1 labels)
    strip.text.x = element_text(face = "bold", size=rel(1.1)), # Slightly larger bold text

    # X-axis text (Method2 labels)
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate labels, adjust vertical alignment

    # --- Color Y Axes ---
    # Left Y Axis (k-value)
    axis.title.y.left = element_text(color = color_k_value, face = "bold", size=rel(1.0)),
    axis.text.y.left = element_text(color = color_k_value, face="bold"),
    axis.ticks.y.left = element_line(color = color_k_value),
    axis.line.y.left = element_line(color = color_k_value), # Add axis line

    # Right Y Axis (Spearman)
    axis.title.y.right = element_text(color = color_spearman, face = "bold", size=rel(1.0)),
    axis.text.y.right = element_text(color = color_spearman, face="bold"),
    axis.ticks.y.right = element_line(color = color_spearman),
    axis.line.y.right = element_line(color = color_spearman), # Add axis line

    # Legend
    legend.position = "bottom"
  )

# Save the plot
ggsave("figures/combined_k_spearman_boxplot_updated.png", plot = p_updated, width = 14, height = 7, dpi = 300) # Adjusted dimensions slightly

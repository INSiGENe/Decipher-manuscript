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
  "SevCOVID" = "results/SevCOVID",
  "MilCOVID" = "results/MilCOVID",
  "SevCOVID_Azimuthl1" = "results/SevCOVID_Azimuthl1",
  "MilCOVID_Azimuthl1" = "results/MilCOVID_Azimuthl1",
  "SevCOVID_Azimuthl2" = "results/SevCOVID_Azimuthl2",
  "MilCOVID_Azimuthl2" = "results/MilCOVID_Azimuthl2"
)

#"cz_ra_pbmc" = "results/cz_ra_pbmc",
#"cz_cz_human_kidney_v1.5" = "results/cz_human_kidney_v1.5"
#"cz_afib_macrophages"	= "results/cz_afib_macrophages"
#"cz_placenta" = "results/cz_placenta_infection"
#"cz_dev_gut_crohns" = "results/cz_dev_gut_crohns"

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


#interaction thresholds
library(ggplot2)
library(tidyr)

analyze_method_thresholds <- function(results_list, method_name) {
  

  # Initialize a list to store summary statistics
  summary_stats <- list()
  
  # Define standard deviation thresholds from 0.1 to 4.0 in increments of 0.1
  std_thresholds <- seq(0.1, 4.0, by = 0.1)
  
  # Process each dataset in results_list
  for (dataset in names(results_list)) {
    
    # Check if the specified method exists in the dataset
    if (!method_name %in% names(results_list[[dataset]])) {
      warning(paste("Skipping dataset", dataset, "- method", method_name, "not found"))
      next
    }
    
    # Extract the method's data frame
    df <- results_list[[dataset]][[method_name]]
    
    # Check if prioritization_score exists in the dataset
    if (!"prioritization_score" %in% colnames(df)) {
      warning(paste("Skipping dataset", dataset, "- no prioritization_score column found"))
      next
    }
    
    # Compute the standard deviation of prioritization_score
    std_dev <- sd(df$prioritization_score, na.rm = TRUE)
    
    # Count total interactions (rows) and set its std_dev to 0
    total_interactions <- nrow(df)
    threshold_counts <- sapply(std_thresholds, function(thresh) {
      sum(abs(df$prioritization_score) > (thresh * std_dev), na.rm = TRUE)
    })
    
    # Store results in a named list
    summary_stats[[dataset]] <- c("> 0σ" = total_interactions, 
                                  setNames(threshold_counts, paste0("> ", std_thresholds, "σ")))
  }

  
  
  # Convert to a data frame for better visualization
  summary_df <- do.call(rbind, summary_stats)
  
  # Print the summary table
  print(summary_df)
  
  # Convert summary_df to long format for ggplot
  summary_long <- as.data.frame(summary_df)
  summary_long$Dataset <- rownames(summary_long)
  summary_long <- gather(summary_long, key = "Threshold", value = "Count", -Dataset)
  
  # Ensure thresholds are ordered numerically
  summary_long$Threshold <- as.numeric(gsub("> ", "", gsub("σ", "", summary_long$Threshold)))
  
  # Set "Total Interactions" threshold to 0 so it appears in the plot
  summary_long$Threshold[summary_long$Threshold == "Total Interactions"] <- 0
  
  # Create the boxplot
  # Get global y-axis limit
    #max_y <- max(summary_long$Count, na.rm = TRUE)
    max_y <- 100000

    # Create the boxplot with a fixed y-axis limit
    p <- ggplot(summary_long, aes(x = factor(Threshold), y = Count)) +
    geom_boxplot(fill = "grey70", color = "black", outlier.shape = NA) +  # Single grey color
    scale_y_log10(limits = c(1, max_y)) +  # Fixed log scale
    labs(x = "Threshold (σ)", y = "Number of Interactions", 
        title = paste("Distribution of Interactions Across Thresholds for", method_name)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
  
  output_filename = paste("interaction_thresholds_boxplot_",method_name,".png",sep="")
  # Save as PNG
  ggsave(output_filename, plot = p, width = 8, height = 6, dpi = 300)
  
  # Return the summary data frame (optional)
  return(summary_df)
}



analyze_method_thresholds(results_preprocessed,"NicheNet")
analyze_method_thresholds(results_preprocessed,"Decipher")
analyze_method_thresholds(results_preprocessed,"LIANA+")
analyze_method_thresholds(results_preprocessed,"NATMI")
analyze_method_thresholds(results_preprocessed,"Connectome")



#overlap box plot - figure X
library(ggplot2)
library(tidyr)

# Initialize an empty list to store overlap results
overlap_data <- list()

# Define n_top for overlap calculation
n_top <- 100

# Loop through each dataset and calculate overlaps
for (dataset in names(results_for_correlation)) {
  
  # Extract method-specific sets for the dataset
  liana_set <- getSet(results_for_correlation[[dataset]]$`LIANA+`, n_top)
  nichenet_set <- getSet(results_for_correlation[[dataset]]$NicheNet, n_top)
  decipher_set <- getSet(results_for_correlation[[dataset]]$Decipher, n_top)
  connectome_set <- getSet(results_for_correlation[[dataset]]$Connectome, n_top)
  natmi_set <- getSet(results_for_correlation[[dataset]]$NATMI, n_top)
  
  # List of all method-specific sets
  list_input <- list(
    "LIANA+" = liana_set,
    "NicheNet" = nichenet_set,
    "Decipher" = decipher_set,
    "Connectome" = connectome_set,
    "NATMI" = natmi_set
  )
  
  # Get overlap results
  dataset_overlap <- getInteractionOverlap(list_input)
  
  # Convert to a data frame and add dataset column
  dataset_overlap_df <- as.data.frame(dataset_overlap)
  dataset_overlap_df$Dataset <- dataset  # Track dataset origin
  
  # Store the results in the list
  overlap_data[[dataset]] <- dataset_overlap_df
}

# Combine results from all datasets
combined_overlap_df <- do.call(rbind, overlap_data)

# Ensure Overlap column is treated as a factor with proper ordering
combined_overlap_df$overlap <- as.factor(combined_overlap_df$overlap)


# Ensure Method Column is a Factor with Custom Ordering
combined_overlap_df$method <- factor(combined_overlap_df$method, 
                                     levels = c("NicheNet", "Decipher", "LIANA+", "Connectome", "NATMI"))

# Generate boxplot
p <- ggplot(combined_overlap_df, aes(x = method, y = Count, fill = overlap)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot with no outliers
  scale_fill_manual(values = c("red", "lightgreen", "green", "darkgreen", "black")) +  # Custom colors
  labs(x = "Method", y = "Overlapping LR pairs", title = "Distribution of Overlapping LR Pairs Across Methods") +
  theme_minimal(base_size = 14) + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"  # Move legend to bottom
  )

# Save the plot
ggsave("overlap_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)

saveRDS(combined_overlap_df,file = "overlap_boxplot.rds")

# Display the plot
print(p)

#spearman

library(ggplot2)
library(reshape2)   # For reshaping matrix to long format
library(gridExtra)  # For arranging multiple heatmaps in a grid

# Initialize list to store Spearman matrices
spearman_matrices <- list()

# Loop through each dataset
for (dataset in names(results_for_correlation)) {
  
  # Compute correlation & search space
  interaction_results_correlation_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(results_for_correlation[[dataset]])
  
  # Extract Spearman correlation matrix
  spearman_matrices[[dataset]] <- interaction_results_correlation_search_space$spearman
}

plot_spearman_heatmap <- function(spearman_matrix, title) {
  
  # Convert matrix to long format for ggplot
  spearman_long <- melt(spearman_matrix)
  
  # Generate heatmap
  p <- ggplot(spearman_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkgreen", na.value = "grey90") +  # Grey for empty values
    labs(x = NULL, y = NULL, fill = "Spearman", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Convert Spearman matrices into heatmap plots
heatmaps <- lapply(names(spearman_matrices), function(dataset) {
  plot_spearman_heatmap(spearman_matrices[[dataset]], title = dataset)
})

# Fill remaining spaces with empty plots if we have fewer than 9 datasets
while (length(heatmaps) < 9) {
  heatmaps <- append(heatmaps, list(ggplot() + theme_void()))
}

# Arrange plots in a 3x3 grid
grid_plot <- grid.arrange(grobs = heatmaps, nrow = 3, ncol = 3)

# Save the grid plot
ggsave("spearman_heatmap_grid.png", plot = grid_plot, width = 10, height = 10, dpi = 300)

#now boxplot of spearman
library(ggplot2)
library(reshape2)
library(dplyr)

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

# Create a boxplot with hierarchical facetting
p <- ggplot(combined_spearman_df, aes(x = Method2, y = value)) +
  geom_boxplot(fill = "blue", color = "black", outlier.shape = NA) +  # Blue for Spearman correlations
  labs(x = "Method Pair", y = "Spearman Correlation", title = "Spearman Correlation Between Methods Across Datasets") +
  facet_grid(. ~ Method1, scales = "free_x", space = "free_x") +  # Hierarchical facetting
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(face = "bold"),  # Bold text for main categories
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels for readability
  )

# Save the plot
ggsave("spearman_boxplot_hierarchical.png", plot = p, width = 12, height = 6, dpi = 300)

# Display the plot
print(p)


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

plot_k_heatmap <- function(k_matrix, title) {
  
  # Convert matrix to long format for ggplot
  k_long <- melt(k_matrix)
  
  # Generate heatmap
  p <- ggplot(k_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "pink", high = "darkred", na.value = "grey90") +  # Red scale
    labs(x = NULL, y = NULL, fill = "k-value", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Convert k-matrices into heatmap plots
k_heatmaps <- lapply(names(k_matrices), function(dataset) {
  plot_k_heatmap(k_matrices[[dataset]], title = dataset)
})

# Fill remaining spaces with empty plots if we have fewer than 9 datasets
while (length(k_heatmaps) < 9) {
  k_heatmaps <- append(k_heatmaps, list(ggplot() + theme_void()))
}

# Arrange plots in a 3x3 grid
grid_plot_k <- grid.arrange(grobs = k_heatmaps, nrow = 3, ncol = 3)

# Save the grid plot
ggsave("k_matrix_heatmap_grid.png", plot = grid_plot_k, width = 10, height = 10, dpi = 300)


#box plot of search space
library(ggplot2)
library(reshape2)
library(dplyr)

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

# Create a boxplot with ordered method pairs
p <- ggplot(combined_k_df, aes(x = factor(Method_Pair, levels = unique(Method_Pair)), y = value)) +
  geom_boxplot(fill = "red", color = "black", outlier.shape = NA) +  # Red boxplots with black outlines
  scale_y_log10() +  # Log scale for better visualization
  labs(x = "Method Pair", y = "Search Space (k-value)", title = "Distribution of Search Space Across Method Pairs") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels for readability
  )

# Save the plot
ggsave("k_matrix_boxplot.png", plot = p, width = 10, height = 6, dpi = 300)

# Display the plot
print(p)

# Split method pairs into two separate columns for hierarchical labeling
combined_k_df <- combined_k_df %>%
  mutate(Method1 = sub("-.*", "", Method_Pair),  # Extract first method
         Method2 = sub(".*-", "", Method_Pair))  # Extract second method

# Generate the boxplot with hierarchical labels
p <- ggplot(combined_k_df, aes(x = Method2, y = value)) +
  geom_boxplot(fill = "red", color = "black", outlier.shape = NA) +  # Boxplot styling
  scale_y_log10() +  # Log scale for better visualization
  labs(x = "Method Pair", y = "Search Space (k-value)", title = "Distribution of Search Space Across Method Pairs") +
  facet_grid(. ~ Method1, scales = "free_x", space = "free_x") +  # Hierarchical facetting
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(face = "bold"),  # Bold text for main categories
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels for readability
  )

# Save the plot
ggsave("k_matrix_boxplot_hierarchical.png", plot = p, width = 12, height = 6, dpi = 300)

# Display the plot
print(p)


#trying to merge -----
library(ggplot2)
library(reshape2)
library(dplyr)

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
  
  scale_fill_manual(values = c("Search Space (k-value)" = "red", "Spearman Correlation" = "green")) +  # Color mapping
  
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
ggsave("combined_k_spearman_boxplot.png", plot = p, width = 12, height = 6, dpi = 300)

dataset_path <- "results/5yr_pic_dup_2"
plotDecipherPrioritizedMap_v3(dataset_path,slice_n=6,output_filename="test_v3_2",log_transform = FALSE)
plotDecipherPrioritizedMap_v4(dataset_path,slice_n=6,output_filename="test_v4",log_transform = FALSE)



#################################
####### Load Cytosig data ############
#################################
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

#################################
####### Cytosig prep data ############
#################################

library(gplots)

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

#################################
####### Cytosig heatmap (single and with colour bar) ############
#################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

combined_mat <- do.call(cbind, lapply(names(valid), function(dataset) {
  mat <- valid[[dataset]]
  colnames(mat) <- paste(dataset, colnames(mat), sep = "_")
  return(mat)
}))

df_long <- as.data.frame(combined_mat) %>%
  tibble::rownames_to_column("ligand") %>%
  tidyr::pivot_longer(-ligand, names_to = "dataset_celltype", values_to = "z") %>%
  mutate(dataset = sub("_.*$", "", dataset_celltype))

df_long <- df_long %>%
  mutate(
    ligand = factor(ligand, levels = sort(unique(ligand))),
    dataset_celltype = factor(dataset_celltype, levels = unique(dataset_celltype))
  )

df_long <- df_long %>%
  mutate(
    formatted_label = gsub("_minus_", "-", dataset_celltype),
    formatted_label = gsub("_plus_", "+", formatted_label),
    #formatted_label = ifelse(nchar(formatted_label) > 20,
    #                         paste0(substr(formatted_label, 1, 17), "..."),
    #                         formatted_label),
    formatted_label = factor(formatted_label, levels = unique(formatted_label))
  )

# Dataset strip
dataset_strip <- df_long %>%
  distinct(dataset_celltype, dataset) %>%
  mutate(y = "Dataset")

dataset_strip <- df_long %>%
  distinct(formatted_label, dataset)

strip_plot <- ggplot(dataset_strip, aes(x = formatted_label, y = "Dataset", fill = dataset)) +
  geom_tile(width = 1) +
  labs(fill = "Dataset") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  )

heatmap_plot <- ggplot(df_long, aes(x = formatted_label, y = ligand, fill = z)) +
  geom_tile(width = 0.96,color = "grey95", size = 0.1) +  # narrower tiles
    scale_fill_gradient2(
    low = "#4B6C8C",      # steel blue
    mid = "white",
    high = "#C44E52",     # deep metallic red
    midpoint = 0,
    name = "Z-score"
)+  # cooler palette
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 7),
    axis.text.y = element_text(face = "bold", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm")  # thinner color bar
  )


library(patchwork)
combined_plot <- strip_plot / heatmap_plot +
  plot_layout(heights = c(0.4, 5), guides = "collect") &
  theme(legend.position = "bottom")



ggsave("figures/cytosig_final_heatmap.png",
       combined_plot, width = 14, height = 8, dpi = 300)

#################################
# Initialize an empty list to collect data frames for each dataset
combined_list <- list()

# Loop through each dataset
for (ds in names(all_cytosig_results)) {
  z_scores <- all_cytosig_results[[ds]]                        # Matrix: ligands as rows, clusters as columns
  cluster_sizes <- cells_per_cluster[[ds]]       # Data frame with Var1 = cluster name, Freq = n_cells

  # Convert z-scores matrix to long format
  z_long <- reshape2::melt(z_scores, varnames = c("ligand", "cluster"), value.name = "z_score")

  # Add cell count information by joining with cluster_sizes
  colnames(cluster_sizes) <- c("cluster", "n_cells")
  merged <- merge(z_long, cluster_sizes, by = "cluster")

  # Add dataset column
  merged$dataset <- ds

  # Append to list
  combined_list[[ds]] <- merged
}

# Combine all datasets into a single data frame
combined_df <- do.call(rbind, combined_list)

# Reorder columns
combined_df <- combined_df[, c("dataset", "cluster", "ligand", "z_score", "n_cells")]
library(dplyr)
library(tidyr)

# Count ligands with |z_score| > 2
ligand_counts <- combined_df %>%
  filter(abs(z_score) > 2) %>%
  group_by(dataset, cluster) %>%
  summarise(n_ligands_above_2 = n(), .groups = "drop")

# Get all dataset–cluster combinations with cell counts
all_clusters <- combined_df %>%
  distinct(dataset, cluster, n_cells)

# Merge and fill in zeros where no ligands pass the threshold
final_df <- all_clusters %>%
  left_join(ligand_counts, by = c("dataset", "cluster")) %>%
  mutate(n_ligands_above_2 = replace_na(n_ligands_above_2, 0))

library(ggplot2)

# Plot
p <- ggplot(final_df, aes(x = n_cells, y = n_ligands_above_2)) +
  geom_point(alpha = 0.7) +
  labs(
    x = "Number of Cells",
    y = "Number of Ligands with |z-score| > 2",
    title = "Ligand Activity vs Cell Abundance"
  ) +
  theme_minimal()

# Save to file
ggsave("figures/ligand_vs_cells_plot.png", plot = p, width = 7, height = 5, dpi = 300)


library(ggplot2)
library(dplyr)

# Create bins of size 1000
final_df <- final_df %>%
  mutate(cell_bin = cut(n_cells,
                        breaks = seq(0, max(n_cells, na.rm = TRUE) + 1000, by = 1000),
                        include.lowest = TRUE,
                        right = FALSE,
                        labels = paste(seq(0, max(n_cells, na.rm = TRUE), by = 1000),
                                       seq(999, max(n_cells, na.rm = TRUE) + 999, by = 1000),
                                       sep = "-")))

# Plot
p <- ggplot(final_df, aes(x = cell_bin, y = n_ligands_above_2)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    x = "Number of Cells (Binned)",
    y = "Number of Ligands with |z-score| > 2",
    title = "Ligand Activity by Cell Count Bin"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
ggsave("ligand_vs_cells_boxplot.png", plot = p, width = 8, height = 5, dpi = 300)


#################################
####### ROC Curves ############
#################################

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


#################################
####### Beeswarm ############
#################################
# Assuming your nested list is called `results`
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
  base_cols <- scales::hue_pal()(length(methods))
  names(base_cols) <- methods
  
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

# Construct the plot:
#   - Use geom_line() to connect points that belong to the same dataset.
#   - Use geom_point() for plotting the points (overlap is fine here).
p <- ggplot(results_df, aes(x = method, y = value)) +
  # Connect points from the same dataset. The 'group' aesthetic uses the 'dataset' column.
  geom_line(aes(group = dataset), color = "gray", size = 1) +
  
  # Plot the points.
  # The color is determined by interaction(method, flagged), which we control with our custom scale.
  geom_point(aes(color = interaction(method, flagged)), size = 4, alpha = 1) +
  
  # Apply the custom manual color scale.
  scale_color_manual(values = create_flag_color_scale(unique(results_df$method))) +
  
  # Add a segment for the median per method.
  stat_summary(fun = median, geom = "segment", 
               aes(xend = after_stat(x), yend = after_stat(y)), 
               size = 2.5, color = "black") +
  
  # Add a horizontal dashed line at y = 0.5.
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  
  # Label the axes.
  labs(y = "AUROC target prediction", x = NULL) +
  
  # Use a minimal theme with some tweaks.
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the plot.
ggsave("figures/beeswarm_auc_plot_points_lines.png", plot = p, width = 4, height = 6, dpi = 300)

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

#boxplot 
p <- ggplot(results_df, aes(x = method, y = value)) +
  geom_line(aes(group = dataset), color = "gray", size = 1, alpha = 0.3) +
  
  # Add boxplot layer
  geom_boxplot(outlier.shape = NA, width = 0.25, alpha = 0.4, color = "black", fill = "lightgray") +
  
  # Beeswarm points
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

ggsave("figures/boxplot_beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)


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

# Function: Create a custom color scale for flagged points
# For each method, we want two colors:
# - A darker version when flagged (flag = TRUE)
# - The regular color when not flagged (flag = FALSE)
create_flag_color_scale <- function(methods) {
  # Get one base color for each method using a hue palette
  base_cols <- scales::hue_pal()(length(methods))
  names(base_cols) <- methods
  
  # Define a helper function to darken a color
  darker <- function(col) {
    adjustcolor(col, red.f = 0.6, green.f = 0.6, blue.f = 0.6)
  }
  
  # For each method, assign a darker color when flagged and the base color when not flagged
  color_values <- unlist(lapply(methods, function(m) {
    c(darker(base_cols[m]), base_cols[m])
  }))
  
  # Names must match the values produced by interaction(method, flagged)
  # This will produce names like "MethodName.TRUE" and "MethodName.FALSE"
  names(color_values) <- unlist(lapply(methods, function(m) {
    c(paste0(m, ".TRUE"), paste0(m, ".FALSE"))
  }))
  
  return(color_values)
}

# Use the custom color scale in your ggplot code:
# 'results_df' must have columns: method, value, and flagged (with either "*" or empty)
p <- ggplot(results_df, aes(x = method, y = value)) +
  # Plot the beeswarm points with jittering based on the computed positions
  geom_beeswarm(aes(color = interaction(method, flagged)),
                size = 4, priority = "random", cex = 3, alpha = 1) +
  # Apply the manual color scale using our custom function
  scale_color_manual(values = create_flag_color_scale(unique(results_df$method))) +
  # Add a segment for the median values per method
  stat_summary(fun = median, geom = "segment", 
               aes(xend = after_stat(x), yend = after_stat(y)), 
               size = 2.5, color = "black") +
  # Add a horizontal dashed red line at y = 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  # Label the axes
  labs(y = "AUROC target prediction", x = NULL) +
  # Use a minimal theme and adjust text
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the plot to a file
ggsave("beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)



#################################
####### Old Cytosig code ########
#################################
plots <- list()

for (dataset in names(valid)) {
  mat <- valid[[dataset]]
  # Set |z| ≤ 2 values to 0
  mat[abs(mat) <= 2] <- 0
  df_long <- as.data.frame(mat) %>%
    rownames_to_column("ligand") %>%
    pivot_longer(-ligand, names_to = "cell_type", values_to = "z") 
  
  p <- ggplot(df_long, aes(x = cell_type, y = ligand, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal(base_size = 10) +
    labs(title = dataset, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  plots[[dataset]] <- p
}
combined_plot <- wrap_plots(plots, nrow = 1) +
  plot_annotation(title = "Cytosig Z-Score Heatmaps (Filtered |z| > 2)")

ggsave("figures/cytosig_heatmaps_combined_geom_tile_padded.png",
       plot = combined_plot,
       width = 4 * length(plots),
       height = 6,
       dpi = 300)



#################################
####### Circos plots ############
#################################
library(circlize)
library(tibble)

# Data with lists of variable length
df <- tibble(
  ligand_receptor = c("VCAN-CD44", "TNFSF13B-TNFSF13C", "SLAMF7-SLAMF7", "MRC-FPC9"),
  sender = list(
    "Mono",
    "B",
    c("Mono", "CD4", "NK"),
    "Mono"
  ),
  receiver = list(
    "B",
    "B",
    c("Mono", "NK"),
    c("B", "CD4", "NK")
  )
)

# Expand into long format
chords <- do.call(rbind, lapply(1:nrow(df), function(i) {
  expand.grid(
    from = df$sender[[i]],
    to = df$receiver[[i]],
    interaction = df$ligand_receptor[i],
    stringsAsFactors = FALSE
  )
}))

# Colors for interactions and cell types
interaction_types <- unique(chords$interaction)
link_colors <- setNames(rainbow(length(interaction_types)), interaction_types)

cell_types <- unique(c(chords$from, chords$to))
grid_colors <- setNames(rainbow(length(cell_types)), cell_types)

# Save to PNG
png("circos_ligand_receptor_labeled.png", width = 1000, height = 1000)
circos.clear()

# Initialize sectors
circos.par(gap.degree = 5)
chordDiagram(
  x = chords[, c("from", "to")],
  grid.col = grid_colors,
  col = link_colors[chords$interaction],
  annotationTrack = "grid",
  directional = 1,
  direction.type = c("arrows"),
  link.arr.type = "big.arrow",
  transparency = 0.3,
  preAllocateTracks = list(track.height = 0.1)
)

# Add sector labels (cell types)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(
    x = mean(get.cell.meta.data("xlim")),
    y = get.cell.meta.data("ylim")[1] + 1,
    labels = sector.name,
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5),
    cex = 1
  )
}, bg.border = NA)

# Add link labels (interaction names)
for (i in 1:nrow(chords)) {
  link <- chords[i, ]
  circos.link(
    sector.index1 = link$from,
    point1 = c(0, 1),
    sector.index2 = link$to,
    point2 = c(0, 1),
    col = link_colors[link$interaction],
    rou = 0.8,
    border = NA
  )
}

title("Ligand-Receptor Circos Plot\n(Colored by Interaction, Labeled by Cell Type)")
dev.off()



#5yr PIC
VCAN-CD44 from Mono to B
TNFSF13B-TNFSF13C rom B to B
SLAMF7-SLAMF7 from Mono, CD4 and NK to Mono and NK
MRC-FPC9 from Mono to B, CD4 and NK


#################
# experimenting with gemini
# Ensure necessary libraries are loaded
library(ggplot2)
library(dplyr)     # Or library(data.table) if you prefer
library(purrr)     # For map functions
library(stringr)   # For string manipulation if needed
library(patchwork) # For combining plots (optional)
library(ggridges)  # For density ridgeline plots

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

# Optional: Check the structure and summary
# print(head(combined_scores_df))
# print(summary(combined_scores_df))
# print(table(combined_scores_df$dataset, combined_scores_df$method))

# Optional: Remove datasets/methods with zero rows if any slipped through
combined_scores_df <- combined_scores_df %>% filter(!is.na(prioritization_score))

# --- Plotting Options ---

# Define a consistent color scheme (adjust colors as needed)
method_colors <- c(
  "Decipher" = "#1f77b4",
  "NicheNet" = "#ff7f0e",
  "LIANA+" = "#2ca02c",
  "Connectome" = "#d62728",
  "NATMI" = "#9467bd"
  # Add more if needed
)

# Option 1: Violin Plots, Faceted by Dataset
# Good for comparing methods within each specific biological context.
# Using scales = "free_y" allows each facet to have its own y-axis range,
# which is essential if score scales differ vastly between datasets.
plot_violin_by_dataset <- ggplot(combined_scores_df, aes(x = method, y = prioritization_score, fill = method)) +
  geom_violin(trim = FALSE, scale = "width") + # trim=FALSE keeps tails, scale="width" makes violins have same max width
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) + # Add boxplot inside, hide outliers as violin shows density
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ dataset, scales = "free_y", ncol = 4) + # Adjust ncol as needed
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom" # Or "none" if redundant
  ) +
  labs(
    title = "Distribution of Interaction Scores by Method (Faceted by Dataset)",
    x = "Method",
    y = "Interaction Prioritization Score",
    fill = "Method"
  )

#print(plot_violin_by_dataset)
ggsave("figures/violin_scores_by_dataset.png", plot_violin_by_dataset, width = 14, height = 10) # Adjust size

# Option 2: Density Plots, Faceted by Dataset
# Good for seeing the shape of the distribution more clearly.
# Using scales = "free" allows both x and y axes to adapt per facet.
plot_density_by_dataset <- ggplot(combined_scores_df, aes(x = prioritization_score, fill = method, color = method)) +
  geom_density(alpha = 0.5) + # Use alpha for transparency if densities overlap
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ dataset, scales = "free", ncol = 4) + # Use 'free' if x-axis ranges also vary significantly
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Density of Interaction Scores by Method (Faceted by Dataset)",
    x = "Interaction Prioritization Score",
    y = "Density",
    fill = "Method",
    color = "Method"
  )

#print(plot_density_by_dataset)
ggsave("figures/density_scores_by_dataset.png", plot_density_by_dataset, width = 14, height = 10)

# Option 3: Ridgeline Density Plots, Faceted by Dataset
# Can be very effective for comparing multiple distributions, especially shapes.
plot_ridgeline_by_dataset <- ggplot(combined_scores_df, aes(x = prioritization_score, y = method, fill = method)) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.5, rel_min_height = 0.01) + # Adjust scale for overlap, rel_min_height removes tiny tails
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ dataset, scales = "free_x", ncol = 4) + # Free x-axis is common here
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none" # Often redundant as y-axis labels methods
  ) +
  labs(
    title = "Ridgeline Density of Interaction Scores by Method (Faceted by Dataset)",
    x = "Interaction Prioritization Score",
    y = "Method",
    fill = "Method"
  ) +
  coord_cartesian(clip = "off") # Allows labels to slightly overflow if needed

#print(plot_ridgeline_by_dataset)
ggsave("figures/ridgeline_scores_by_dataset.png", plot_ridgeline_by_dataset, width = 14, height = 10)


# --- Addressing Log Scales ---
# If your scores span many orders of magnitude OR are heavily skewed,
# plotting on a log scale might be beneficial.

# Check score ranges
 summary(combined_scores_df$prioritization_score)
# If scores include zero or negative values, use a pseudo-log or add a constant.
# Example: Using pseudo-log for density plot x-axis
# Note: Requires the 'scales' package
 library(scales)
 plot_density_log_x <- ggplot(combined_scores_df, aes(x = prioritization_score, fill = method, color = method)) +
   geom_density(alpha = 0.5) +
   scale_fill_manual(values = method_colors) +
   scale_color_manual(values = method_colors) +
   scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1), # Adjust sigma as needed
                      breaks = c(0, 0.1, 1, 10, 100, max(combined_scores_df$prioritization_score, na.rm=TRUE))) + # Adjust breaks
   facet_wrap(~ dataset, scales = "free", ncol = 4) +
   theme_bw(base_size = 12) +
   theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom") +
   labs(title = "Density of Interaction Scores (Pseudo-Log Scale)")
# print(plot_density_log_x)
ggsave("figures/plot_density_log_x.png", plot_density_log_x, width = 14, height = 10)

# Example: Using log10 for violin plot y-axis (if all scores > 0)
# Need to filter scores <= 0 or add a small constant first if necessary
# combined_scores_positive <- combined_scores_df %>% filter(prioritization_score > 0)
# plot_violin_log_y <- ggplot(combined_scores_positive, aes(x = method, y = prioritization_score, fill = method)) +
#   geom_violin(trim = FALSE, scale = "width") +
#   geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) +
#   scale_fill_manual(values = method_colors) +
#   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   facet_wrap(~ dataset, scales = "free_y", ncol = 4) +
#   theme_bw(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), ...) +
#   labs(title = "Distribution of Interaction Scores (Log10 Scale)", ...)
# print(plot_violin_log_y)

# Option 1: Violin Plots, Faceted by Dataset (Linear Scale - Recommended as primary)
plot_violin_by_dataset_linear <- ggplot(combined_scores_df, aes(x = method, y = prioritization_score, fill = method)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, coef = 0) + # coef=0 shows whiskers extending to min/max
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ dataset, scales = "free_y", ncol = 4) + # *** free_y is crucial ***
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Distribution of Interaction Scores by Method (Linear Scale)",
    subtitle = "Faceted by Dataset (Note: Y-axes vary between facets)",
    x = "Method",
    y = "Interaction Prioritization Score",
    fill = "Method"
  )

#print(plot_violin_by_dataset_linear)
ggsave("figures/violin_scores_by_dataset_linear.png", plot_violin_by_dataset_linear, width = 14, height = 10)

# Option 2: Density Plots, Faceted by Dataset (Pseudo-Log Scale)
library(scales) # Make sure scales package is loaded

plot_density_by_dataset_pseudolog <- ggplot(combined_scores_df, aes(x = prioritization_score, fill = method, color = method)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  # Apply pseudo-log transformation to the x-axis
  scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.1, base = 10), # Adjust sigma if needed (controls linear part around 0)
      # Define breaks manually for better readability on pseudo-log scale
      breaks = c(-10, -1, 0, 1, 10, 100) # Adjust breaks based on your data range and interest
  ) +
  facet_wrap(~ dataset, scales = "free", ncol = 4) + # free allows density axis (y) and transformed x-axis to adapt
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) # Angle might be needed for breaks
  ) +
  labs(
    title = "Density of Interaction Scores by Method (Pseudo-Log Scale)",
    subtitle = "Faceted by Dataset (Note: Axes vary between facets)",
    x = "Interaction Prioritization Score (Pseudo-Log Scale)",
    y = "Density",
    fill = "Method",
    color = "Method"
  )

#print(plot_density_by_dataset_pseudolog)
ggsave("figures/density_scores_by_dataset_pseudolog.png", plot_density_by_dataset_pseudolog, width = 14, height = 10)

# Option 3: Ridgeline Plots (Pseudo-Log Scale) - Often very effective
plot_ridgeline_by_dataset_pseudolog <- ggplot(combined_scores_df, aes(x = prioritization_score, y = method, fill = method)) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.5, rel_min_height = 0.01) +
  scale_fill_manual(values = method_colors) +
  # Apply pseudo-log transformation to the x-axis
  scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
      breaks = c(-10, -1, 0, 1, 10, 100) # Adjust breaks
  ) +
  facet_wrap(~ dataset, scales = "free_x", ncol = 4) + # free_x allows transformed x-axis to adapt
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none"
  ) +
  labs(
    title = "Ridgeline Density of Interaction Scores by Method (Pseudo-Log Scale)",
    subtitle = "Faceted by Dataset (Note: X-axes vary between facets)",
    x = "Interaction Prioritization Score (Pseudo-Log Scale)",
    y = "Method",
    fill = "Method"
  ) +
  coord_cartesian(clip = "off")

#print(plot_ridgeline_by_dataset_pseudolog)
ggsave("figures/ridgeline_scores_by_dataset_pseudolog.png", plot_ridgeline_by_dataset_pseudolog, width = 14, height = 10)



# Ensure necessary libraries are loaded: ggplot2, dplyr, scales, ggridges

# (Assuming combined_scores_df and method_colors are already created)

plot_ridgeline_overlaid_datasets <- ggplot(combined_scores_df,
                                          aes(x = prioritization_score,
                                              y = method,
                                              # Define fill by method
                                              fill = method,
                                              # Group by the combination of method AND dataset
                                              # This ensures a separate ridge is drawn for each dataset within each method's y-position
                                              group = interaction(method, dataset))) +
  # Use lower alpha for more transparency; adjust scale if ridges overlap too much vertically
  ggridges::geom_density_ridges(alpha = 0.15,  # Significantly reduced alpha
                                  scale = 1.5,   # Adjust vertical scaling if needed (smaller value reduces overlap)
                                  rel_min_height = 0.01) + # Remove tiny tails
  scale_fill_manual(values = method_colors) + # Use consistent colors
  # Apply pseudo-log transformation to the common x-axis
  scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.1, base = 10), # Same transformation
      breaks = c(-10, -1, 0, 1, 10, 100) # Same breaks
      # You might want to manually set x-axis limits if the auto-range isn't ideal
      # limits = c(-20, 110)
  ) +
  theme_bw(base_size = 12) +
  theme(
    # Ensure legend is shown
    legend.position = "bottom" 
  ) +
  labs(
    title = "Overlaid Ridgeline Densities of Interaction Scores by Method",
    subtitle = "Each ridge represents one dataset (15 datasets overlaid per method)",
    x = "Interaction Prioritization Score (Pseudo-Log Scale)",
    y = "Method",
    fill = "Method" # Legend title
  ) +
  coord_cartesian(clip = "off") # Keep allowing drawing outside plot area slightly if needed

print(plot_ridgeline_overlaid_datasets)

# Adjust dimensions for saving - might need more height than width now
ggsave("figures/ridgeline_scores_overlaid_datasets.png", plot_ridgeline_overlaid_datasets, width = 10, height = 5.5) # Adjust as needed

#n3xt iteration
# Ensure necessary libraries are loaded: ggplot2, dplyr, scales, ggridges

# (Assuming combined_scores_df and method_colors are already created)

# --- Modification 1: Set Method Order ---
# Define the desired order for methods on the y-axis, with Decipher first
# The order here determines the order in the legend and calculation,
# but we use rev() below for plotting order.
desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# Convert the 'method' column to a factor with the specified levels.
# We use rev() because ggplot plots factors bottom-up on the y-axis.
# To have "Decipher" at the top visually, it needs to be the last level.
combined_scores_df <- combined_scores_df %>%
  mutate(method = factor(method, levels = rev(desired_method_order)))

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
    legend.position = "right" # Set legend position to the right
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



##########################
#an overlap alternative
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr) # May not be needed depending on getSet, but good practice

# --- Placeholder for getSet Function ---
# This function needs to extract the top 'n' interaction identifiers
# from the input object (which is results_for_correlation[[ds]][[method_name]]).
# You MUST adapt this function based on the actual structure of your data.
# Assumption: The input object is a dataframe-like structure, already sorted
#             by rank, with a column named "interaction" containing the L-R pair ID.
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
  interactions <- as.character(method_result[[interaction_col]])
  return(head(interactions, n))
}


# --- Helper function to calculate sizes of specific intersections ---
# Input: list_input = named list of character vectors (sets) for 5 methods
# Output: dataframe with Intersection_Name and Count
calculate_specific_intersections <- function(list_input) {
  method_names <- names(list_input)
  expected_methods <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

  # Ensure all 5 methods are potentially present, use empty set if NULL or missing
  sets <- list()
  for(m in expected_methods) {
    sets[[m]] <- if (!is.null(list_input[[m]]) && length(list_input[[m]]) > 0) list_input[[m]] else character(0)
  }

  # Assign to named variables for easier reading
  D_set <- sets[["Decipher"]]
  N_set <- sets[["NicheNet"]]
  L_set <- sets[["LIANA+"]]
  A_set <- sets[["NATMI"]] # A for NATMI
  C_set <- sets[["Connectome"]]

  # --- Set Operations ---
  # Union of all sets NOT including the target method (for calculating uniques)
  union_not_D <- Reduce(union, sets[c("NicheNet", "LIANA+", "NATMI", "Connectome")])
  union_not_N <- Reduce(union, sets[c("Decipher", "LIANA+", "NATMI", "Connectome")])
  union_not_L <- Reduce(union, sets[c("Decipher", "NicheNet", "NATMI", "Connectome")])
  union_not_A <- Reduce(union, sets[c("Decipher", "NicheNet", "LIANA+", "Connectome")])
  union_not_C <- Reduce(union, sets[c("Decipher", "NicheNet", "LIANA+", "NATMI")])

  # Union of all sets NOT including the specific pairs (for calculating pair-only)
  union_not_DL <- Reduce(union, sets[c("NicheNet", "NATMI", "Connectome")])
  union_not_CA <- Reduce(union, sets[c("Decipher", "NicheNet", "LIANA+")])

  results <- list()

  # Calculate unique to each
  results[["Unique: Decipher"]]   <- length(setdiff(D_set, union_not_D))
  results[["Unique: NicheNet"]]   <- length(setdiff(N_set, union_not_N))
  results[["Unique: LIANA+"]]     <- length(setdiff(L_set, union_not_L))
  results[["Unique: NATMI"]]      <- length(setdiff(A_set, union_not_A))
  results[["Unique: Connectome"]] <- length(setdiff(C_set, union_not_C))

  # Calculate intersection for specific pairs (ONLY these two)
  intersect_D_L <- intersect(D_set, L_set)
  results[["Shared: Decipher & LIANA+"]] <- length(setdiff(intersect_D_L, union_not_DL))

  intersect_C_A <- intersect(C_set, A_set)
  results[["Shared: Connectome & NATMI"]] <- length(setdiff(intersect_C_A, union_not_CA))

  # Calculate intersection for ALL methods
  all_5 <- Reduce(intersect, sets)
  results[["Shared: All"]] <- length(all_5)

  # Convert list to dataframe
  results_df <- data.frame(
    Intersection_Name = names(results),
    Count = unlist(results)
  )
  rownames(results_df) <- NULL # Clean up row names
  return(results_df)
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

# Define n_top for overlap calculation
n_top <- 100
# Define expected method names
method_names_expected <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")
# Initialize list to store results
all_intersection_data <- list()

# Ensure 'results_for_correlation' list is populated by your previous data loading loop
if (!exists("results_for_correlation") || length(results_for_correlation) == 0) {
    stop("The 'results_for_correlation' list was not found or is empty. ",
         "Please ensure the data loading and preprocessing loop has run successfully.")
}

print(paste("Processing", length(results_for_correlation), "datasets for overlap analysis..."))

# Loop through each dataset IN THE POPULATED results_for_correlation list
dataset_names_to_process <- names(results_for_correlation)

for (ds in dataset_names_to_process) {
    print(paste("Calculating overlaps for:", ds))
    dataset_results <- results_for_correlation[[ds]]

    # Check if dataset_results entry is valid
    if (is.null(dataset_results) || !is.list(dataset_results) || length(dataset_results) == 0) {
        warning(paste("Skipping dataset", ds, "- entry in results_for_correlation is NULL, not a list, or empty."))
        next # Skip to the next dataset
    }

    # Create the list_input (list of top N sets) for the current dataset
    list_input_current_ds <- list()
    available_methods_in_ds <- names(dataset_results)

    for (method_name in method_names_expected) {
        if (method_name %in% available_methods_in_ds) {
             # Use tryCatch for robustness when calling getSet
             list_input_current_ds[[method_name]] <- tryCatch({
                  getSet(dataset_results[[method_name]], n_top)
             }, error = function(e) {
                  warning(paste("Error in getSet for", method_name, "in dataset", ds, ":", e$message))
                  return(character(0)) # Return empty set on error
             })
        } else {
            # Method results not found for this dataset in results_for_correlation
             list_input_current_ds[[method_name]] <- character(0) # Use empty set
             # Optional: Add a warning if you expect all methods in all datasets
             # warning(paste("Method", method_name, "not found in results_for_correlation for dataset", ds))
        }
         # Sanity check the set length (optional)
         # print(paste("  ", ds, "-", method_name, "set size:", length(list_input_current_ds[[method_name]])))
    }

    # Calculate specific intersections for the current dataset
    dataset_summary <- tryCatch({
         calculate_specific_intersections(list_input_current_ds)
    }, error = function(e) {
         warning(paste("Error calculating intersections for dataset", ds, ":", e$message))
         return(NULL) # Return NULL on error
    })

    # Store results if calculation was successful
    if (!is.null(dataset_summary)) {
        dataset_summary$Dataset <- ds # Add dataset identifier
        all_intersection_data[[ds]] <- dataset_summary
    } else {
         warning(paste("No intersection summary generated for dataset:", ds))
    }
}


# --- Combine Results and Plot ---
if (length(all_intersection_data) > 0) {
    combined_intersection_df <- bind_rows(all_intersection_data)

    # Define logical order for plotting intersections on the x-axis
    intersection_order <- c(
        "Unique: Decipher", "Unique: NicheNet", "Unique: LIANA+", "Unique: NATMI", "Unique: Connectome",
        "Shared: Decipher & LIANA+", "Shared: Connectome & NATMI",
        "Shared: All"
    )

    # Convert Intersection_Name to factor with specified order
    combined_intersection_df <- combined_intersection_df %>%
        filter(!is.na(Intersection_Name)) %>% # Remove rows if name is NA
        mutate(Intersection_Name = factor(Intersection_Name, levels = intersection_order)) %>%
        filter(!is.na(Intersection_Name)) # Remove rows where factor conversion failed (name not in levels)

    if(nrow(combined_intersection_df) == 0) {
        stop("No valid intersection data remained after filtering.")
    }

    print("Generating overlap summary plot...")
    # --- Plotting ---
    plot_overlap_summary <- ggplot(combined_intersection_df, aes(x = Intersection_Name, y = Count)) +
      # Boxplot layer: fill by intersection type, hide legend as x-axis is clear
      geom_boxplot(aes(fill = Intersection_Name), outlier.shape = NA, show.legend = FALSE) +
      # Jitter layer: show individual data points (datasets)
      geom_jitter(width = 0.25, height = 0, alpha = 0.6, size = 2, shape = 16) +
      theme_bw(base_size = 12) +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Angled text for readability
          panel.grid.major.x = element_blank(), # Cleaner look
          panel.grid.minor.y = element_blank()  # Cleaner look
      ) +
      labs(
          title = "Distribution of Specific Interaction Overlaps (Top 100)",
          subtitle = paste("Across", length(unique(combined_intersection_df$Dataset)), "datasets"),
          x = "Type of Overlap",
          y = "Number of L-R Pairs"
      )

    print(plot_overlap_summary)

    # Save the plot
    ggsave("figures/overlap_specific_intersections_boxplot.png", plot = plot_overlap_summary, width = 9, height = 6, dpi = 300)
    print("Overlap plot saved to figures/overlap_specific_intersections_boxplot.png")

    # Save the combined data used for the plot
    saveRDS(combined_intersection_df, file = "figures/overlap_specific_intersections_data.rds")
    print("Overlap data saved to figures/overlap_specific_intersections_data.rds")

} else {
    warning("No intersection data was successfully calculated for any dataset. Cannot create plot.")
}


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# --- Placeholder for getSet Function ---
# Make sure this function is correctly defined based on your data structure
getSet <- function(method_result, n) {
  if (is.null(method_result)) return(character(0))
  if (!is.data.frame(method_result) && !is(method_result, "DataFrame") && !is.matrix(method_result)) {
     if(is.character(method_result)) return(head(method_result, n))
     warning("Input to getSet not recognized. Returning empty set.")
     return(character(0))
  }
   if (nrow(method_result) == 0) return(character(0))
  interaction_col <- "interaction" # <--- ADJUST IF NEEDED
  if (!interaction_col %in% colnames(method_result)) {
     warning(paste("Column '", interaction_col, "' not found. Returning empty set.", sep=""))
     return(character(0))
  }
  interactions <- as.character(method_result[[interaction_col]])
  return(head(interactions, n))
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

    # Save the combined data used for the plot
    saveRDS(combined_intersection_df, file = "figures/overlap_all_intersections_data.rds")
    print("Full overlap data saved to figures/overlap_all_intersections_data.rds")

} else {
    warning("No intersection data was successfully calculated for any dataset. Cannot create plot.")
}

# --- Add this filtering step after calculating medians and ordering ---

# (Assuming combined_intersection_df and median_counts exist from previous step)

# Define the threshold for median count
median_threshold <- 0 # Or 1, or other small number

# Filter the median counts table first
median_counts_filtered <- median_counts %>%
  filter(median_count > median_threshold)

# Get the names of intersections meeting the threshold
intersection_order_filtered <- median_counts_filtered$Intersection_Name

# Filter the main dataframe AND apply factor levels
combined_intersection_filtered_df <- combined_intersection_df %>%
  filter(Intersection_Name %in% intersection_order_filtered) %>%
  mutate(Intersection_Name = factor(Intersection_Name, levels = intersection_order_filtered)) %>%
  filter(!is.na(Intersection_Name)) # Ensure factor conversion worked

# --- Plotting (using the filtered data) ---
if(nrow(combined_intersection_filtered_df) > 0) {

  print(paste("Generating filtered overlap summary plot (Median Count >", median_threshold, ")..."))

  # Recalculate Degree if needed for coloring
  combined_intersection_filtered_df <- combined_intersection_filtered_df %>%
        mutate(Degree = sapply(strsplit(as.character(Intersection_Name), " & "), length))


  plot_overlap_filtered_summary <- ggplot(combined_intersection_filtered_df,
                                            aes(x = Intersection_Name, y = Count)) +
    geom_boxplot(aes(fill = factor(Degree)), outlier.shape = NA, show.legend = TRUE) +
    geom_jitter(width = 0.25, height = 0, alpha = 0.4, size = 1.5, shape = 16) +
    scale_fill_viridis_d(name = "Intersection Degree", guide = guide_legend(reverse = TRUE)) +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9), # Adjust angle/size as needed
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom"
    ) +
    labs(
        title = "Distribution of Major Interaction Overlaps (Top 100)",
        subtitle = paste("Across", length(unique(combined_intersection_filtered_df$Dataset)), "datasets (Intersections with Median Count >", median_threshold, ")"),
        x = "Method(s) Defining Intersection (Ordered High-to-Low by Median)",
        y = "Number of L-R Pairs"
    )

  print(plot_overlap_filtered_summary)

  # Save the filtered plot
  ggsave("figures/overlap_filtered_intersections_boxplot.png", plot = plot_overlap_filtered_summary, width = 10, height = 6, dpi = 300) # Adjust size
  print("Filtered overlap plot saved.")

} else {
   warning("No intersections met the filtering threshold. No plot generated.")
}


#############################################
##### ok let's try upsetR
# Load necessary libraries
library(dplyr)
library(UpSetR) # Make sure UpSetR is installed install.packages("UpSetR")

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
png("figures/overlap_upset_median_plot.png", width = 10, height = 6, units = "in", res = 300)

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
  text.scale = c(intersection_size=1.1, # Adjust text sizes
                 tick_labels=1.1,
                 set_size=1.1,
                 main_bar_text=1.1,    # Use main_bar_text for y-axis label size if needed
                 sets_names=1.1        # Use sets_names for set name size if needed
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

# Optional: Save the median counts data used for the plot
# saveRDS(median_counts_all, file = "figures/overlap_median_counts_data.rds")
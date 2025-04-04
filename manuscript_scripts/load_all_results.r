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
  "MilCOVID" = "results/MilCOVID"
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



ggsave("cytosig_final_heatmap.png",
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
ggsave("ligand_vs_cells_plot.png", plot = p, width = 7, height = 5, dpi = 300)


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
#library(ggbeeswarm)
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
ggsave("beeswarm_auc_plot_points_lines.png", plot = p, width = 4, height = 6, dpi = 300)

library(ggbeeswarm)  # Make sure it's installed

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

ggsave("beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)

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

ggsave("boxplot_beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)


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
ggsave("variance_w_boxplot_beeswarm_auc_plot.png", plot = p, width = 4, height = 6, dpi = 300)

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

ggsave("cytosig_heatmaps_combined_geom_tile_padded.png",
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

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the data frame
df <- data.frame(
  Comparison = c("PIC 5yr old", "PIC Cord", "BCG", "COVID", "Lupus", "Sepsis", "TNBC", "ERP"),
  Reported_D = c(679, 570, 618, 683, 158, 273, 594, 466),
  Reported_L = c(1821, 2302, 4335, 2348, 467, 1147, 2142, 2031),
  Reported_N = c(5566, 6648, 11710, 9659, 1948, 5683, 8863, 8119),
  Reported_C = c(4166, 5088, 14947, 12075, 4147, 8476, 9960, 9357),
  Reported_T = c(4166, 5088, 14947, 16899, 3818, 8476, 9960, 9357),
  `>1σ_D` = c(40, 38, 25, 41, 13, 7, 16, 11),
  `>1σ_L` = c(208, 257, 451, 380, 73, 160, 331, 342),
  `>1σ_N` = c(1024, 1179, 2246, 1661, 339, 1006, 1548, 1641),
  `>1σ_C` = c(390, 497, 1961, 1318, 474, 1011, 1065, 1082),
  `>1σ_T` = c(681, 1052, 1869, 2634, 677, 1739, 1714, 1641),
  `>2σ_D` = c(20, 25, 10, 16, 3, 3, 8, 3),
  `>2σ_L` = c(73, 51, 148, 27, 1, 15, 59, 43),
  `>2σ_N` = c(108, 166, 214, 330, 47, 159, 125, 93),
  `>2σ_C` = c(166, 239, 830, 544, 213, 406, 460, 452),
  `>2σ_T` = c(5, 35, 114, 171, 93, 200, 161, 133),
  `>3σ_D` = c(9, 14, 6, 6, 1, 0, 5, 2),
  `>3σ_L` = c(36, 25, 38, 0, 0, 0, 8, 5),
  `>3σ_N` = c(0, 4, 0, 5, 0, 2, 0, 0),
  `>3σ_C` = c(76, 122, 312, 256, 93, 183, 218, 191),
  `>3σ_T` = c(0, 0, 6, 10, 3, 7, 11, 3)
)

# Reshape data to long format
df_long <- df %>%
  pivot_longer(cols = -Comparison, names_to = "Metric", values_to = "Interactions") %>%
  mutate(
    Method = case_when(
      grepl("_D$", Metric) ~ "D",
      grepl("_L$", Metric) ~ "L",
      grepl("_N$", Metric) ~ "N",
      grepl("_C$", Metric) ~ "C",
      grepl("_T$", Metric) ~ "T"
    ),
    Significance = case_when(
      grepl("^Reported", Metric) ~ "Total",
      grepl("^>1σ", Metric) ~ ">1σ",
      grepl("^>2σ", Metric) ~ ">2σ",
      grepl("^>3σ", Metric) ~ ">3σ"
    )
  )

# Plot with ggplot2
ggplot(df_long, aes(x = Method, y = Interactions, fill = Significance)) +
  geom_boxplot() +
  scale_y_log10() +  # Logarithmic scale
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dashed", color = "black") +  # Vertical separators
  labs(
    title = "Box Plot of Reported and Significant Interactions by Method (Log Scale)",
    x = "Method",
    y = "Total Interactions (Log Scale)",
    fill = "Significance Level"
  ) +
  theme_minimal() +
  theme(legend.position = "right")







# Define the results directory
results_dir <- "results"

# Get the list of folders, excluding those with "dup" or "up2" in the name
folders <- list.files(results_dir, full.names = TRUE)
folders <- folders[!grepl("dup|dup2", folders)]

# Initialize an empty list
results_list <- list()

# Loop through each folder and read the RDS file
for (folder in folders) {
  file_path <- file.path(folder, "data", "decipher_scores_by_cluster.rds")
  
  if (file.exists(file_path)) {
    folder_name <- basename(folder)  # Get the folder name
    results_list[[folder_name]] <- readRDS(file_path)
  }
}

# Print summary of the stored results
print(names(results_list))


# Initialize a new list to store combined data frames per dataset
summarized_results <- list()

# Loop through each dataset
for (dataset in names(results_list)) {
  dataset_list <- results_list[[dataset]]
  
  # Convert list of data frames into a single data frame, adding a cell_type column
  combined_df <- do.call(rbind, lapply(names(dataset_list), function(cell_type) {
    df <- dataset_list[[cell_type]]
    df$cell_type <- cell_type  # Add cell type column
    return(df)
  }))
  
  # Store in the summarized list
  summarized_results[[dataset]] <- combined_df
}

# Print summary
print(names(summarized_results))

# Initialize a list to store summary statistics
summary_stats <- list()

# Define standard deviation thresholds from 0.1 to 4.0 in increments of 0.1
std_thresholds <- seq(0.1, 4.0, by = 0.1)

# Process each dataset
for (dataset in names(summarized_results)) {
  df <- summarized_results[[dataset]]
  
  # Compute the standard deviation of decipher_score
  std_dev <- sd(df$decipher_score, na.rm = TRUE)
  
  # Count total interactions (rows)
  total_interactions <- nrow(df)
  
  # Compute number of interactions above absolute thresholds
  threshold_counts <- sapply(std_thresholds, function(thresh) {
    sum(abs(df$decipher_score) > (thresh * std_dev), na.rm = TRUE)
  })
  
  # Store results in a named list
  summary_stats[[dataset]] <- c("Total Interactions" = total_interactions, setNames(threshold_counts, paste0("> ", std_thresholds, "σ")))
}

# Convert to a data frame for better visualization
summary_df <- do.call(rbind, summary_stats)

# Print the summary table
print(summary_df)

library(ggplot2)
library(tidyr)

# Convert summary_df to long format for ggplot
summary_long <- as.data.frame(summary_df)
summary_long$Dataset <- rownames(summary_long)
summary_long <- gather(summary_long, key = "Threshold", value = "Count", -Dataset)

# Ensure thresholds are ordered numerically
summary_long$Threshold <- as.numeric(gsub("> ", "", gsub("σ", "", summary_long$Threshold)))

# Create the boxplot
p <- ggplot(summary_long, aes(x = factor(Threshold), y = Count)) +
  geom_boxplot(aes(fill = factor(Threshold)), outlier.shape = NA) +  # Boxplot with no outliers
  scale_y_log10() +  # Log scale for better visualization
  labs(x = "Threshold (σ)", y = "Number of Interactions", title = "Distribution of Interactions Across Thresholds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Save as PNG
ggsave("interaction_thresholds_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)



library(ggplot2)
library(tidyr)

# Convert summary_df to long format for ggplot
summary_long <- as.data.frame(summary_df)
summary_long$Dataset <- rownames(summary_long)
summary_long <- gather(summary_long, key = "Threshold", value = "Count", -Dataset)

# Ensure thresholds are ordered numerically
summary_long$Threshold <- as.numeric(gsub("> ", "", gsub("σ", "", summary_long$Threshold)))

# Plot
p <- ggplot(summary_long, aes(x = Threshold, y = Count, color = Dataset)) +
  geom_line(size = 1) +
  scale_y_log10() +  # Log scale for better visualization
  labs(x = "Threshold (σ)", y = "Number of Interactions", title = "Interactions at Different Thresholds") +
  theme_minimal()

ggsave("interaction_thresholds.png", plot = p, width = 8, height = 6, dpi = 300)


#connectome

# Define the results directory
results_dir <- "results"

# Get the list of folders, excluding those with "dup" or "up2" in the name
folders <- list.files(results_dir, full.names = TRUE)
folders <- folders[!grepl("dup|dup2", folders)]

# Initialize an empty list
results_list <- list()

# Loop through each folder and read the RDS file
for (folder in folders) {
  file_path <- file.path(folder,"connectome/data", "connectome_results.rds")
  
  if (file.exists(file_path)) {
    folder_name <- basename(folder)  # Get the folder name
    results_list[[folder_name]] <- readRDS(file_path)
  }
}

# Print summary of the stored results
print(names(results_list))

#functions
preProcessConnectome <- function(connectome_df){
  result <- connectome_df%>%
    #filter
    filter(score!=Inf) %>%
    #rename
    rename(sender = source,
           receiver = target,
           prioritization_score = score) %>%
    #add interaction column
    mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
    #select columns of interest for downstream analysis
    select(sender, receiver, interaction, prioritization_score) %>%
    #organized based on score
    arrange(receiver, desc(prioritization_score))

  #scale prioritization score to make it comparable
  result <- scale_prioritization_score(result,"prioritization_score")

  return(result)
}

scale_prioritization_score <- function(df, score_column) {
  # Check if any value in the specified score column is negative
  has_negatives <- any(df[[score_column]] < 0)

  # Determine the scaling method based on the presence of negative values
  scaled_score <- if (!has_negatives) {
    # Scale between 0 and 1 for non-negative values
    (df[[score_column]] - min(df[[score_column]])) /
      (max(df[[score_column]]) - min(df[[score_column]]))
  } else {
    # Scale between -1 and 1 for ranges that include negative values
    ifelse(
      df[[score_column]] >= 0,
      df[[score_column]] / max(df[[score_column]]),
      df[[score_column]] / abs(min(df[[score_column]]))
    )
  }

  # Add the scaled score to the data frame
  df <- df %>% mutate(scaled_score = scaled_score)

  return(df)
}


#replace infinite values in the score column
results_list <- lapply(results_list, FUN = "preProcessConnectome")


summarized_results <- results_list
# Print summary
print(names(summarized_results))

# Initialize a list to store summary statistics
summary_stats <- list()

# Define standard deviation thresholds from 0.1 to 4.0 in increments of 0.1
std_thresholds <- seq(0.1, 4.0, by = 0.1)

# Process each dataset
for (dataset in names(summarized_results)) {
  df <- summarized_results[[dataset]]
  
  # Compute the standard deviation of decipher_score
  std_dev <- sd(df$prioritization_score, na.rm = TRUE)
  
  # Count total interactions (rows)
  total_interactions <- nrow(df)
  
  # Compute number of interactions above absolute thresholds
  threshold_counts <- sapply(std_thresholds, function(thresh) {
    sum(abs(df$prioritization_score) > (thresh * std_dev), na.rm = TRUE)
  })
  
  # Store results in a named list
  summary_stats[[dataset]] <- c("Total Interactions" = total_interactions, setNames(threshold_counts, paste0("> ", std_thresholds, "σ")))
}

# Convert to a data frame for better visualization
summary_df <- do.call(rbind, summary_stats)

# Print the summary table
print(summary_df)


library(ggplot2)
library(tidyr)

# Convert summary_df to long format for ggplot
summary_long <- as.data.frame(summary_df)
summary_long$Dataset <- rownames(summary_long)
summary_long <- gather(summary_long, key = "Threshold", value = "Count", -Dataset)

# Ensure thresholds are ordered numerically
summary_long$Threshold <- as.numeric(gsub("> ", "", gsub("σ", "", summary_long$Threshold)))

# Create the boxplot
p <- ggplot(summary_long, aes(x = factor(Threshold), y = Count)) +
  geom_boxplot(aes(fill = factor(Threshold)), outlier.shape = NA) +  # Boxplot with no outliers
  scale_y_log10() +  # Log scale for better visualization
  labs(x = "Threshold (σ)", y = "Number of Interactions", title = "Distribution of Interactions Across Thresholds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Save as PNG
ggsave("connectome_interaction_thresholds_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)

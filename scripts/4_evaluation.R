# Visualise the results of 100 runs of Decipher without setting a seed

# Load libraries
library(dplyr)
library(purrr)
library(ggplot2)

# Initialize empty lists to store the results
all_runs_dat_1 <- list()
all_runs_dat_2 <- list()

N <- 4
# Loop through all 100 runs
for (i in 1:N) {
  # Load the data
  dat_1 <- readRDS(paste0("sample_analysis/validity/data/Run_", i, "_decipher_scores_by_cluster.rds"))
  dat_2 <- readRDS(paste0("sample_analysis/validity/data/Run_", i, "_decipher_scores_by_regulon_and_cluster.rds"))

  # Add the run number and cell type columns, and combine the results
  combined_tibble_1 <- bind_rows(
    lapply(names(dat_1), function(cell_type) {
      dat_1[[cell_type]] %>%
        mutate(Cells = cell_type,
               run_number = i)
    })
  )

  combined_tibble_2 <- bind_rows(
    lapply(names(dat_2), function(cell_type) {
      dat_2[[cell_type]] %>%
        mutate(Cells = cell_type,
               run_number = i)
    })
  )

  # Store the results in the lists
  all_runs_dat_1[[i]] <- combined_tibble_1
  all_runs_dat_2[[i]] <- combined_tibble_2
}

# Combine all runs into single tibbles
final_combined_tibble_1 <- bind_rows(all_runs_dat_1)
final_combined_tibble_2 <- bind_rows(all_runs_dat_2)

#===================================================
# Step 2: Identify Top 10 Interactions Across Runs
#===================================================

# Identify top 10 interactions for each run and cell type
top_10_interactions <- final_combined_tibble_1 %>%
  group_by(Cells, run_number) %>%
  arrange(desc(decipher_score)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Count how often each interaction appears in the top 10
interaction_consistency <- top_10_interactions %>%
  group_by(Cells, interaction) %>%
  summarise(count = n(), .groups = "drop")

# Calculate the proportion of runs each interaction was in the top 10
interaction_consistency <- interaction_consistency %>%
  mutate(proportion = count / 100)

#===================================================
# Step 3: Visualize the Consistency of Interactions
#===================================================
# Plot the proportion of runs where each interaction was in the top 10
ggplot(interaction_consistency, aes(x = interaction, y = proportion, fill = Cells)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Proportion of Runs Where Interactions are in the Top 10",
       x = "Interaction",
       y = "Proportion of Runs (Top 10)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#===================================================
# Step 4: Track Decipher Scores Over Runs
#===================================================
# Filter for interactions that are consistently in the top 10 across runs
consistent_interactions <- interaction_consistency %>%
  #filter(proportion == 1) %>%
  select(Cells, interaction)

# Join with the original data to track the decipher score changes
score_trends <- final_combined_tibble_1 %>%
  semi_join(consistent_interactions, by = c("Cells", "interaction")) %>%
  group_by(Cells, interaction, run_number) %>%
  summarise(mean_decipher_score = mean(decipher_score), .groups = "drop")

# Plot the trends
ggplot(score_trends, aes(x = run_number, y = mean_decipher_score, color = interaction, group = interaction)) +
  geom_line() +
  facet_wrap(~Cells, scales = "free_y") +
  theme_minimal() +
  labs(title = "Decipher Score Trends for Consistent Top 10 Interactions",
       x = "Run Number",
       y = "Mean Decipher Score")


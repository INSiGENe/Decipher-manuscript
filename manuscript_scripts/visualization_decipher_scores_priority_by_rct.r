library(ggplot2)
library(dplyr)
library(tidyr)

# Paths to your datasets
severe_path <- "results/SevCOVID_Azimuthl2/data"
mild_path   <- "results/MilCOVID_Azimuthl2/data"

cell_type <- "cDC2"
k_param <- 1
# Load data
severe_scores <- readRDS(file.path(severe_path, "decipher_scores_by_cluster.rds")) %>%
  bind_rows() %>%
  mutate(condition = "severe")

mild_scores <- readRDS(file.path(mild_path, "decipher_scores_by_cluster.rds")) %>%
  bind_rows() %>%
  mutate(condition = "mild")

# Combine
combined_scores <- bind_rows(severe_scores, mild_scores)

# Join with interaction metadata and filter for CD8 TEM
decipher_df <- combined_scores %>%
  filter(receiver_cluster == cell_type) %>%
  select(interaction, condition, decipher_score) %>%
  mutate(condition = factor(condition, levels = c("mild", "severe")),
         decipher_score = replace_na(decipher_score, 0)) %>%
  distinct()

# Apply signed log1p transform
decipher_df <- decipher_df %>%
  mutate(decipher_score_log = sign(decipher_score) * log1p(abs(decipher_score)))

# Select top N interactions for plotting
top_interactions <- decipher_df %>%
  group_by(interaction) %>%
  summarize(max_score = max(abs(decipher_score))) %>%
  arrange(desc(max_score)) %>%
  slice_head(n = 30) %>%
  pull(interaction)

heatmap_data <- decipher_df %>%
  filter(interaction %in% top_interactions)

# Plot
p <- ggplot(heatmap_data, aes(x = condition, y = interaction, fill = decipher_score_log)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Log1p Score") +
  theme_minimal(base_size = 14) +
  labs(title = paste("Decipher Scores for",cell_type), x = "Condition", y = "Interaction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_filename <- paste0("differential_decipher_signalling_heatmap_",cell_type,"_",k_param,".png")
plot_path <- file.path("figures",plot_filename)

ggsave(plot_path,p)
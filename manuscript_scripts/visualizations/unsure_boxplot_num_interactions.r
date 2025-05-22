
# ==== Box-plot number of interactions V April 2025 ====
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
  geom_jitter(width = 0.2, size = 2.8, alpha = 0.8) +  # increased point size by ~40%
  scale_fill_manual(values = method_colors) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Number of Interactions"
  ) +
  theme_bw(base_size = 12) +  # using theme_bw as a base
  theme(
    legend.position = "none",
    # Updated tick labels with ~20% size increase and bold formatting:
    axis.text.x = element_text(angle = 0, vjust = 1, size = 17, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    # Updated axis titles with ~50% size increase and bold formatting:
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    # Optionally, you can enforce a bold base for all text if desired:
    text = element_text(face = "bold")
  )


ggsave(file.path(figures_folder,"boxplot_number_of_interactions_by_method.png"), plot = p, width = 8, height = 4, dpi = 300)

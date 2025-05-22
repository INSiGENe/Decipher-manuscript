library(dplyr)
library(ggplot2)
install.packages("viridis")
library(tidyr)
library(viridis)

# Start with the filtered data
library(dplyr)
library(tidyr)

# 1. Take your full table
summary_df <- combined_intersection_df %>%
  # 2. explode "A & B & C" into separate rows
  separate_rows(Intersection_Name, sep = " & ") %>%
  # 3. rename the exploded column to "method"
  rename(method = Intersection_Name) %>%
  # 4. sum up counts by dataset, method, and degree
  group_by(Dataset, method, Degree) %>%
  summarise(total_Count = sum(Count), .groups = "drop")

summary_df$method <- factor(summary_df$method,levels=desired_method_order)

# Step 4: Plot
p <- ggplot(summary_df, aes(x = method, y = total_Count, fill = as.factor(Degree))) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6, size = 1.5) +
  scale_fill_viridis_d(name = NULL, labels = degree_labels, option = "D") +
  labs(y = "Overlap", x = NULL) +
  theme_minimal(base_size = 17) +
  scale_y_continuous(position = "right")+
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    axis.text.x = element_text(size = 17, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold")
  )

# Save plot
ggsave(filename = file.path(figures_folder, "overlap_plot_degree.png"),
       plot = p, width = 8, height = 4)

library(dplyr)
library(ggplot2)

# Create 'figures' folder if it doesn't exist
if (!dir.exists("figures")) dir.create("figures")

# Extract sample-level metadata
sample_metadata <- test@meta.data %>%
  select(sample_id, days_since_symptom_onset, severity_group, age) %>%
  distinct()

# 1. Number of samples by day since symptom onset and severity group
sample_counts_by_day <- sample_metadata %>%
  group_by(days_since_symptom_onset, severity_group) %>%
  summarise(n_samples = n(), .groups = "drop")

  # Replace NA with 0 for days_since_symptom_onset
sample_metadata <- test@meta.data %>%
  select(sample_id, days_since_symptom_onset, severity_group, age) %>%
  distinct() %>%
  mutate(days_since_symptom_onset = ifelse(is.na(days_since_symptom_onset), 0, days_since_symptom_onset))

# Group and count
sample_counts_by_day <- sample_metadata %>%
  group_by(days_since_symptom_onset, severity_group) %>%
  summarise(n_samples = n(), .groups = "drop")


# Plot: Sample count by day
plot1 <- ggplot(sample_counts_by_day, aes(x = days_since_symptom_onset, y = n_samples, fill = severity_group)) +
  geom_col(position = "dodge") +
  labs(title = "Samples by Day After Symptom Onset", x = "Days Since Onset", y = "Number of Samples") +
  theme_minimal()

ggsave("figures/samples_by_day.png", plot = plot1, width = 3, height = 5, dpi = 300)

# 2. Age distribution by severity group
plot2 <- ggplot(sample_metadata, aes(x = severity_group, y = age, fill = severity_group)) +
  geom_boxplot() +
  labs(title = "Age Distribution by Severity Group", x = "Severity Group", y = "Age") +
  theme_minimal()

ggsave("figures/age_by_severity_group.png", plot = plot2, width = 3, height = 5, dpi = 300)

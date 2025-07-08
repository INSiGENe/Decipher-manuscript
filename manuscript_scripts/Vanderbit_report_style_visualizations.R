# ---- Libraries ----
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(scales)
library(RColorBrewer)
library(tibble)
library(patchwork)

#  Functions ----
calculate_pct_change <- function(new_val, old_val) {
  # ifelse(condition, value_if_true, value_if_false)
  change <- ifelse(old_val == 0,                     # Test condition (vectorized)
                   ifelse(new_val > 0, Inf, 0),      # Value if old_val is 0 (handle 0 new_val too)
                   ((new_val - old_val) / old_val) * 100) # Value if old_val is not 0
  return(change)
}
calculate_pct_change <- function(new_val, old_val) {
  change <- ifelse(old_val == 0,
                   ifelse(new_val > 0, Inf, 0),
                   ((new_val - old_val) / old_val) * 100)
  return(change)
}

load_seurat <- function(path) {
  rds_path <- file.path(path, "pre_processing", "seurat_object_oi.rds")
  if (!file.exists(rds_path)) stop("Missing file:", rds_path)
  readRDS(rds_path)
}
de_by_cluster <- function(seurat_obj, case_group, control_group, logfc_threshold, up_col, down_col) {
  seurat_obj$cluster <- as.character(seurat_obj$predicted.celltype.l2)
  clusters <- unique(seurat_obj$cluster)
  results <- list()

  for (cluster in clusters) {
    cells <- WhichCells(seurat_obj, cells = rownames(seurat_obj@meta.data)[seurat_obj$cluster == cluster])
    if (length(cells) < 2) next

    subset_obj <- subset(seurat_obj, cells = cells)
    Idents(subset_obj) <- subset_obj$severity_group
    if (!all(c(case_group, control_group) %in% levels(Idents(subset_obj)))) next
    if (sum(Idents(subset_obj) == case_group) < 3 || sum(Idents(subset_obj) == control_group) < 3) next

    markers <- FindMarkers(
      subset_obj,
      ident.1 = case_group,
      ident.2 = control_group,
      test.use = "wilcox",
      logfc.threshold = logfc_threshold,
      min.pct = 0.1,
      random.seed = 123
    )

    if (!nrow(markers)) {
      results[[cluster]] <- tibble(Cluster = cluster, !!up_col := 0, !!down_col := 0)
    } else {
      results[[cluster]] <- tibble(
        Cluster = cluster,
        !!up_col := sum(markers$avg_log2FC > 0),
        !!down_col := sum(markers$avg_log2FC < 0)
      )
    }
  }

  bind_rows(results)
}

run_all_comparisons <- function(dataset_paths, comparisons, logfc_threshold) {
  result_list <- list()
  seurat_objects <- list()

  for (dataset_name in names(dataset_paths)) {
    message("Processing ", dataset_name)
    seurat_obj <- load_seurat(dataset_paths[[dataset_name]])
    seurat_objects[[dataset_name]] <- seurat_obj

    comp <- comparisons[[dataset_name]]
    df <- de_by_cluster(seurat_obj, comp$case, comp$control, logfc_threshold, comp$up_col, comp$down_col)
    result_list[[dataset_name]] <- df
  }

  list(results = result_list, seurat_objects = seurat_objects)
}
calculate_condition_normalized_proportions <- function(seurat_objects) {
  proportions_list <- list()

  for (dataset_name in names(seurat_objects)) {
    obj <- seurat_objects[[dataset_name]]
    meta <- obj@meta.data

    # Step 1: Count per cell type × severity_group
    counts <- meta %>%
      count(predicted.celltype.l2, severity_group, name = "RawCount")

    # Step 2: Normalize by total cells per severity_group
    total_by_condition <- meta %>%
      count(severity_group, name = "TotalCells")

    counts <- counts %>%
      left_join(total_by_condition, by = "severity_group") %>%
      mutate(NormalizedCount = RawCount / TotalCells)

    # Step 3: Normalize within each cell type to sum to 100%
    proportions <- counts %>%
      group_by(predicted.celltype.l2) %>%
      mutate(Proportion = NormalizedCount / sum(NormalizedCount) * 100) %>%
      ungroup() %>%
      mutate(Dataset = dataset_name)

    proportions_list[[dataset_name]] <- proportions
  }

  bind_rows(proportions_list)
}
combine_deg_counts <- function(df_severe, df_mild) {
  full_join(df_severe, df_mild, by = "Cluster") %>%
    replace_na(list(
      Severe_Up = 0, Severe_Down = 0,
      Mild_Up = 0, Mild_Down = 0
    )) %>%
    arrange(Cluster) %>%
    mutate(
      Delta = Severe_Up + Severe_Down - Mild_Up - Mild_Down,
      Total = Severe_Up + Severe_Down + Mild_Up + Mild_Down,
      Severe_pct = 100 * (Severe_Up + Severe_Down) / Total,
      Mild_pct = 100 * (Mild_Up + Mild_Down) / Total
    )
}



# ---- reproducible seed ----
set.seed(123)

# load data ----
combined <- readRDS("pre_processing_test/data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")

# umap plot ----
# Base UMAP plot without labels or legend
p <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "severity_group",
  group.by = "predicted.celltype.l2",
  pt.size = 0.5
) + theme(legend.position = "bottom")

# Save the plot
ggsave("Manuscript_jan_2025/figures/umap_split_by_condition.png", plot = p, width = 10, height = 6, dpi = 300)

# umap CD63 expression ----
p <- FeaturePlot(
  combined,
  features = "CD63",
  reduction = "umap",
  split.by = "condition",
  pt.size = 0.5
) + theme(legend.position = "right")
ggsave("figures/umap_CD63_expression.png", plot = p, width = 10, height = 6, dpi = 300)

# abundance V1 ----
metadata <- combined@meta.data %>%
  select(severity_group, predicted.celltype.l2)

# Calculate counts per group
cell_counts <- metadata %>%
  group_by(severity_group, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  # Ensure all severity levels are factors for pivoting
  mutate(severity_group = factor(severity_group, levels = c("Healthy", "Moderate", "Severe")))

# Pivot wider to easily compare conditions
counts_wide <- cell_counts %>%
  pivot_wider(names_from = severity_group, values_from = n_cells, values_fill = 0) # Fill missing with 0 count


# --- Corrected mutate block ---
abundance_changes <- counts_wide %>%
  mutate(
    pct_change_M_vs_H = calculate_pct_change(Moderate, Healthy),
    pct_change_S_vs_H = calculate_pct_change(Severe, Healthy),
    pct_change_S_vs_M = calculate_pct_change(Severe, Moderate),

    # Determine the label text based on the logic - Use & instead of &&
    label_Moderate = case_when(
      Healthy == 0 & Moderate > 0 & Severe == 0 ~ "em (M)", # Emergent only in Moderate (Use &)
      Healthy > 0 ~ paste0(ifelse(pct_change_M_vs_H >= 0, "\u2191", "\u2193"),
                           # Handle Inf case explicitly if needed, otherwise it might show "Inf%"
                           ifelse(is.infinite(pct_change_M_vs_H), "Inf", round(abs(pct_change_M_vs_H))),
                           "%"),
      TRUE ~ "" # No label otherwise
    ),
    label_Severe = case_when(
       Healthy == 0 & Moderate == 0 & Severe > 0 ~ "em (S)", # Emergent only in Severe (Use &)
       Healthy == 0 & Moderate > 0 & Severe > 0 ~ paste0("SvM ", # Compare Severe vs Moderate (Use &)
                                                           ifelse(pct_change_S_vs_M >= 0, "\u2191", "\u2193"),
                                                           # Handle potential Inf for SvM
                                                           ifelse(is.infinite(pct_change_S_vs_M), "Inf", round(abs(pct_change_S_vs_M))),
                                                           "%"),
       Healthy > 0 ~ paste0(ifelse(pct_change_S_vs_H >= 0, "\u2191", "\u2193"),
                            # Handle potential Inf for S vs H
                            ifelse(is.infinite(pct_change_S_vs_H), "Inf", round(abs(pct_change_S_vs_H))),
                            "%"),
       TRUE ~ "" # No label otherwise
    )
  ) %>%
  # Select relevant columns and pivot longer for easy joining
  select(predicted.celltype.l2, label_Moderate, label_Severe) %>%
  pivot_longer(
      cols = starts_with("label_"),
      names_to = "severity_group_label", # Temp column name
      values_to = "abundance_label"
  ) %>%
  # Map back to original severity group names for faceting
  mutate(severity_group = case_when(
      severity_group_label == "label_Moderate" ~ "Moderate",
      severity_group_label == "label_Severe" ~ "Severe",
      TRUE ~ NA_character_
  )) %>%
  # Filter out rows where no label needed or mapping failed
  filter(!is.na(severity_group), abundance_label != "", !is.na(abundance_label)) %>%
  # Select final columns needed for annotation data frame
  select(predicted.celltype.l2, severity_group, abundance_label)


# --- Step 3 (Revised): Determine Label Positions ---

# Get UMAP coordinates along with metadata
umap_coords <- FetchData(combined, vars = c("UMAP_1", "UMAP_2", "severity_group", "predicted.celltype.l2"))

# Calculate median coordinates for each cell type within the 'Healthy' group ONLY
healthy_label_positions <- umap_coords %>%
  filter(severity_group == "Healthy") %>%
  group_by(predicted.celltype.l2) %>%
  summarise(
    healthy_median_umap1 = median(UMAP_1),
    healthy_median_umap2 = median(UMAP_2),
    .groups = 'drop'
  )

# Calculate median coordinates within EACH group (as a fallback for missing cell types in Healthy)
all_label_positions <- umap_coords %>%
  group_by(severity_group, predicted.celltype.l2) %>%
  summarise(
    group_median_umap1 = median(UMAP_1),
    group_median_umap2 = median(UMAP_2),
    .groups = 'drop'
  )

# --- Step 4 (Revised): Create Annotation Data ---

# Start with the labels calculated earlier
annotation_data <- abundance_changes %>%
  # Join with the Healthy positions (based on cell type only)
  left_join(healthy_label_positions, by = "predicted.celltype.l2") %>%
  # Join with the group-specific positions (based on cell type and the group the label is for)
  left_join(all_label_positions, by = c("predicted.celltype.l2", "severity_group")) %>%
  # Determine final coordinates: Use Healthy position if available, otherwise use group-specific position
  mutate(
    final_umap1 = coalesce(healthy_median_umap1, group_median_umap1),
    final_umap2 = coalesce(healthy_median_umap2, group_median_umap2)
  ) %>%
  # Keep only necessary columns and filter out any rows where coordinates couldn't be determined
  # (Shouldn't happen if cell types exist in at least one group, but good practice)
  select(predicted.celltype.l2, severity_group, abundance_label, final_umap1, final_umap2) %>%
  filter(!is.na(final_umap1) & !is.na(final_umap2))

# --- Step 5: Generate Base Plot and Add Annotations ---
# (Use the annotation_data frame as created above)

# Generate the base DimPlot
p <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "severity_group",
  group.by = "predicted.celltype.l2",
  pt.size = 0.5
) + theme(legend.position = "bottom")

# Add the annotation layer using ggrepel, now using final_umap1/2
p_annotated <- p +
  ggrepel::geom_text_repel(
    data = annotation_data,
    # Use the final calculated coordinates
    aes(x = final_umap1, y = final_umap2, label = abundance_label),
    color = "black", # Adjust color as needed
    size = 3,       # Adjust size as needed
    fontface = "bold",
    bg.color = "white", # Add background color for readability
    bg.r = 0.15,        # Background radius
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf # Or set a limit
  )

# Print the plot
ggsave("Manuscript_jan_2025/figures/umap_split_by_condition_ct_abundance.png", plot = p_annotated, width = 10, height = 6, dpi = 300)

#V3


# --- Step 0: Ensure Seurat object 'combined' exists ---
# combined <- ... your Seurat object ...

# --- Step 1: Calculate Cell Counts ---
metadata <- combined@meta.data %>%
  select(severity_group, predicted.celltype.l2)

cell_counts <- metadata %>%
  group_by(severity_group, predicted.celltype.l2) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  mutate(severity_group = factor(severity_group, levels = c("Healthy", "Moderate", "Severe")))

counts_wide <- cell_counts %>%
  pivot_wider(names_from = severity_group, values_from = n_cells, values_fill = 0)

# --- Step 2 (Modified): Calculate Changes, Check Criteria, Create Labels ---



# Calculate changes, determine labels AND highlighting status
abundance_calcs_and_labels <- counts_wide %>%
  mutate(
    # Calculate changes
    pct_change_M_vs_H = calculate_pct_change(Moderate, Healthy),
    pct_change_S_vs_H = calculate_pct_change(Severe, Healthy),
    pct_change_S_vs_M = calculate_pct_change(Severe, Moderate),

    # Check minimum count criterion
    meets_min_count = Healthy >= 500 | Moderate >= 500 | Severe >= 500,

    # --- Determine Labels and Highlight Status for Moderate ---
    label_Moderate = case_when(
      Healthy == 0 & Moderate > 0 & Severe == 0 ~ "em (M)",
      Healthy > 0 ~ paste0(ifelse(pct_change_M_vs_H >= 0, "\u2191", "\u2193"),
                           ifelse(is.infinite(pct_change_M_vs_H), "Inf", round(abs(pct_change_M_vs_H))),
                           "%"),
      TRUE ~ ""
    ),
    # Check highlight condition for Moderate label
    highlight_Moderate = meets_min_count &
                         label_Moderate != "" & # Only highlight if there is a label
                         (pct_change_M_vs_H >= 100 | pct_change_M_vs_H <= -50),

    # --- Determine Labels and Highlight Status for Severe ---
    # First, determine which comparison is relevant for the Severe label
    severe_comparison_type = case_when(
        Healthy == 0 & Moderate == 0 & Severe > 0 ~ "EmergentS",
        Healthy == 0 & Moderate > 0 & Severe > 0 ~ "SvM",
        Healthy > 0 ~ "SvH",
        TRUE ~ "None"
    ),
    # Create Severe label based on comparison type
    label_Severe = case_when(
      severe_comparison_type == "EmergentS" ~ "em (S)",
      severe_comparison_type == "SvM" ~ paste0("SvM ",
                                                ifelse(pct_change_S_vs_M >= 0, "\u2191", "\u2193"),
                                                ifelse(is.infinite(pct_change_S_vs_M), "Inf", round(abs(pct_change_S_vs_M))),
                                                "%"),
      severe_comparison_type == "SvH" ~ paste0(ifelse(pct_change_S_vs_H >= 0, "\u2191", "\u2193"),
                                                ifelse(is.infinite(pct_change_S_vs_H), "Inf", round(abs(pct_change_S_vs_H))),
                                                "%"),
      TRUE ~ ""
    ),
    # Check highlight condition for Severe label based on the relevant comparison
    highlight_Severe = meets_min_count &
                       label_Severe != "" & # Only highlight if there is a label
                       case_when(
                         # Emergent S automatically meets criteria if min count is met (Inf increase)
                         severe_comparison_type == "EmergentS" ~ TRUE,
                         # Check SvM change
                         severe_comparison_type == "SvM" ~ (pct_change_S_vs_M >= 100 | pct_change_S_vs_M <= -50),
                         # Check SvH change
                         severe_comparison_type == "SvH" ~ (pct_change_S_vs_H >= 100 | pct_change_S_vs_H <= -50),
                         # Don't highlight otherwise
                         TRUE ~ FALSE
                       )
  )

# --- Step 2b: Pivot longer and combine labels/highlights ---

# Select and pivot labels
labels_long <- abundance_calcs_and_labels %>%
  select(predicted.celltype.l2, label_Moderate, label_Severe) %>%
  pivot_longer(
    cols = starts_with("label_"),
    names_to = "severity_group_label",
    values_to = "abundance_label",
    names_prefix = "label_"
  )

# Select and pivot highlights
highlights_long <- abundance_calcs_and_labels %>%
  select(predicted.celltype.l2, highlight_Moderate, highlight_Severe) %>%
  pivot_longer(
    cols = starts_with("highlight_"),
    names_to = "severity_group_label",
    values_to = "highlight_status",
    names_prefix = "highlight_"
  )

# Combine labels and highlights
abundance_changes_long <- labels_long %>%
  left_join(highlights_long, by = c("predicted.celltype.l2", "severity_group_label")) %>%
  # Map back to original severity group names
  mutate(severity_group = severity_group_label) %>% # Simplified mapping
  # Assign color based on highlight status
  mutate(
    highlight_color = ifelse(highlight_status, "red", "black") # Set colors here
  ) %>%
  # Filter out rows with no labels or failed mapping
  filter(abundance_label != "", !is.na(abundance_label), !is.na(severity_group)) %>%
  select(predicted.celltype.l2, severity_group, abundance_label, highlight_color)

# --- Step 3: Determine Label Positions (using Healthy baseline + fallback) ---

umap_coords <- FetchData(combined, vars = c("UMAP_1", "UMAP_2", "severity_group", "predicted.celltype.l2"))

healthy_label_positions <- umap_coords %>%
  filter(severity_group == "Healthy") %>%
  group_by(predicted.celltype.l2) %>%
  summarise(
    healthy_median_umap1 = median(UMAP_1),
    healthy_median_umap2 = median(UMAP_2),
    .groups = 'drop'
  )

all_label_positions <- umap_coords %>%
  group_by(severity_group, predicted.celltype.l2) %>%
  summarise(
    group_median_umap1 = median(UMAP_1),
    group_median_umap2 = median(UMAP_2),
    .groups = 'drop'
  )

# --- Step 4: Create Final Annotation Data ---

annotation_data <- abundance_changes_long %>%
  left_join(healthy_label_positions, by = "predicted.celltype.l2") %>%
  left_join(all_label_positions, by = c("predicted.celltype.l2", "severity_group")) %>%
  mutate(
    final_umap1 = coalesce(healthy_median_umap1, group_median_umap1),
    final_umap2 = coalesce(healthy_median_umap2, group_median_umap2)
  ) %>%
  select(predicted.celltype.l2, severity_group, abundance_label, highlight_color, final_umap1, final_umap2) %>%
  filter(!is.na(final_umap1) & !is.na(final_umap2))

# --- Step 5: Generate Base Plot and Add Annotations with Color ---

# --- Add this section before your DimPlot call ---

# 1. Get all unique cell type names from your Seurat object
#cell_types <- levels(combined$predicted.celltype.l2) # If it's a factor
# or if it's not a factor:
cell_types <- unique(combined$predicted.celltype.l2)
# print(cell_types) # Good to check if "Eryth" is indeed present

# 2. Create a named vector mapping each cell type to a valid color
# Option A: Using scales::hue_pal for distinct colors

n_cell_types <- length(cell_types)
my_cell_type_colors <- setNames(hue_pal()(n_cell_types), cell_types)
p <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "severity_group",
  group.by = "predicted.celltype.l2",
  pt.size = 0.5,
  cols = my_cell_type_colors # <--- Add this argument
) + theme(legend.position = "bottom")

p_annotated <- p +
  ggrepel::geom_text_repel(
    data = annotation_data,
    aes(x = final_umap1,
        y = final_umap2,
        label = abundance_label,
        color = highlight_color), # Map color aesthetic
    size = 3,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf,
    seed = 42 # for reproducibility of label placement
  ) +
  scale_color_manual(values = c("black" = "black",
                                "red" = "red")) # Add other colors if necessary

# Print the plot
ggsave("Manuscript_jan_2025/figures/umap_split_by_condition_ct_abundance_highlight.png", plot = p_annotated, width = 10, height = 6, dpi = 300)


#V4

# --- Step 1: Determine Overall Highlight Status per Cell Type ---

print("Determining overall highlight status for each cell type...")
highlight_status_overall <- abundance_calcs_and_labels %>%
  # Determine if highlight criteria met in EITHER Moderate OR Severe comparison
  mutate(
    highlight_overall = meets_min_count & (highlight_Moderate | highlight_Severe)
    ) %>%
  # Select distinct status per cell type
  select(predicted.celltype.l2, highlight_overall) %>%
  distinct() # Ensure one row per cell type

# Check how many types are highlighted
# print(table(highlight_status_overall$highlight_overall))

# --- Step 2: Create Modified Color Palette ---

print("Creating modified color palette (grey for non-highlighted)...")
# Start with the original unique colors
modified_colors <- my_cell_type_colors
# Define the grey color for non-highlighted types
grey_color <- "grey85" # Or "grey80", "grey90" etc.

# Identify cell types that should NOT be highlighted
types_to_grey_out <- highlight_status_overall %>%
  filter(!highlight_overall) %>%
  pull(predicted.celltype.l2)

# Set the color to grey for these types in the modified palette
if (length(types_to_grey_out) > 0) {
  print(paste("Setting color to", grey_color, "for:", paste(types_to_grey_out, collapse=", ")))
  modified_colors[types_to_grey_out] <- grey_color
} else {
  print("No cell types met the criteria to be greyed out.")
}

# Verify the modified palette (optional)
# print(modified_colors)

# --- Step 3: Create DimPlot with the Modified Palette ---

print("Generating DimPlot with modified colors...")
p_modified_colors <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "severity_group",
  group.by = "predicted.celltype.l2",
  cols = modified_colors, # <-- Use the palette with greyed-out types
  pt.size = 0.5
) +
theme(legend.position = "bottom") +
ggtitle("UMAP by Cell Type (Highlighted if Abundance Change Criteria Met)")

# --- Step 4: Add Text Labels (Optional - e.g., All Black) ---

# Assuming 'annotation_data' exists with label positions and text:
# predicted.celltype.l2, severity_group, abundance_label, final_umap1, final_umap2

print("Adding text labels (black)...")
# If you want red/black labels, keep highlight_color in annotation_data
# and use aes(color=highlight_color) + scale_color_manual(...) below.
# For all black labels:
p_final_annotated <- p_modified_colors +
  ggrepel::geom_text_repel(
    data = annotation_data,
    aes(x = final_umap1,
        y = final_umap2,
        label = abundance_label), # Label text
    color = "black",             # Set all labels to black
    size = 3,
    fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.5, point.padding = 0.5,
    max.overlaps = Inf, seed = 42
  )

# --- Step 5: Save the Final Plot ---
print("Saving final plot...")
ggsave("Manuscript_jan_2025/figures/umap_split_condition_clusters_highlighted_overall.png", # Updated filename maybe
       plot = p_final_annotated, width = 10, height = 6, dpi = 300)
print("Plot saved.")


#v5
# Load necessary libraries


# --- Prerequisite Steps (Assumed to be done already) ---
# 1. 'combined' Seurat object exists
# 2. 'abundance_calcs_and_labels' exists (contains highlight_Moderate, etc.)
# 3. 'cell_types' exists (vector of all unique cell type names)
#       cell_types <- unique(combined$predicted.celltype.l2)

# --- Step 1: Determine Overall Highlight Status per Cell Type ---

print("Determining overall highlight status for each cell type...")
highlight_status_overall <- abundance_calcs_and_labels %>%
  mutate(
    highlight_overall = meets_min_count & (highlight_Moderate | highlight_Severe)
  ) %>%
  select(predicted.celltype.l2, highlight_overall) %>%
  distinct()

# --- Step 2: Identify Highlighted Types and Generate Distinct Palette ---

# Get the names of cell types to be highlighted
highlighted_types <- highlight_status_overall %>%
  filter(highlight_overall) %>%
  pull(predicted.celltype.l2)

n_highlighted <- length(highlighted_types)
print(paste("Number of cell types to highlight with distinct colors:", n_highlighted))

# Generate distinct colors ONLY for the highlighted types
distinct_colors <- character(0) # Initialize empty vector
if (n_highlighted > 0) {
    print("Generating distinct color palette...")
    # Choose a suitable RColorBrewer palette (e.g., Set1, Set3, Paired)
    # Handle cases where n_highlighted is small or larger than palette limits
    palette_choice <- "Set1" # Good distinct colors (max 9)
    max_colors_brewer <- brewer.pal.info[palette_choice, "maxcolors"]

    if (n_highlighted <= max_colors_brewer) {
      # Use brewer.pal, ensuring n is at least 3 for some palettes if necessary
      distinct_colors <- brewer.pal(max(3, n_highlighted), palette_choice)[1:n_highlighted]
    } else {
      # Fallback if more colors are needed than the chosen palette offers
      palette_choice <- "Set3" # Try another one (max 12)
      max_colors_brewer <- brewer.pal.info[palette_choice, "maxcolors"]
      if (n_highlighted <= max_colors_brewer) {
         distinct_colors <- brewer.pal(n_highlighted, palette_choice)
      } else {
        # Absolute fallback if still too many categories
        print(paste("Warning: Number of highlighted types (", n_highlighted, ") exceeds tested RColorBrewer limits. Using scales::hue_pal as fallback."))
        distinct_colors <- scales::hue_pal()(n_highlighted)
        # distinct_colors <- viridisLite::viridis(n_highlighted)
      }
    }
} else {
    print("No types met the criteria for highlighting.")
}

# --- Step 3: Create Final Merged Color Palette ---

print("Creating final merged color palette...")
grey_color <- "grey85"

# Start with all types assigned grey
final_color_palette <- setNames(rep(grey_color, length(cell_types)), cell_types)

# If there are highlighted types, map them to the distinct colors
if (n_highlighted > 0) {
  highlighted_color_map <- setNames(distinct_colors, highlighted_types)
  # Update the main palette: colors for highlighted types will be replaced
  final_color_palette[names(highlighted_color_map)] <- highlighted_color_map
}

# Verify the final palette (optional)
# print(final_color_palette)
# print(paste("Colors assigned to highlighted types:", paste(final_color_palette[highlighted_types], collapse=", ")))
# print(paste("Colors assigned to non-highlighted types:", paste(unique(final_color_palette[!names(final_color_palette) %in% highlighted_types]), collapse=", ")))


# --- Step 4: Create DimPlot with the Final Palette ---

print("Generating DimPlot with final colors...")
p_final_colors <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "severity_group",
  group.by = "predicted.celltype.l2",
  cols = final_color_palette, # <-- Use the final palette
  pt.size = 0.5
) +
theme(legend.position = "bottom") +
ggtitle("UMAP by Cell Type (Highlighted Clusters with Distinct Colors)")


# --- Step 5: Add Text Labels (Optional - e.g., All Black) ---

# Assuming 'annotation_data' exists with label positions and text

print("Adding text labels (black)...")
p_final_annotated <- p_final_colors +
  ggrepel::geom_text_repel(
    data = annotation_data,
    aes(x = final_umap1,
        y = final_umap2,
        label = abundance_label), # Label text
    color = "black",             # Set all labels to black
    size = 3,
    fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.5, point.padding = 0.5,
    max.overlaps = Inf, seed = 42
  )

# --- Step 6: Save the Final Plot ---
print("Saving final plot...")
ggsave("Manuscript_jan_2025/figures/umap_split_condition_clusters_highlighted_distinct_colors.png", # Updated filename
       plot = p_final_annotated, width = 10, height = 6, dpi = 300)
print("Plot saved.")

#V6


# --- Prerequisite Steps (Assumed to be done already) ---
# 1. 'combined' Seurat object exists
# 2. 'abundance_calcs_and_labels' exists
# 3. 'cell_types' exists
# 4. 'my_cell_type_colors' may exist (though we generate a new one)
# 5. 'annotation_data' exists: This data frame should contain the text labels
#    and positions calculated earlier. It needs columns:
#    predicted.celltype.l2, severity_group, abundance_label, final_umap1, final_umap2

# --- Step 1: Determine Overall Highlight Status per Cell Type ---
# (Same as before)
print("Determining overall highlight status for each cell type...")
highlight_status_overall <- abundance_calcs_and_labels %>%
  mutate(
    highlight_overall = meets_min_count & (highlight_Moderate | highlight_Severe)
  ) %>%
  select(predicted.celltype.l2, highlight_overall) %>%
  distinct()

# --- Step 2: Identify Highlighted Types and Generate Distinct Palette ---
# (Same as before)
print("Identifying highlighted types and generating palette...")
highlighted_types <- highlight_status_overall %>%
  filter(highlight_overall) %>%
  pull(predicted.celltype.l2)
n_highlighted <- length(highlighted_types)
distinct_colors <- character(0)
# ... (palette generation logic as before) ...
if (n_highlighted > 0) {
    palette_choice <- "Set1"
    max_colors_brewer <- brewer.pal.info[palette_choice, "maxcolors"]
    if (n_highlighted <= max_colors_brewer) {
      distinct_colors <- brewer.pal(max(3, n_highlighted), palette_choice)[1:n_highlighted]
    } else {
      palette_choice <- "Set3"
      max_colors_brewer <- brewer.pal.info[palette_choice, "maxcolors"]
      if (n_highlighted <= max_colors_brewer) {
         distinct_colors <- brewer.pal(n_highlighted, palette_choice)
      } else {
        print(paste("Warning: Using scales::hue_pal fallback for", n_highlighted, "colors."))
        distinct_colors <- scales::hue_pal()(n_highlighted)
      }
    }
}

    # --- Step 3: Create Final Merged Color Palette ---
    # (Same as before)
    print("Creating final merged color palette...")
    grey_color <- "grey85"
    final_color_palette <- setNames(rep(grey_color, length(cell_types)), cell_types)
    if (n_highlighted > 0) {
    highlighted_color_map <- setNames(distinct_colors, highlighted_types)
    final_color_palette[names(highlighted_color_map)] <- highlighted_color_map
    }

    # --- Step 4: Create DimPlot with the Final Palette ---
    # (Same as before)
    print("Generating DimPlot with final colors...")
    p_final_colors <- DimPlot(
    combined,
    reduction = "umap",
    split.by = "severity_group",
    group.by = "predicted.celltype.l2",
    cols = final_color_palette,
    pt.size = 0.5
    ) +
    theme(legend.position = "bottom") +
    ggtitle("UMAP by Cell Type (Highlighted Clusters with Distinct Colors)")

# --- Step 5 (MODIFIED): Filter Annotation Data and Add Labels ---

print("Preparing annotation data for highlighted clusters only...")
# We need the 'annotation_data' from previous steps
# Ensure it has: predicted.celltype.l2, severity_group, abundance_label, final_umap1, final_umap2

# Join the overall highlight status to the annotation data
annotation_data_with_highlight <- annotation_data %>%
  left_join(highlight_status_overall, by = "predicted.celltype.l2")

# Filter to keep labels ONLY for highlighted cell types
filtered_annotation_data <- annotation_data_with_highlight %>%
  filter(highlight_overall == TRUE)

print(paste("Number of labels to plot after filtering:", nrow(filtered_annotation_data)))

# Add the filtered text labels to the plot
print("Adding text labels for highlighted clusters only (black)...")
p_final_annotated <- p_final_colors +
  ggrepel::geom_text_repel(
    # Use the FILTERED data frame
    data = filtered_annotation_data,
    aes(x = final_umap1,
        y = final_umap2,
        label = abundance_label), # Label text
    color = "black",             # Set labels to black
    size = 3,
    fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.5, point.padding = 0.5,
    max.overlaps = Inf, seed = 42
  )

# --- Step 6: Save the Final Plot ---
# (Same as before)
print("Saving final plot...")
ggsave("Manuscript_jan_2025/figures/umap_split_condition_clusters_highlighted_labels_filtered.png", # Updated filename
       plot = p_final_annotated, width = 10, height = 6, dpi = 300)
print("Plot saved.")

# Display plot (optional)
# print(p_final_annotated)
# Display plot (optional)
# print(p_final_annotated)

#V7
# V8 - Using Combined scCODA results for highlighting and labeling


# Load scCODA results
files <- c(
  moderate = "pre_processing_test/data/SevMilCOVID/results_sccoda_moderate_vs_healthy.csv",
  severe   = "pre_processing_test/data/SevMilCOVID/results_sccoda_severe_vs_healthy.csv"
)

sccoda_results <- bind_rows(
  lapply(names(files), function(name) {
    df <- read.csv(files[[name]])
    df$comparison_source <- name
    df
  })
) %>%
  rename_with(~gsub("[^[:alnum:]_]", "_", .)) %>%
  rename(CellType = Cell_Type) %>%
  mutate(
    Final_Parameter = as.numeric(Final_Parameter),
    log2_fold_change = as.numeric(log2_fold_change)
  )

# Get significant cell types
sig_covariates <- c("severity_group[T.Moderate]", "severity_group[T.Severe]")
significant_types <- sccoda_results %>%
  filter(Covariate %in% sig_covariates, Final_Parameter != 0) %>%
  pull(CellType) %>%
  unique()

highlight_status <- data.frame(predicted.celltype.l2 = cell_types) %>%
  mutate(highlight = predicted.celltype.l2 %in% significant_types)

# Generate color palette
highlighted <- highlight_status %>% filter(highlight) %>% pull(predicted.celltype.l2)
n <- length(highlighted)
palette_name <- "Set1"
colors <- if (n <= brewer.pal.info[palette_name, "maxcolors"]) {
  brewer.pal(max(3, n), palette_name)[1:n]
} else {
  print("Too many highlights for Set1, using hue palette")
  hue_pal()(n)
}
final_colors <- setNames(rep("grey85", length(cell_types)), cell_types)
final_colors[highlighted] <- colors

# Plot
p <- DimPlot(combined, reduction = "umap", split.by = "severity_group",
             group.by = "predicted.celltype.l2", cols = final_colors,
             pt.size = 0.5) +
  theme(legend.position = "bottom") +
  ggtitle("UMAP by Cell Type (Highlighted if Significant in Mod OR Sev vs Healthy)")


# Label data
labels <- annotation_data %>%
  left_join(highlight_status, by = "predicted.celltype.l2") %>%
  filter(highlight, severity_group %in% c("Moderate", "Severe")) %>%
  mutate(Covariate_map = paste0("severity_group[T.", severity_group, "]")) %>%
  left_join(
    sccoda_results %>%
      filter(Covariate %in% sig_covariates, Final_Parameter != 0) %>%
      select(Covariate, CellType, log2_fold_change),
    by = c("predicted.celltype.l2" = "CellType", "Covariate_map" = "Covariate")
  ) %>%
  filter(!is.na(log2_fold_change))

# Add labels if available
if (nrow(labels) > 0) {
  p <- p + geom_text_repel(data = labels,
    aes(x = final_umap1, y = final_umap2, label = round(log2_fold_change, 2)),
    color = "black", size = 3, fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.5, point.padding = 0.5,
    max.overlaps = Inf, seed = 42
  )
}

# Save plot
ggsave("Manuscript_jan_2025/figures/umap_split_sccoda_highlighted_ModOrSev_labels.png",
       plot = p, width = 10, height = 6, dpi = 300)


#V10
# Recode non-highlighted cell types as "Other"
combined$celltype_highlighted <- ifelse(
  combined$predicted.celltype.l2 %in% highlighted,
  combined$predicted.celltype.l2,
  "Other"
)

# Update color palette: assign grey to "Other"
plot_levels <- c(sort(highlighted), "Other")
plot_colors <- c(setNames(colors, highlighted), Other = "grey85")

# Plot with updated grouping
p <- DimPlot(combined, reduction = "umap", split.by = "severity_group",
             group.by = "celltype_highlighted", cols = plot_colors,
             pt.size = 0.5) +
  theme(legend.position = "bottom") +
  ggtitle("UMAP by Cell Type (Highlighted if Significant in Mod OR Sev vs Healthy)") +
  scale_color_manual(values = plot_colors, breaks = plot_levels)


# Label data
labels <- annotation_data %>%
  left_join(highlight_status, by = "predicted.celltype.l2") %>%
  filter(highlight, severity_group %in% c("Moderate", "Severe")) %>%
  mutate(Covariate_map = paste0("severity_group[T.", severity_group, "]")) %>%
  left_join(
    sccoda_results %>%
      filter(Covariate %in% sig_covariates, Final_Parameter != 0) %>%
      select(Covariate, CellType, log2_fold_change),
    by = c("predicted.celltype.l2" = "CellType", "Covariate_map" = "Covariate")
  ) %>%
  filter(!is.na(log2_fold_change))

# Add labels if available
if (nrow(labels) > 0) {
  p <- p + geom_text_repel(data = labels,
    aes(x = final_umap1, y = final_umap2, label = round(log2_fold_change, 2)),
    color = "black", size = 3, fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.5, point.padding = 0.5,
    max.overlaps = Inf, seed = 42
  )
}

# Save plot
ggsave("Manuscript_jan_2025/figures/umap_split_sccoda_highlighted_ModOrSev_labels.png",
       plot = p, width = 10, height = 6, dpi = 300)



# cell-type funnel ----

#params
logfc_threshold <- log2(1.2)

dataset_paths <- list(
  SevCOVID = "Manuscript_jan_2025/results/SevCOVID_Azimuthl2",
  MilCOVID = "Manuscript_jan_2025/results/MilCOVID_Azimuthl2"
)

comparisons <- list(
  SevCOVID = list(case = "Severe", control = "Healthy", up_col = "Severe_Up", down_col = "Severe_Down"),
  MilCOVID = list(case = "Moderate", control = "Healthy", up_col = "Mild_Up", down_col = "Mild_Down")
)

#run
output <- run_all_comparisons(dataset_paths, comparisons, logfc_threshold)
df_sev <- output$results$SevCOVID
df_mild <- output$results$MilCOVID
seurat_objs <- output$seurat_objects  # Retain these for downstream analysis

celltype_props <- calculate_celltype_proportions(seurat_objs)
celltype_props_normalized <- calculate_condition_normalized_proportions(seurat_objs)

celltype_props_normalized <- celltype_props_normalized %>%
  filter(!(severity_group == "Healthy" & Dataset == "MilCOVID")) %>%
  select(-Dataset) %>%
  group_by(predicted.celltype.l2) %>%
  mutate(Proportion = NormalizedCount / sum(NormalizedCount) * 100) %>%
  ungroup() %>%
  mutate(severity_group = factor(severity_group, levels = c("Severe", "Moderate", "Healthy")))

# Plot
celltype_props_normalized_plot <- ggplot(celltype_props_normalized, aes(x = Proportion, y = predicted.celltype.l2, fill = severity_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Healthy" = "#A6CEE3",   # Light Blue
      "Moderate" = "#FFFACD",  # Light Yellow
      "Severe" = "#FFD700"     # Dark Yellow
    )
  ) +
  labs(x = "Proportion (%)", y = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

ggsave("Manuscript_jan_2025/figures/celltype_props_normalized_plot.png",
       plot = celltype_props_normalized_plot, width = 5, height = 6, dpi = 300)


#used to be called "final_res"
deg_summary_by_cluster <- combine_deg_counts(df_sev, df_mild)

write.csv(deg_summary_by_cluster, "final_deg_table.csv")

message("Done! Wrote final_deg_table.csv with columns: Cluster, Severe_Up, Severe_Down, Mild_Up, Mild_Down.")

# Reshape
abundance_degs <- deg_summary_by_cluster %>%
  select(Cluster, Severe_Up, Severe_Down, Mild_Up, Mild_Down) %>%
  pivot_longer(cols = -Cluster, names_to = "Group", values_to = "Count") %>%
  mutate(
    Condition = case_when(
      grepl("Severe", Group) ~ "Severe",
      grepl("Mild", Group) ~ "Mild"
    ),
    Direction = case_when(
      grepl("Up", Group) ~ "Up",
      grepl("Down", Group) ~ "Down"
    ),
    ColorGroup = paste(Condition, Direction, sep = "_")
  )

abundance_degs <- abundance_degs %>%
  group_by(Cluster) %>%
  mutate(Percent = 100 * Count / sum(Count)) %>%
  ungroup()

# Order clusters by total DEGs if you want
abundance_degs$Cluster <- factor(
  abundance_degs$Cluster,
  levels = deg_summary_by_cluster %>% arrange(desc(Total)) %>% pull(Cluster)
)

# Set Cluster order based on total Moderate (Mild) proportion
cluster_order <- abundance_degs %>%
  filter(Condition == "Mild") %>%
  group_by(Cluster) %>%
  summarise(TotalMild = sum(Percent)) %>%
  arrange(TotalMild) %>%
  pull(Cluster)

# Apply ordering to the factor
abundance_degs <- abundance_degs %>%
  mutate(Cluster = factor(Cluster, levels = cluster_order))

# Custom colors
custom_colors <- c(
  Severe_Up = "#B2182B",     # dark red
  Severe_Down = "#FDBBA0",   # light red
  Mild_Up = "#2166AC",       # dark blue
  Mild_Down = "#B2CDE3"      # light blue
)

# Plot
degs_Severe_moderate_plot <- ggplot(abundance_degs, aes(x = Cluster, y = Percent, fill = ColorGroup)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +  # <- Add this line
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(
    y = "Percentage of DEGs",
    x = "Cluster"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
    )

ggsave("Manuscript_jan_2025/figures/degs_Severe_moderate.png",
       plot = degs_Severe_moderate_plot, width = 4, height = 6, dpi = 300)


shared_cluster_order <- levels(abundance_degs$Cluster)

celltype_props_normalized <- celltype_props_normalized %>%
  mutate(predicted.celltype.l2 = factor(predicted.celltype.l2, levels = shared_cluster_order))

celltype_props_normalized_plot <- ggplot(celltype_props_normalized, aes(x = Proportion, y = predicted.celltype.l2, fill = severity_group)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 33.3, linetype = "dashed", color = "gray40") +  # <- Add this line
  geom_vline(xintercept = 66.6, linetype = "dashed", color = "gray40") +  # <- Add this line
  scale_fill_manual(
    values = c(
      "Healthy" = "#A6CEE3",   # Light Blue
      "Moderate" = "#FFFACD",  # Light Yellow
      "Severe" = "#FFD700"     # Dark Yellow
    )
  ) +
  labs(x = "Proportion (%)", y = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )



combined_plot <- celltype_props_normalized_plot + 
                 degs_Severe_moderate_plot +
                 plot_layout(ncol = 2, widths = c(0.4, 0.6)) & 
                 theme(legend.position = "bottom")

ggsave("Manuscript_jan_2025/figures/combined_composition_and_degs.png",
       plot = combined_plot, width = 10, height = 10, dpi = 300)


# Ensure ordering matches shared_cluster_order
text_data <- celltype_props_normalized %>%
  group_by(predicted.celltype.l2) %>%
  summarise(TotalNormalized = sum(NormalizedCount), .groups = "drop") %>%
  mutate(TotalNormalized = 100*TotalNormalized/3) %>%
  mutate(
    predicted.celltype.l2 = factor(predicted.celltype.l2, levels = shared_cluster_order),
    Label = sprintf("%.2f", TotalNormalized)
  )


text_plot <- ggplot(text_data, aes(y = predicted.celltype.l2, x = 1, label = Label)) +
  geom_text(size = 4.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = "Total cells (%)", y = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )


combined_plot <- text_plot +
                 celltype_props_normalized_plot + 
                 degs_Severe_moderate_plot 
                 plot_layout(ncol = 3, widths = c(0.12,0.4, 0.48)) & 
                 theme(legend.position = "bottom")

combined_plot <- (
  text_plot +
  celltype_props_normalized_plot +
  degs_Severe_moderate_plot +
  plot_layout(ncol = 3, widths = c(0.12, 0.4, 0.48))
) & 
  theme(legend.position = "bottom")


ggsave("Manuscript_jan_2025/figures/combined_composition_and_degs_and_text.png",
       plot = combined_plot, width = 10, height = 10, dpi = 300)


#####################
##### ok back to the heatmaps#
#####################


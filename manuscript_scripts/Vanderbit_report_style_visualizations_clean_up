# abundance v2

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



# abundance V3 ----
# --- Summary of Changes from V2 to V3 ---
# - Introduced `highlight_overall`: cell type is highlighted if *either* Moderate or Severe meets criteria
# - Cell types failing highlight criteria are automatically greyed out in UMAP (`grey85`)
# - Custom color palette (`modified_colors`) created by overriding original colors for non-highlighted types
# - Plot title updated to reflect highlighting based on abundance change criteria
# - Label colors simplified: all text labels are now black for clarity (optional override of red/black logic from V2)
# - Focus of V3: visualize which *cell types overall* are worth attention by dimming irrelevant ones

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


#v4
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
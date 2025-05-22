
# load data ----
combined <- readRDS("pre_processing_test/data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")


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




PCA_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # First scale: point color by Condition
  geom_point(aes(color = Condition), size = 4) +
  scale_color_manual(values = condition_colors, name = NULL,
                     guide = guide_legend(override.aes = list(size = 4))) +

  # Reset color scale
  new_scale_color() +

  # Add line segments with a new color scale (Euclidean distance)
  geom_segment(data = segment_df, 
               aes(x = x, y = y, xend = xend, yend = yend, color = dist),
               inherit.aes = FALSE, linetype = "dashed", size = 1) +
  scale_color_gradient(
  name = "Euclidean distance",
  low = "#DCE6F2",  # light gray-blue
  high = "#08306B"  # deep navy blue
) +
  geom_text_repel(aes(label = label_text), color = "black", size = 4, max.overlaps = Inf) +

  # Axes and theme
  scale_x_continuous(trans = pseudo_log_trans()) +
  scale_y_continuous(trans = pseudo_log_trans()) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",            # <-- this puts legends side by side
    legend.spacing.x = unit(1, "cm"),     # spacing between the two
    legend.title = element_text(hjust = 0.5),
    plot.title = element_blank()
  )


# Save
ggsave(file.path(figures_folder,"pca_plot_custom.png"), PCA_plot, width = 3, height = 6)


p1 <- ggplot(top_PC1, aes(x = reorder(TF, PC1), y = PC1)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top TFs contributing to PC1", x = "TF", y = "Loading") +
  theme_minimal()

p2 <- ggplot(top_PC2, aes(x = reorder(TF, PC2), y = PC2)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top TFs contributing to PC2", x = "TF", y = "Loading") +
  theme_minimal()

# Save
ggsave(file.path(figures_folder,"p1.png"), p1, width = 8, height = 6)
ggsave(file.path(figures_folder,"p2.png"), p2, width = 8, height = 6)
# ==== libraries ====
library(ggplot2)
library(dplyr)    
library(purrr)     # For map functions
library(stringr)   
library(patchwork) 
library(ggridges)  
library(tidyr)
library(Seurat)
library(tibble)
library(reshape2)   # For reshaping matrix to long format
library(gridExtra)  # For arranging multiple heatmaps in a grid
library(pROC)
library(data.table)
library(scales)          # for pseudo_log_trans
library(ggrepel)
library(viridisLite)
library(purrr)
library(ggnewscale)
library(ggbeeswarm) 
library(viridis)

##########################
## FIGURE 2a
##########################

set.seed(1)
figures_folder <- "figures_04_08_2025"
#supp_figures_folder <- "figures_14_06_2025/supp"
dir.create(figures_folder,recursive = TRUE)
#dir.create(supp_figures_folder,recursive = TRUE)

# ==== clean up results from load_all_results.r ====
all_scores_list <- imap(results_preprocessed, ~{
  # For each dataset, iterate through its methods
  dataset_name <- .y
  methods_list <- .x

  # Map over the methods within the current dataset
  imap_dfr(methods_list, ~{
    method_name <- .y
    method_data <- .x

    if (!is.null(method_data) && nrow(method_data) > 0 && "prioritization_score" %in% names(method_data)) {
      # Select the score and add dataset/method identifiers
      method_data %>%
        select(prioritization_score) %>%
        mutate(
          method = method_name,
          dataset = dataset_name
        )
    } else {
      # Return an empty tibble/dataframe if data is missing or invalid
      tibble(prioritization_score = numeric(0), method = character(0), dataset = character(0))
    }
  })
})

# Combine the list of dataframes into one single dataframe
combined_scores_df <- bind_rows(all_scores_list)

# Define a consistent color scheme (adjust colors as needed)
method_colors <- c(
  "Decipher" = "#1f77b4",
  "NicheNet" = "#ff7f0e",
  "LIANA+" = "#2ca02c",
  "Connectome" = "#d62728",
  "NATMI" = "#9467bd"
)

# Optional: Remove datasets/methods with zero rows if any slipped through
combined_scores_df <- combined_scores_df %>% filter(!is.na(prioritization_score))

desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# Convert the 'method' column to a factor with the specified levels.
combined_scores_df <- combined_scores_df %>%
  mutate(method = factor(method, levels = rev(desired_method_order)))

# ==== Box-plot number of interactions ====
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
interaction_counts$method <- factor(interaction_counts$method, levels = rev(desired_method_order))

# Define the breaks you want to keep (excluding the smallest and largest ones)
custom_breaks <- c(1e3, 1e4, 1e5)

# Plot
p <- ggplot(interaction_counts, aes(y = method, x = interaction_count, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.2, size = 2.8, alpha = 0.8) +  # changed from width to height since x/y are flipped
  scale_fill_manual(values = method_colors) +
  scale_x_log10(position = "bottom",breaks = custom_breaks) +  # place x-axis labels at the top
  labs(
    x = "Number of Interactions",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    # X-axis at the top, horizontal text
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 17, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    # Y-axis on left, horizontal (perpendicular to axis)
    axis.text.y = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    text = element_text(face = "bold")
    
  )

ggsave(file.path(figures_folder,"figure_2a.png"), plot = p, width = 5, height = 4, dpi = 300)

write.csv(
  interaction_counts,
  file = file.path(figures_folder, "figure_2a.csv"),
  row.names = TRUE
)

##########################
## FIGURE 2b
##########################

# Generate the Violin Plot
plot_violin_scores_by_method <- ggplot(
  combined_scores_df,
  aes(y = method,                     # flip: categorical axis on Y
      x = prioritization_score,       # numeric axis on X
      fill = method)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7,
    scale = "width",
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) +

  # same colour palette / order
  scale_fill_manual(values = method_colors, breaks = rev(desired_method_order)) +

  # numeric axis at the top, pseudo-log like in your box plot
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(sigma = 0.1, base = 10),
    breaks = c(-10, -1, 0, 1, 10, 100),
    position = "bottom",
    name = "Prioritization Score"
  ) +

  # categorical axis order, but hide ticks and title
  scale_y_discrete(limits = rev(desired_method_order), name = NULL) +

  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    # X-axis (now across the top)
    axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 17, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    # Suppress Y-axis ticks/labels/title
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(face = "bold")
  ) +
  labs(fill = "Method")


# Save the Plot
output_filename_violin <- "figure_2b.png"
ggsave(file.path(figures_folder,output_filename_violin), plot = plot_violin_scores_by_method, width = 3.5, height = 4, dpi = 300)

write.csv(
  combined_scores_df,
  file = file.path(figures_folder, "figure_2b.csv"),
  row.names = TRUE
)


##########################
## FIGURE 2c
##########################

n_top <- 100
method_names_expected <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")
all_intersection_data <- list()

print(paste("Processing", length(results_for_correlation), "datasets for full overlap analysis..."))
dataset_names_to_process <- names(results_for_correlation)

for (ds in dataset_names_to_process) {
    print(paste("Calculating overlaps for:", ds))
    dataset_results <- results_for_correlation[[ds]]

    if (is.null(dataset_results) || !is.list(dataset_results) || length(dataset_results) == 0) {
        warning(paste("Skipping dataset", ds, "- invalid entry."))
        next
    }

    list_input_current_ds <- list()
    available_methods_in_ds <- names(dataset_results)
    for (method_name in method_names_expected) {
        if (method_name %in% available_methods_in_ds) {
             list_input_current_ds[[method_name]] <- tryCatch({
                  getSet(dataset_results[[method_name]], n_top)
             }, error = function(e) {
                  warning(paste("Error in getSet for", method_name, ds, ":", e$message)); return(character(0))
             })
        } else {
             list_input_current_ds[[method_name]] <- character(0)
        }
    }

    # *** Use the new function here ***
    dataset_summary <- tryCatch({
         calculate_all_intersections(list_input_current_ds) # Call the new function
    }, error = function(e) {
         warning(paste("Error calculating all intersections for", ds, ":", e$message)); return(NULL)
    })

    if (!is.null(dataset_summary)) {
        dataset_summary$Dataset <- ds
        all_intersection_data[[ds]] <- dataset_summary
    } else {
         warning(paste("No intersection summary generated for dataset:", ds))
    }
}

# --- Combine Results and Plot ---
if (length(all_intersection_data) > 0) {
    combined_intersection_df <- bind_rows(all_intersection_data)

    # --- Define Logical Order for Plotting (by Degree, then Name) ---
    # Calculate degree (number of methods in the intersection)
    combined_intersection_df <- combined_intersection_df %>%
      mutate(Degree = sapply(strsplit(Intersection_Name, " & "), length))

    # Order the dataframe to get factor levels easily
    ordered_df <- combined_intersection_df %>%
      arrange(desc(Degree), Intersection_Name) # Order by Degree (high to low), then alphabetically

    # Get the unique names in the desired order
    intersection_order_all <- unique(ordered_df$Intersection_Name)

    # Apply the factor levels
    combined_intersection_df <- combined_intersection_df %>%
        mutate(Intersection_Name = factor(Intersection_Name, levels = intersection_order_all)) %>%
        filter(!is.na(Intersection_Name)) # Clean up just in case

    if(nrow(combined_intersection_df) == 0) {
        stop("No valid intersection data remained after ordering.")
    }

    # Save the combined data used for the plot
    saveRDS(combined_intersection_df, file = file.path(figures_folder,"overlap_all_intersections_data.rds"))
    print("Full overlap data saved to figures/overlap_all_intersections_data.rds")

} else {
    warning("No intersection data was successfully calculated for any dataset. Cannot create plot.")
}


summary_df <- combined_intersection_df %>%
  # 2. explode "A & B & C" into separate rows
  separate_rows(Intersection_Name, sep = " & ") %>%
  # 3. rename the exploded column to "method"
  rename(method = Intersection_Name) %>%
  # 4. sum up counts by dataset, method, and degree
  group_by(Dataset, method, Degree) %>%
  summarise(total_Count = sum(Count), .groups = "drop")

summary_df$method <- factor(summary_df$method,levels=desired_method_order)

degree_labels <- c("1" = "unique", "2" = "two-way", "3" = "three-way", "4" = "four-way", "5" = "five-way")
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
ggsave(filename = file.path(figures_folder, "figure_2c.png"),
       plot = p, width = 8, height = 4)

write.csv(
  summary_df,
  file = file.path(figures_folder, "figure_2c.csv"),
  row.names = TRUE
)
##########################
## FIGURE 2d
##########################
spearman_matrices <- list()
k_matrices <- list()

for (dataset in names(results_for_correlation)) {
  # Compute correlation & search space
  interaction_results_correlation_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(results_for_correlation[[dataset]])

  # Extract Spearman correlation matrix
  spearman_matrices[[dataset]] <- interaction_results_correlation_search_space$spearman

  # Extract k-matrix
  k_matrices[[dataset]] <- interaction_results_correlation_search_space$k_matrix
}

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

#search space -----
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

# Split method pairs into two separate columns for hierarchical labeling
combined_k_df <- combined_k_df %>%
  mutate(Method1 = sub("-.*", "", Method_Pair),  # Extract first method
         Method2 = sub(".*-", "", Method_Pair))  # Extract second method

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

# --- Define colors and desired order ---
color_k_value <- "#E69F00" # Gold/Orange
color_spearman <- "#56B4E9" # Light Blue
# Make sure this exists in your environment
desired_method_order <- c("Decipher", "NicheNet", "LIANA+", "NATMI", "Connectome")

# --- 1. Filter out self-comparisons AND Prepare Method Columns ---
# Ensure Method1/Method2 exist, using Var1/Var2 if necessary
if (!all(c("Method1", "Method2") %in% names(combined_df))) {
   combined_df <- combined_df %>% mutate(Method1 = Var1, Method2 = Var2)
   warning("Recreated Method1/Method2 columns from Var1/Var2.", call. = FALSE)
}

plot_data <- combined_df %>%
  filter(Method1 != Method2) %>% # Keep only pairs where methods are different
  # Ensure methods are within the desired order list
  filter(Method1 %in% desired_method_order & Method2 %in% desired_method_order)

if(nrow(plot_data) == 0) {
  stop("Filtering removed all data. Check input 'combined_df' and filtering conditions.")
}


# --- 2. & 3. Order Facets (Method1) and X-axis (Method2) ---
# Convert to factors WITH the desired levels
plot_data$Method1 <- factor(plot_data$Method1, levels = desired_method_order)
plot_data$Method2 <- factor(plot_data$Method2, levels = desired_method_order)

# proper plot below #######

## 1. define your orders
method_order  <- c("Decipher","NicheNet","LIANA+","NATMI","Connectome")
dataset_order <- c(
  "5yr_pic","bcg","cord_pic","covid",
  "erp","lupus","sepsis","tnbc",
  "cz_influenza","cz_hpap_t1d_islets","cz_hnscc_hpv",
  "cz_human_kidney_v1.5","cz_cf_bronchial_biopsy",
  "SevCOVID_Azimuthl2","MilCOVID_Azimuthl2"
)

## 2. build the “full grid” of every (Method1,Method2) × (tile_row,tile_col)
grid_df <- expand.grid(
  Method1  = method_order,
  Method2  = method_order,
  tile_row = 1:4,
  tile_col = 1:4,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Method1 = factor(Method1, levels = method_order),
    Method2 = factor(Method2, levels = method_order),
    # figure out which dataset should live here (1–15 → real names; 16 → NA)
    dataset_index = (tile_row - 1)*4 + tile_col,
    Dataset       = ifelse(dataset_index <= length(dataset_order),
                           dataset_order[dataset_index],
                           NA_character_)
  )

## 3. left-join your real values on to that grid
heat_df <- grid_df %>%
  left_join(combined_df, by = c("Method1","Method2","Dataset")) %>%
  mutate(
    big_row = match(Method1, method_order) - 1,
    big_col = match(Method2, method_order) - 1,
    x       = big_col*4 + tile_col,
    y       = (4 - tile_row) + big_row*4
  )

# 4. build a small data.frame of the 5×5 block‐centres
border_df <- expand.grid(
  Method1 = method_order,
  Method2 = method_order,
  stringsAsFactors = FALSE
) %>% 
  mutate(
    big_row = match(Method1, method_order) - 1,
    big_col = match(Method2, method_order) - 1,
    # centre of each 4×4 block:
    x = big_col*4 + 2.5,
    y = big_row*4 + 1.5
  )

## 5. plot—with one fill scale for Spearman (upper triangle)
##    and a second for k_value (lower triangle)
p <- ggplot() +
  # Diagonal override
  geom_tile(
    data = heat_df %>% filter(big_row == big_col),
    aes(x = x, y = y),
    fill = "lightgray",
    inherit.aes = FALSE
  ) +
  # upper tri: Spearman
  geom_tile(
    data = heat_df %>% filter(big_row < big_col),
    aes(x, y, fill = Spearman)
  ) +
  scale_fill_gradient2(
    name    = "Spearman",
    low     = "#b2182b", mid = "white", high = "#008837",
    limits  = c(-1,1),
    na.value= "grey80"
  ) +

  # start a fresh fill mapping
  new_scale_fill() +

  # lower tri: k_value
  geom_tile(
    data = heat_df %>% filter(big_row > big_col),
    aes(x, y, fill = k_value)
  ) +
  scale_fill_viridis_c(
    name    = "Search-space\n(k value)",
    option  = "B", end = 0.95,
    limits  = c(100, max(heat_df$k_value, na.rm = TRUE)),
    na.value= "grey80"
  ) +
  # --- border layer ---
  geom_tile(
    data      = border_df,
    aes(x, y),
    width     = 4,        # span 4 tiles
    height    = 4,
    fill      = NA,       # transparent
    color     = "white",  # white outline
    size      = 0.8,      # line thickness
    inherit.aes = FALSE
  ) +
  # tidy up axes so each 4×4 block is labelled by method
  coord_fixed() +
  scale_x_continuous(
    expand = c(0,0),
    breaks = (0:4)*4 + 2.5,
    labels = method_order,
    position = "bottom"
  ) +
  scale_y_reverse(
    expand = c(0,0),
    breaks = (0:4)*4 + 1.5,
    labels = method_order
  ) +

  theme_minimal(base_size = 12) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(
      size   = 17,
      face   = "bold",
      angle  = 45,
      hjust  = 1              # right-justify along the diagonal
    ),
    axis.text.y  = element_text(face = "bold", size = 17,vjust=1),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(figures_folder,"figure_2d.png"),p)
write.csv(
  heat_df,
  file = file.path(figures_folder, "figure_2d_heatmap.csv"),
  row.names = TRUE
)
write.csv(
  border_df,
  file = file.path(figures_folder, "figure_2d_border.csv"),
  row.names = TRUE
)


##########################
## FIGURE 2e
##########################
#load cytosig data
cytosig_significance   <- list() 

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
  
  # Load mapping table
  mapping_table <- read.csv(file.path(reference_filepath,"cytosig_mapping_table_ligands_genes.csv"),header=TRUE)

  # Process Cytosig significance
  cytosig_significance[[ds]] <- summarizeZScores(z_score_files, z_score_folder, mapping_table)
}

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


#==== ROC curves ====
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

#ggbeeswarm
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

# Filter flagged points n_true lt 10
flagged_points <- filter(results_df, flag == "*")

# Plot
results_df$flagged <- results_df$flag == "*"

#color line by variance
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

# --- Prepare Data: Order Methods ---
available_methods_plot <- unique(results_df$method)
valid_order_plot <- intersect(desired_method_order, available_methods_plot)

# Filter data if some methods are not in the desired order list
if(length(valid_order_plot) < length(available_methods_plot)){
    warning("Some methods found in 'results_df$method' are not in 'desired_method_order'. Plotting only those included in the defined order and filtering data.")
    results_df <- results_df %>% filter(method %in% valid_order_plot)
}
# Stop if no data remains after filtering
if(nrow(results_df) == 0) {
    stop("No data remains after filtering 'results_df' based on 'desired_method_order'. Check method names.")
}

# Convert 'method' column to an ordered factor
results_df$method <- factor(results_df$method, levels = valid_order_plot)


# --- Dynamically Create Line Color Mapping ---
# Identify the unique values in the 'line_color' column to map them
line_color_values_in_data <- sort(unique(results_df$line_color)) # Sort for consistent mapping order
if (length(line_color_values_in_data) == 0) {
    stop("No values found in the 'line_color' column.")
} else if (length(line_color_values_in_data) > 2) {
    warning(paste("Found more than 2 unique values in 'line_color':",
                  paste(line_color_values_in_data, collapse=", "),
                  "- Only the first two will be mapped to Gold/Blue."))
    # Map only the first two values found
    line_color_map <- setNames(c(color_k_value, color_spearman), line_color_values_in_data[1:2])
} else if (length(line_color_values_in_data) == 1) {
     warning(paste("Found only 1 unique value in 'line_color':", line_color_values_in_data[1],
                  "- Mapping it to Gold. All lines will have the same color."))
     line_color_map <- setNames(color_k_value, line_color_values_in_data[1])
} else {
   # Exactly two values found, map them
   line_color_map <- setNames(c(color_k_value, color_spearman), line_color_values_in_data)
}


####################
# FIGURE 2e
####################
p_updated <- ggplot(results_df, aes(x = method, y = value)) +

  # 1. Lines - Mapped to 'line_color', using first color scale
  #geom_line(aes(group = dataset, color = line_color), size = 1, alpha = 0.6) +
  geom_line(aes(group = dataset, color = "lightgray"), size = 1, alpha = 0.6) +

  # 2. Boxplot - Neutral colors (Plot after lines, before points?)
  geom_boxplot(outlier.shape = NA, width = 0.25, alpha = 0.4, color = "black", fill = "lightgray") +

  # 3. Beeswarm Points - Mapped to 'method', using second color scale
  geom_beeswarm(
    aes(color = method),   # Color points by method
    size = 4.0,            # Adjust point size if needed
    cex = 5,               # Use 'size', not 'cex' for ggplot point size control
    priority = "density",
    groupOnX = TRUE,        # Ensure beeswarm groups correctly on discrete x-axis
    alpha = 0.5
  ) +
  # Second color scale specifically for the beeswarm points
  scale_color_manual(
    # name = "Method",     # Optional legend name
    values = method_colors, # Use the method colors for points
    guide = "none"         # Hide legend for points
  ) +

  # Median line (same as before)
  stat_summary(fun = median, geom = "segment",
               aes(xend = after_stat(x), yend = after_stat(y)),
               size = 2.5, color = "black") +

  # Horizontal line at 0.5 (same as before)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +

  # Ensure X-axis order respects the factor levels
  scale_x_discrete(limits = levels(results_df$method), # Use levels from the factor
                   drop = FALSE) + # Prevent dropping unused levels

  # Labels and Theme
  labs(y = "AUROC target prediction", x = NULL) + # x = NULL removes x-axis title
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=17,face="bold"),
        panel.grid = element_blank()) # Rotate labels

# --- Save the Plot ---
output_filename_updated <- file.path(figures_folder,"figure_2e.png")
ggsave(output_filename_updated, plot = p_updated, width = 4, height = 8, dpi = 300) # Adjusted width slightly

write.csv(
  results_df,
  file = file.path(figures_folder, "figure_2e.csv"),
  row.names = TRUE
)

####################
# FIGURE 2f
####################

# Initialize empty lists to store the results
all_runs_dat_1 <- list()
all_runs_dat_2 <- list()

N <- 100
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

# Calculate total count per interaction across all cell types
interaction_totals <- interaction_consistency %>%
  group_by(interaction) %>%
  summarise(total_count = sum(count), .groups = "drop")

# Join the total_count back into your interaction_consistency table
interaction_consistency <- interaction_consistency %>%
  left_join(interaction_totals, by = "interaction")

# Base R reordering 
interaction_consistency$interaction <- factor(
  interaction_consistency$interaction,
  levels = interaction_totals$interaction[order(interaction_totals$total_count,decreasing = FALSE)]
)

top20_interactions <- interaction_consistency %>%
  arrange(desc(proportion)) %>%
  slice_head(n = 30)

top20_mapped <- top20_interactions %>%
  mutate(Cells = case_when(
    grepl("^B_", Cells)                     ~ "B",
    grepl("^CD14_plus_Monocytes", Cells)    ~ "CD14+ M",
    grepl("^CD4_T", Cells)                  ~ "CD4 T",
    grepl("^CD8", Cells)                    ~ "CD8 T",
    TRUE                                    ~ Cells  # fallback if unmatched
  ))
  
p <- ggplot(top20_mapped, aes(x = Cells, y = interaction, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", name = "Count") +  # Viridis scale
  scale_y_discrete(position = "right") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 17,face="bold"),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 13), 
    panel.grid = element_blank(),
    legend.position = "bottom"     # move legend to bottom
  )

ggsave(file.path(figures_folder,"figure_2f.png"), p, width = 4.2, height = 8)
write.csv(
  top20_mapped,
  file = file.path(figures_folder, "figure_2f.csv"),
  row.names = TRUE
)

####################
# FIGURE 2g
####################
#logic to clean this plot sits inside the function, not optimal but ok
#result figure written in sample_analysis/validity/data/for_plotting/sample_1
plotDecipherPrioritizedMap("sample_analysis/validity/data/for_plotting",top_n=4,dataset_name="sample_1", abs_decipher_plot_limit = 20,width=21,height=9)

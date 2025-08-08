# ── libraries ────────────────────────────────────────────────────────────────
library(dplyr)
library(stringr)
library(gplots)          # heatmap.2
library(RColorBrewer)
library(matrixStats)     # rowVars()

set.seed(1)
figures_folder <- "figures_04_08_2025"

##### supplementary figure 2 ########

# ── 1. read the full interaction-potential matrix ───────────────────────────
output_data_filepath <- "results/cord_pic/data"
pot_file <- file.path(output_data_filepath,
                      "interaction_potential_by_clusters.rds")

all_potentials <- readRDS(pot_file)

# ── 2. slice to B-cell receiver only ────────────────────────────────────────
receiver_cluster <- "B"    
mat <- all_potentials[[receiver_cluster]]

# ── 3. tidy up NA / zero-variance rows --------------------------------------
# keep LR pairs that have at least two non-NA values (needed for correlation)
keep <- rowSums(!is.na(mat)) >= 2 & rowVars(mat, na.rm = TRUE) > 0
mat  <- mat[keep, , drop = FALSE]

# ── 4. correlation matrix (pairwise, Spearman) ------------------------------
cor_mat <- cor(t(mat),                             # rows → variables
               use = "pairwise.complete.obs",
               method = "spearman")

# cor() may still return NA if two rows never overlap → set those to 0
cor_mat[is.na(cor_mat)] <- 0

# ── 5. receptor vector & colour side-bar ------------------------------------
receptor_vec <- str_extract(rownames(cor_mat), "(?<=-).*$")
## how often does each receptor show up?
rec_freq <- sort(table(receptor_vec), decreasing = TRUE)

top_k    <- 10                       # <- tweak here
keepers  <- names(rec_freq)[seq_len(top_k)]

set2_big <- colorRampPalette(brewer.pal(8, "Set2"))
receptor_colours <- setNames(set2_big(length(keepers)), keepers)

# map every row’s receptor to a colour (grey for non-keepers)
row_side <- ifelse(receptor_vec %in% keepers,
                   receptor_colours[receptor_vec],
                   "grey80")
col_side <- row_side  

# ── 6. dendrogram (distance = 1-ρ) -------------------------------
hc   <- hclust(as.dist(1 - cor_mat), method = "average")   
dend <- as.dendrogram(hc)                                 
ord  <- hc$order                                           

png(file.path(figures_folder,"supp_figure_2.png"),
    width  = 3000,     
    height = 3000,     
    res    = 300)      
    
    # ── 7. draw with heatmap.2() -------------------------------------------------
heatmap.2(
  cor_mat,
  Rowv = dend, Colv = dend,
  dendrogram = "both",
  trace = "none",
  col   = colorRampPalette(c("royalblue4","white","firebrick3"))(101),
  RowSideColors = row_side,
  ColSideColors = col_side,
  key.title = "Spearman ρ",
  key       = TRUE,
  labRow = rownames(cor_mat),
  labCol = colnames(cor_mat),
  cexRow = 0.6,                
  cexCol = 0.6,
  density.info = "none",
  margins = c(5,5),
  main = "LR interaction-potential correlation – B-cells"
)
dev.off()



## ---- 3. save to csv -----------------------------------------------------
cor_mat_ordered <- cor_mat[ord, ord]

write.csv(
  cor_mat_ordered,
  file = file.path(figures_folder,
                   "supp_figure_2.csv"),
  row.names = TRUE
)

################## Supplementary Figure 1
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
L_set <- readRDS(file.path(dataset_path,"data/L_set.rds"))
# assume your helper bubble function is already loaded:
# source("R/plotBubble.R")  

# 1) load & preprocess all three methods for one dataset
dataset_path <- "results/cord_pic"
decipher_raw    <- readRDS(file.path(dataset_path, "data/decipher_scores_by_cluster.rds"))
nichenet_raw    <- readRDS(file.path(dataset_path, "nichenet/data/prior_table_all_clusters.rds"))
liana_raw       <- read.csv(file.path(dataset_path, "liana/data/liana_p_interaction_results.csv"),
                            row.names=1, check.names=FALSE)

# your existing preprocessors
decipher_df <- preProcessDecipher(decipher_raw)   # from your load_all() script
nichenet_df <- preProcessNicheNet(nichenet_raw)
liana_df    <- preProcessLIANA(liana_raw)

liana_df <- liana_df %>%                # ← your original LIANA table
  left_join(
    L_set %>% select(interaction, ligand, receptor),
    by = "interaction"
  ) %>%
  relocate(ligand, receptor, .after = interaction)  # nicer column order

# pick top‐8 (4 pos, 4 neg) interactions per method
select_top_n <- function(df, top_n=4, score_col){
  score_col <- rlang::ensym(score_col)
  df %>%
    group_by(interaction) %>%
    summarize(score = mean(!!score_col, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(score)) %>%
    # take top_n highest AND top_n lowest
    slice(c(1:top_n, (n()-top_n+1):n())) %>%
    pull(interaction)
}

decipher_top <- select_top_n(decipher_df, 4, scaled_score)
nichenet_top <- select_top_n(nichenet_df, 4, prioritization_score)
liana_top    <- select_top_n(liana_df,    4, prioritization_score)

# focus on one receiver cell at a time; for example “Monocyte”
receiver <- "Mono"


##### second version
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ---- 1. pick top +/- per receiver -------------------------------------------
# top_n here is PER receiver; split_by_direction=TRUE gives top_n pos & top_n neg per receiver
select_top_per_receiver <- function(df, top_n = 4, score_col, split_by_direction = TRUE) {
  score_col <- rlang::ensym(score_col)
  tmp <- df %>%
    mutate(.score = !!score_col,
           .sign  = if_else(.score >= 0, "pos", "neg")) %>%
    group_by(receiver, interaction, .sign) %>%
    summarise(score = mean(.score, na.rm = TRUE), .groups = "drop")

  if (split_by_direction) {
    top_pos <- tmp %>%
      filter(.sign == "pos") %>%
      group_by(receiver) %>%
      slice_max(order_by = score, n = top_n, with_ties = FALSE) %>%
      ungroup()

    top_neg <- tmp %>%
      filter(.sign == "neg") %>%
      group_by(receiver) %>%
      slice_min(order_by = score, n = top_n, with_ties = FALSE) %>%
      ungroup()

    bind_rows(top_pos, top_neg) %>%
      distinct(receiver, interaction)
  } else {
    tmp %>%
      group_by(receiver) %>%
      slice_max(order_by = abs(score), n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      distinct(receiver, interaction)
  }
}

# ---- 2. coerce each method to a "Decipher-like" plotting schema -------------
# Expected output columns:
# interaction, ligand, receptor, sender, receiver,
# method_score (center color), method_score_abs (for bubble size),
# ligand.diff.expr, receptor.diff.expr,
# ligand.frac, receptor.frac
#
# NicheNet has rich LFC/pct info; LIANA+ often has just scores.

coerce_nichenet_schema <- function(nichenet_df) {
  nichenet_df %>%
    mutate(
      method_score     = coalesce(scaled_activity, prioritization_score, scaled_score),
      method_score_abs = abs(method_score),
      ligand.diff.expr   = lfc_ligand,
      receptor.diff.expr = lfc_receptor,
      # use % expressed as a proxy for bubble stroke/size thresholds
      ligand.frac   = pct_expressed_sender/100,
      receptor.frac = pct_expressed_receiver/100
    ) %>%
    select(interaction, ligand, receptor, sender, receiver,
           method_score, method_score_abs,
           ligand.diff.expr, receptor.diff.expr,
           ligand.frac, receptor.frac)
}

coerce_liana_schema <- function(liana_df) {
  # LIANA often lacks per-gene LFC; keep NA for color panels on the sides,
  # rely on method_score for center, and scaled_score to gate bubble size a bit.
  liana_df %>%
    mutate(
      method_score     = coalesce(scaled_score, prioritization_score),
      method_score_abs = abs(method_score),
      ligand.diff.expr   = NA_real_,
      receptor.diff.expr = NA_real_,
      ligand.frac   = NA_real_,
      receptor.frac = NA_real_
    ) %>%
    select(interaction, ligand, receptor, sender, receiver,
           method_score, method_score_abs,
           ligand.diff.expr, receptor.diff.expr,
           ligand.frac, receptor.frac)
}

# plot limits helper
.range_or_zero <- function(x) {
  if (all(is.na(x))) list(min = 0, max = 0) else list(min = min(x, na.rm=TRUE), max = max(x, na.rm=TRUE))
}

plotMethodPrioritizedMap <- function(method_name,
                                     df_std,          # output of coerce_*_schema()
                                     top_tbl,         # receiver/interaction table from select_top_per_receiver()
                                     selected_receivers = NULL,  # optional subset of receivers
                                     abs_center_limit = NULL,    # optional symmetric cap for center color
                                     width_cm = 21, height_cm = 11,
                                     out_prefix = NULL) {

  # keep only selected receivers/interactions
  df_plot <- df_std %>%
    inner_join(top_tbl, by = c("receiver","interaction"))

  if (!is.null(selected_receivers)) {
    df_plot <- df_plot %>% filter(receiver %in% selected_receivers)
  }

  # bubbles: size from method_score_abs (center), and from ligand/receptor.frac on sides if present
  df_plot <- df_plot %>%
    mutate(
      size_center   = if_else(method_score_abs > quantile(method_score_abs, 0.1, na.rm=TRUE), 1, NA_real_),
      stroke_center = 0.5,
      stroke_ligand = if_else(!is.na(ligand.frac) & ligand.frac > 0.05, 0.5, NA_real_),
      size_ligand   = if_else(!is.na(ligand.frac) & ligand.frac > 0.05, ligand.frac, NA_real_),
      stroke_recept = if_else(!is.na(receptor.frac) & receptor.frac > 0.05, 0.5, NA_real_),
      size_recept   = if_else(!is.na(receptor.frac) & receptor.frac > 0.05, receptor.frac, NA_real_)
    )

  # color ranges
  lim_lig  <- .range_or_zero(df_plot$ligand.diff.expr)
  lim_rec  <- .range_or_zero(df_plot$receptor.diff.expr)

  if (is.null(abs_center_limit)) {
    lim_ctr <- .range_or_zero(df_plot$method_score)
  } else {
    # cap center colors symmetrically (like your Decipher code)
    eps <- 0.01 * abs_center_limit
    df_plot <- df_plot %>%
      mutate(method_score = pmax(pmin(method_score,  abs_center_limit - eps),
                                 -abs_center_limit + eps))
    lim_ctr <- list(min = -abs_center_limit, max = abs_center_limit)
  }

  # ---- three panels ----
  p_lig <- plotBubble(
    df = df_plot,
    x_var = "receiver",
    color.var = "ligand.diff.expr",
    size.var  = "size_ligand",
    stroke.var= "stroke_ligand",
    plot.position = "left",
    col.min.val = lim_lig$min, col.max.val = lim_lig$max,
    plot.title = "Ligand",
    x_lab = "SCT", y_lab = "Interaction"
  )

  p_ctr <- plotBubble(
    df = df_plot,
    x_var = "receiver",
    color.var = "method_score",
    size.var  = "size_center",
    stroke.var= "stroke_center",
    plot.position = "middle",
    col.min.val = lim_ctr$min, col.max.val = lim_ctr$max,
    plot.title = paste0(method_name, " score"),
    x_lab = "RCT", y_lab = ""
  )

  p_rec <- plotBubble(
    df = df_plot,
    x_var = "receiver",
    color.var = "receptor.diff.expr",
    size.var  = "size_recept",
    stroke.var= "stroke_recept",
    plot.position = "right",
    col.min.val = lim_rec$min, col.max.val = lim_rec$max,
    plot.title = "Receptor",
    x_lab = "RCT", y_lab = ""
  )

  # compose & save (optional)
  composed <- p_lig + p_ctr + p_rec + patchwork::plot_layout(widths = c(2,1,1))

  if (!is.null(out_prefix)) {
    dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
    png(paste0(out_prefix, ".png"), width = width_cm, height = height_cm, units = "cm", res = 600)
    print(composed)
    dev.off()
    write.csv(df_plot, paste0(out_prefix, ".csv"), row.names = FALSE)
  }

  composed
}

# already loaded:
# dataset_path <- "results/cord_pic"
# decipher_df <- preProcessDecipher(decipher_raw)   # you already have this
# nichenet_df <- preProcessNicheNet(nichenet_raw)
# liana_df    <- preProcessLIANA(liana_raw)

# choose how many per receiver (4 pos + 4 neg)
top_n_per_receiver <- 4

# DECIPHER (you already have a dedicated function, but for completeness:)
decipher_top_tbl <- select_top_per_receiver(decipher_df, top_n = top_n_per_receiver, score_col = scaled_score, split_by_direction = TRUE)

# NICHE NET
nichenet_std <- coerce_nichenet_schema(nichenet_df)
nichenet_top_tbl <- select_top_per_receiver(
  nichenet_std,
  top_n = top_n_per_receiver,
  score_col = method_score,            # uses scaled_activity/prioritization
  split_by_direction = TRUE
)

p_nn <- plotMethodPrioritizedMap(
  method_name = "NicheNet",
  df_std = nichenet_std,
  top_tbl = nichenet_top_tbl,
  selected_receivers = NULL,           # or e.g. c("B","Mono","NK",...)
  abs_center_limit = NULL,             # or give a symmetric cap like 2
  out_prefix = file.path("figures_04_08_2025", "supp_fig1_nichenet")
)

# LIANA+
liana_std <- coerce_liana_schema(liana_df)
liana_top_tbl <- select_top_per_receiver(
  liana_std,
  top_n = top_n_per_receiver,
  score_col = method_score,
  split_by_direction = TRUE
)

p_li <- plotMethodPrioritizedMap(
  method_name = "LIANA+",
  df_std = liana_std,
  top_tbl = liana_top_tbl,
  selected_receivers = NULL,
  abs_center_limit = NULL,
  out_prefix = file.path("figures_04_08_2025", "supp_fig1_liana")
)

# If you also want a Decipher panel with the *same* logic:
# First coerce Decipher to the same schema quickly:
decipher_std <- decipher_df %>%
  mutate(
    method_score     = coalesce(scaled_score, prioritization_score),
    method_score_abs = abs(method_score),
    ligand.diff.expr = NA_real_,   # you have this in your full pipeline; put NA if not available
    receptor.diff.expr = NA_real_,
    ligand.frac = NA_real_,
    receptor.frac = NA_real_
  ) %>%
  select(interaction, ligand, receptor, sender, receiver,
         method_score, method_score_abs,
         ligand.diff.expr, receptor.diff.expr, ligand.frac, receptor.frac)

p_dc <- plotMethodPrioritizedMap(
  method_name = "Decipher",
  df_std = decipher_std,
  top_tbl = decipher_top_tbl,
  selected_receivers = NULL,
  abs_center_limit = NULL,
  out_prefix = file.path("figures_04_08_2025", "supp_fig1_decipher")
)
png(file.path(figures_folder,"supp_figure_1.png"),width  = 3000,     
    height = 6000,     
    res    = 300)
# Arrange three method panels vertically if you want one composite figure:
(p_dc / p_nn / p_li) + plot_annotation(
  title = "Comparison of systems-level CCC maps for Cord PIC vs Unstimulated Cord CBMC",
  subtitle = "Top eight (four positive, four negative) LR interactions per receiver; left=ligand stats, middle=method score, right=receptor stats."
)
dev.off()

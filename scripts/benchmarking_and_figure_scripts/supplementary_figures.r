# ── libraries ────────────────────────────────────────────────────────────────
library(dplyr)
library(stringr)
library(gplots)         
library(RColorBrewer)
library(matrixStats)     
library(tidyr)
library(ggplot2)
library(patchwork)


#functions

# ---- 1. pick top +/- per receiver -------------------------------------------
# top_n here is PER receiver; split_by_direction=TRUE gives top_n pos & top_n neg per receiver
select_top_per_receiver <- function(df, top_n = 4, score_col, split_by_direction = TRUE) {
  score_col <- rlang::ensym(score_col)
  tmp <- df %>%
    mutate(.score = !!score_col,
           .sign  = if_else(.score >= 0, "pos", "neg")) %>%
    group_by(receiver, sender, interaction, .sign) %>%
    summarise(score = mean(.score, na.rm = TRUE), .groups = "drop") #%>%
    #distinct(receiver, sender, interaction)  

  if (split_by_direction) {
    top_pos <- tmp %>%
      filter(.sign == "pos") %>%
      group_by(receiver,sender) %>%
      slice_max(order_by = score, n = top_n, with_ties = FALSE) %>%
      ungroup()

    top_neg <- tmp %>%
      filter(.sign == "neg") %>%
      group_by(receiver,sender) %>%
      slice_min(order_by = score, n = top_n, with_ties = FALSE) %>%
      ungroup()

    bind_rows(top_pos, top_neg) %>%
      distinct(receiver, interaction)
  } else {
    tmp %>%
      group_by(receiver,sender) %>%
      slice_max(order_by = abs(score), n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      distinct(receiver, interaction)
  }
}

coerce_nichenet_schema <- function(nichenet_df) {
  nichenet_df %>%
    mutate(
      method_score     = prioritization_score,
      method_score_abs = abs(method_score),
      ligand.diff.expr   = scaled_p_val_adapted_ligand,
      receptor.diff.expr = scaled_p_val_adapted_receptor,
      # use % expressed as a proxy for bubble stroke/size thresholds
      ligand.frac   = scaled_avg_exprs_ligand,
      receptor.frac = scaled_avg_exprs_receptor
    ) %>%
    select(interaction, ligand, receptor, sender, receiver,
           method_score, method_score_abs,
           ligand.diff.expr, receptor.diff.expr,
           ligand.frac, receptor.frac)
}

invert_scores <- function(score, eps = .Machine$double.eps) {
  -log10(score + eps)
}

coerce_liana_schema <- function(liana_df) {
  # LIANA often lacks per-gene LFC; keep NA for color panels on the sides,
  # rely on method_score for center, and scaled_score to gate bubble size a bit.
  liana_df %>%
    mutate(
      method_score     = prioritization_score,
      method_score_abs = abs(method_score),
      ligand.diff.expr   = invert_scores(ligand_padj),
      receptor.diff.expr = invert_scores(receptor_padj),
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

LocalPlotBubble <- function(df,color.var,size.var,stroke.var,plot.position,col.min.val,col.max.val,plot.title,x_lab,y_lab,x_var,legend.title = NULL){
  legend.title <- legend.title %||% color.var      

  this.plot <- ggplot(df,aes_string(y="interaction", x=x_var,fill = color.var)) +
    geom_point(aes_string(size = size.var,stroke=stroke.var), shape = 21) +
    labs(x = NULL, y = y_lab) +
    theme_bw()+
    scale_fill_gradient2(
        low  = "blue", mid = "white", high = "red", midpoint = 0,
        name   = legend.title,  
        limits   = c(col.min.val, col.max.val),
        n.breaks = 3,                             # ← FEWER TICKS
        guide = guide_colorbar(
            title.position = "top", title.hjust = .5,
            barwidth  = unit(3, "cm"),          
            barheight = unit(.4, "cm")
        )
        )+
    ggtitle(label = plot.title)+ guides(size = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size=10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, 'cm'),           # Smaller legend keys
          legend.text = element_text(size = 10),        # Smaller legend text
          plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 2, unit = "pt"),
          plot.title = element_text(hjust = 0.5))

  if(plot.position == "middle"){
    this.plot <- this.plot + theme(axis.text.y=element_blank(),
                                   axis.ticks.y=element_blank())
  } else if (plot.position == "right"){
    this.plot <- this.plot + scale_y_discrete(position = "right")+
      theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 10, unit = "pt"))
  }
  return(this.plot)
}

make_df_plot <- function(df_std, top_tbl, selected_receivers) {
  df_std %>%
    dplyr::semi_join(top_tbl, by = "interaction") %>%
    { if (!is.null(selected_receivers)) dplyr::filter(., receiver %in% selected_receivers) else . } %>%
    dplyr::mutate(
      size_center  = if_else(method_score_abs > quantile(method_score_abs, 0.1, na.rm = TRUE), 1, NA_real_),
      stroke_center= 0.5,
      stroke_ligand= if_else(!is.na(ligand.frac)   & ligand.frac   > 0.05, 0.5, NA_real_),
      size_ligand  = if_else(is.na(ligand.frac),    1, ligand.frac),
      stroke_recept= if_else(!is.na(receptor.frac) & receptor.frac > 0.05, 0.5, NA_real_),
      size_recept  = if_else(is.na(receptor.frac),  1, receptor.frac)
    )
}

plotMethodPrioritizedMap <- function(method_name,
                                     df_std,          # output of coerce_*_schema()
                                     top_tbl,         # receiver/interaction table from select_top_per_receiver()
                                     selected_receivers = NULL,  # optional subset of receivers
                                     abs_center_limit = NULL,    # optional symmetric cap for center color
                                     width_cm = 21, height_cm = 11,
                                     out_prefix = NULL,
                                     ligand_col_min = NULL,ligand_col_max = NULL,method_score_min = NULL,method_score_max = NULL,receptor_col_min = NULL,receptor_col_max = NULL) {

  # keep only selected receivers/interactions
  df_plot <- make_df_plot(df_std, top_tbl, selected_receivers)

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

  # choose defaults only when user hasn’t supplied a manual limit
  col_min_lig <- if (is.null(ligand_col_min))    (lim_lig$min-1) else ligand_col_min
  col_max_lig <- if (is.null(ligand_col_max))    (lim_lig$max+1) else ligand_col_max

  col_min_ctr <- if (is.null(method_score_min))  (lim_ctr$min-1) else method_score_min
  col_max_ctr <- if (is.null(method_score_max))  (lim_ctr$max+1) else method_score_max

  col_min_rec <- if (is.null(receptor_col_min))  (lim_rec$min-1) else receptor_col_min
  col_max_rec <- if (is.null(receptor_col_max))  (lim_rec$max+1) else receptor_col_max

  lig_title <- if (method_name == "NicheNet") "scaled pval adapted (lig)" else "-log10(lig pval adj)"
  rec_title <- if (method_name == "NicheNet") "scaled pval adapted (rec)" else "-log10(rec pval adj)"


  # ---- three panels ----
  p_lig <- LocalPlotBubble(
    df = df_plot,
    x_var = "sender",
    color.var = "ligand.diff.expr",
    size.var  = "size_ligand",
    stroke.var= "stroke_ligand",
    plot.position = "left",
    col.min.val = col_min_lig, 
    col.max.val = col_max_lig,
    plot.title = "Ligand",
    x_lab = "SCT", y_lab = "Interaction",
    legend.title = lig_title
  )

  p_ctr <- LocalPlotBubble(
    df = df_plot,
    x_var = "sender",
    color.var = "method_score",
    size.var  = "size_center",
    stroke.var= "stroke_center",
    plot.position = "middle",
    col.min.val = col_min_ctr, 
    col.max.val = col_max_ctr,
    plot.title = paste0(method_name, " score"),
    x_lab = "RCT", y_lab = "",
    legend.title = "method_score"
  )

  p_rec <- LocalPlotBubble(
    df = df_plot,
    x_var = "receiver",
    color.var = "receptor.diff.expr",
    size.var  = "size_recept",
    stroke.var= "stroke_recept",
    plot.position = "middle",
    col.min.val = col_min_rec, 
    col.max.val = col_max_rec,
    plot.title = "Receptor",
    x_lab = "RCT", y_lab = "",
    legend.title = rec_title
  )

  # compose & save (optional)
  composed <- p_lig + p_ctr + p_rec + patchwork::plot_layout(widths = c(2,1,1))

  if (!is.null(out_prefix)) {
    dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
    #png(paste0(out_prefix, ".png"), width = width_cm, height = height_cm, units = "cm", res = 600)
    #print(composed)
    #dev.off()
    write.csv(df_plot, paste0(out_prefix, ".csv"), row.names = FALSE)
  }

  composed
}
##################################
#analysis & parameters
##################################
set.seed(1) 
figures_folder <- "figures_04_08_2025"
output_data_filepath <- "results/cord_pic/data"
dataset_path <- "results/cord_pic"




##################################
# ===== supplementary Figure 1 ==== 
##################################


L_set <- readRDS(file.path(dataset_path,"data/L_set.rds"))

# 1) load & preprocess all three methods for one dataset

plotDecipherPrioritizedMap(dataset_path,top_n=4,dataset_name="supp_figure_1_decipher", abs_decipher_plot_limit = 20,width=21,height=9)
#top_interactions from decipher plot above
top_interactions <- c("SPN-SIGLEC1","IL10-IL10RA","CCL4-CCR1","CD80-CD274","LAMB2-RPSA","IL7-IL7R","ICAM4-ITGB2","LRPAP1-SORL1")
#top_interactions <- c("IFNG-IFNGR1","IFNG-IFNGR2","IL27-IL27RA","NECTIN2-CD96","TNFSF13-TNFRSF14")

decipher_raw    <- readRDS(file.path(dataset_path, "data/decipher_scores_by_cluster.rds"))
nichenet_raw    <- readRDS(file.path(dataset_path, "nichenet/data/prior_table_all_clusters.rds"))
liana_raw       <- read.csv(file.path(dataset_path, "liana/data/liana_p_interaction_results.csv"),
                            row.names=1, check.names=FALSE)

# your existing preprocessors
decipher_df <- preProcessDecipher(decipher_raw)   # from your load_all() script
nichenet_df <- preProcessNicheNet(nichenet_raw)
liana_df    <- preProcessLIANA(liana_raw)

liana_df <- liana_df %>%               
  left_join(
    L_set %>% select(interaction, ligand, receptor),
    by = "interaction"
  ) %>%
  relocate(ligand, receptor, .after = interaction)  

liana_padj <- liana_raw %>% 
  transmute(
    sender   = source,            # same renaming you did in preProcessLIANA
    receiver = target,
    interaction = str_replace(interaction, "\\^", "-"),  # make delimiter match
    ligand_padj,
    receptor_padj
  ) %>% 
  distinct()      # avoid accidental duplicates

liana_df <- liana_df %>% 
  left_join(liana_padj,
            by = c("sender", "receiver", "interaction"))

# DECIPHER 
top_tbl <- tibble::tibble(interaction = top_interactions)

#CD4_T, B, Mono, NK
# focus on one receiver cell at a time; for example “Monocyte”
for(receiver_sel in c("CD4_T","B","Mono","NK")){

    # NICHE NET
    nichenet_std <- coerce_nichenet_schema(nichenet_df)

    p_nn <- plotMethodPrioritizedMap(
        method_name = "NicheNet",
        df_std = nichenet_std,
        top_tbl = top_tbl,
        selected_receivers = receiver_sel,          
        abs_center_limit = NULL,            
        out_prefix = file.path("figures_04_08_2025", paste("supp_fig1_nichenet_",receiver_sel,sep="")),
        method_score_min = 0,
        method_score_max = 1.002,
        ligand_col_min = 0,
        ligand_col_max = 3,
        receptor_col_min = 0,
        receptor_col_max = 3
    )

    # LIANA+
    liana_std <- coerce_liana_schema(liana_df)

    p_li <- plotMethodPrioritizedMap(
        method_name = "LIANA+",
        df_std = liana_std,
        top_tbl = top_tbl,
        selected_receivers = receiver_sel,
        abs_center_limit = NULL,
        out_prefix = file.path("figures_04_08_2025", paste("supp_fig1_liana_",receiver_sel,sep="")),
        method_score_min = -5,
        method_score_max = 5,
        ligand_col_min = 0,
        ligand_col_max = 3,
        receptor_col_min = 0,
        receptor_col_max = 10
    )

    png_filename <- paste("supp_figure_1_",receiver_sel,".png",sep="")
    png(file.path(figures_folder,png_filename),width  = 2000,     
        height = 2200,     
        res    = 300)
    # Arrange three method panels vertically if you want one composite figure:
    print((p_nn / p_li))
    dev.off()
}

#####################################
##### supplementary figure 2 ########
#####################################

# ── 1. read the full interaction-potential matrix ───────────────────────────
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


#attempt at overlap
per_sign_n <- 100                    # ← bump this until the overlap is rich enough
get_top_pool <- function(df, score_col, k = per_sign_n) {
  score_col <- rlang::ensym(score_col)
  df %>%
    mutate(.score = !!score_col,
           .sign  = if_else(.score >= 0, "pos", "neg")) %>%
    group_by(receiver, sender, interaction, .sign) %>%
    summarise(score = mean(.score, na.rm = TRUE), .groups = "drop") %>%
    group_by(receiver, .sign) %>%
    slice_max(order_by = score * if_else(.sign == "pos",  1, -1),
              n = k, with_ties = FALSE) %>%
    ungroup() %>%
    select(receiver, interaction, .sign) %>%
    unique()
}

top_dec   <- get_top_pool(decipher_df,  prioritization_score,  per_sign_n)
top_niche <- get_top_pool(nichenet_df,  prioritization_score, per_sign_n)
top_liana <- get_top_pool(liana_df,     prioritization_score,    per_sign_n)

overlap <- reduce(list(top_dec, top_niche, top_liana),
                  \(x, y) inner_join(x, y,
                     by = c("receiver", "interaction", ".sign")))

# fallback: if overlap still < 8, just keep whatever is there
final_n <- 8
result <- bind_rows(
  overlap %>% filter(.sign == "pos") %>%
    slice_max(order_by = row_number(), n = final_n / 2),
  overlap %>% filter(.sign == "neg") %>%
    slice_max(order_by = row_number(), n = final_n / 2)
)

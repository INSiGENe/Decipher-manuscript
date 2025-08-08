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


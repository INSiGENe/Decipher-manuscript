# ── libraries ────────────────────────────────────────────────────────────────
library(dplyr)
library(stringr)
library(gplots)          # heatmap.2
library(RColorBrewer)
library(matrixStats)     # rowVars()

# ── 1. read the full interaction-potential matrix ───────────────────────────
output_data_filepath <- "results/cord_pic/data"
pot_file <- file.path(output_data_filepath,
                      "interaction_potential_by_clusters.rds")

all_potentials <- readRDS(pot_file)
# `all_potentials` should be a matrix or data.frame
#   rows   = LR pairs   (e.g. "VIM-CD44")
#   cols   = receiver clusters (e.g. "B_cells", "pDC", …)

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
receptor_pal <- brewer.pal(8, "Set2")
receptor_colours <- setNames(
  receptor_pal[seq_along(unique(receptor_vec))],
  unique(receptor_vec)
)

row_side <- receptor_colours[receptor_vec]
col_side <- row_side  # same order, symmetric matrix

# ── 6. dendrogram (optional – distance = 1-ρ) -------------------------------
dend <- as.dendrogram(
  hclust(as.dist(1 - cor_mat), method = "average")
)

png(file.path("figures_04_08_2025","supp_figure_2.png"))
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
  density.info = "none",
  labRow = FALSE, labCol = FALSE,
  cexRow = 0.4,   cexCol = 0.4,
  margins = c(5,5),
  main = "LR interaction-potential correlation – B-cells"
)
dev.off()

# wrap in png()/pdf() if you want a file:
# png("fig_LR_corr_Bcells_cordpic.png", 1800, 1800, res = 300)
# heatmap.2(…)
# dev.off()


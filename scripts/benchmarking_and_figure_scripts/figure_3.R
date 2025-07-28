#Load Manuscript package ----
library(devtools)
document()
load_all()

#Required libraries ----
library(ggplot2)
library(dplyr)
library(gplots)
library(circlize)
library(reshape2)
library(igraph)
library(scales)
library(Seurat)
library(nichenetr)
library(stringr)
library(patchwork)
library(networkD3)

#Parameters ----
##SELECT parameters ----
dataset_path <- "manuscript_analysis/flt3_100_v_1"

pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- file.path("reference_data")
output_data_filepath <- file.path(dataset_path,"data")
output_figures_filepath <- file.path(dataset_path,"figures")

#read data ----
#used to be called regulon_scores now called significant_regulons_by_cluster
significant_regulons_by_cluster <- readRDS(file.path(output_data_filepath,"significant_regulons_by_cluster.rds"))
lr_markers_by_cluster <- readRDS(file.path(output_data_filepath,"lr_markers_by_cluster.rds"))
#also called decipher_scores in some parts of this code
decipher_scores_by_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
#also called decipher_tf and lr_tf in some parts of this code
decipher_scores_by_regulon_and_cluster <- readRDS(file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
#also called test in slamf7_signature.R script
de_markers_by_cluster <- readRDS(file.path(output_data_filepath,"de_markers_by_cluster.rds"))
#called L_set before
L_set <- readRDS(file.path(output_data_filepath,"L_set.rds"))
feature_statistics <- readRDS(file.path(output_data_filepath,"feature_statistics.rds"))
decipher_seurat_lr <- readRDS(file.path(output_data_filepath,"decipher_seurat_lr.rds"))
regulon_deltas_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
regulons_scores_by_clusters <- readRDS(file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
#also called decipher_seurat
pseudobulk_seurat <- readRDS(file.path(output_data_filepath,"pseudobulk_seurat.rds"))
#also called regulon_grns_by_cluster
regulon_grns_by_cluster <- readRDS(file.path(output_data_filepath,"regulon_grns_by_cluster.rds"))
slamf7_high <- read.csv("SLAMF7_high_signature.csv")


#PANEL A - volcano plot ----
## Data Wrangling ----
regulons_scores_by_clusters_c8 <- regulons_scores_by_clusters$Tumour_Cells
condition_match <- match(colnames(regulons_scores_by_clusters_c8),names(decipher_seurat_lr$condition))
group_vector <- decipher_seurat_lr$condition[condition_match]

regulon_deltas_c8 <- regulon_deltas_by_cluster$Tumour_Cells %>%
  filter(class == "real")
regulon_deltas_c8


group_vector[group_vector == "control"] <- 0
group_vector[group_vector == "case"] <- 1
group_factor <- factor(c(group_vector), levels = c(0, 1), labels = c("control", "case"))

diff_regulon_scores_p_values <- do_t_test_by_feature_by_grouping_factor(regulons_scores_by_clusters_c8,group_factor)

regulon_deltas_c8$p_value <- diff_regulon_scores_p_values[regulon_deltas_c8$name]
regulon_deltas_c8$log_10 <-  -1*log(regulon_deltas_c8$p_value,base=10)

#threshold used to be 2.883
#not sure what threshold to pick here, but let's say 5 for this code
regulon_signature <- regulon_deltas_c8 %>%
  filter(log_10 > -log(0.01,base=10) & abs(deltaPagoda) > 2) %>%
  pull(name)

##visualization parameters ----
abs_max_regulon_delta <- max(abs(regulon_deltas_c8$deltaPagoda))

##visualization ----
p <- ggplot(regulon_deltas_c8, aes(x = deltaPagoda, y = log_10)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = 2.883, linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = -2.883, linetype = "dashed", colour = "grey") +
  geom_point(shape = 16, alpha = 0.5) +
  scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
  scale_x_continuous(limits = c(-1*6,abs_max_regulon_delta), breaks = seq(-6,6,2)) +
  scale_y_continuous(limits = c(0,300)) +
  xlab("delta TF activity") +
  ylab("-log10(P)") +
  theme_classic(9) +
  theme(legend.position = "none",
        axis.title.x= element_text(size = 20),
        axis.title.y= element_text(size = 20))

p <- p + geom_text(data = subset(regulon_deltas_c8, log_10 > -log(0.01,base=10) & abs(deltaPagoda) > 2.883), aes(label = name),
                   vjust = "inward", hjust = "inward", check_overlap = TRUE,size=4)

# Print the plot
# png(file.path(output_data_filepath,"volcano_plot_C8.png"),
#     height = 20,
#     width = 12,
#     units = "cm",
#     res=500)
print(p)
# dev.off()

# write.csv(regulon_deltas_c8,file.path("figures","figure_5_panel_a.csv"))

#PANEL B (real) ----
decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)
decipher_scores_by_cluster_bound_filtered <- decipher_scores_by_cluster_bound %>%
  mutate(decipher_score = sign(decipher_score)*log10(abs(decipher_score)+1)) %>%
  filter(abs(decipher_score) > 0.4)

circos_data  <- decipher_scores_by_cluster_bound_filtered %>%
  select(interaction,receiver_cluster,decipher_score) %>%
  dplyr::rename(from=interaction, to=receiver_cluster, value=decipher_score)
circos_data$to <- clean_names(circos_data$to)

# Use the sub function to retain only the ligand from the ligand-receptor pair
circos_data$from <- sub("-.*", "", circos_data$from)
circos_data_matrix <- reshape2::acast(circos_data,from~to,value.var = "value",fill=0,fun.aggregate = sum)
colnames(circos_data_matrix) <- convert_text_patterns(colnames(circos_data_matrix) )
#circos_data_matrix_subset <- circos_data_matrix[,c("pDC","cDC2","CD14+BDCA1+PD-L1+cells","CD14+monocytes","CD16+monocytes")]
#colnames(circos_data_matrix_subset) <- c("pDC","cDC2","C8","CD14+ Mono","CD16+ Mono")
#desired_order <- c("B", "CD4_T", "CD8_T","Naive_CD8_T", "NK","cDC2","C8","CD14+Mono","CD16+Mono") # Example order
#circos_data_matrix <- circos_data_matrix[,c(desired_order)]

df <- melt(circos_data_matrix)
colnames(df) <- c("from", "to", "weight")

# png(file.path("figures","fig5_panel_b_circos_prioritized_ligands_into_all_cts.png"),width=15,height = 15,units = "cm",res=400)
circos.par(start.degree = 90, track.margin = c(0, 0))
chordDiagram(
  x = df,
  transparency = 0.5,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight  = -0.04,
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  preAllocateTracks = 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.9) # Adjust 'cex' for smaller text
}, bg.border = NA)

circos.clear()
# dev.off()

# write.csv(df,file.path("figures","fig5_panel_b_circos_prioritized_ligands_into_all_cts.csv"))


#PANEL C ----
##relies on Panel B ----
ligands_expressed <- rownames(circos_data_matrix)[abs(rowSums(circos_data_matrix))>0]

##data wrangling ----
ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)

#normalize feature statistics
normalized_feature_statistics <- feature_statistics
normalized_feature_statistics$normalized.counts <- normalized_feature_statistics$sum.counts/normalized_feature_statistics$n.cell
#total counts feature condition
normalized_feature_statistics <- normalized_feature_statistics %>%
  group_by(condition,feature) %>%
  mutate(total.normalized.counts = sum(normalized.counts)) %>%
  ungroup() %>%
  mutate(frac.normalized.counts.features.condition = normalized.counts/total.normalized.counts)

bubble_plot_data <- normalized_feature_statistics %>%
  dplyr::left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("cluster"="cluster","feature"="gene")) %>%
  mutate(avg_log2FC = replaceNAw0(avg_log2FC)) %>%
  select(cluster,condition,feature,frac.normalized.counts.features.condition,avg_log2FC) %>%
  filter(feature %in% ligands_expressed & condition == "case")

x_lab <- "CT"
y_lab <- "Ligand"
col.min.val <- -1*max(abs(bubble_plot_data$avg_log2FC))
col.max.val <- max(abs(bubble_plot_data$avg_log2FC))
plot.title <- "test"

bubble_plot_data$cluster <- clean_names(bubble_plot_data$cluster)

# Convert 'feature' to a factor and order it (if not already done)
bubble_plot_data$feature <- factor(bubble_plot_data$feature, levels = sort(unique(bubble_plot_data$feature)))

# Reverse the factor levels so A is at the top
bubble_plot_data$feature <- forcats::fct_rev(bubble_plot_data$feature)

data_for_circos <- bubble_plot_data %>%
  filter(frac.normalized.counts.features.condition > 0.1)

bubble_plot_data_matrix <- reshape2::acast(data_for_circos,cluster~feature,value.var = "frac.normalized.counts.features.condition",fill=0,fun.aggregate = sum)
rownames(bubble_plot_data_matrix)[2]<- c("C8")
rownames(bubble_plot_data_matrix)[3]<- c("CD14+ Mono")
rownames(bubble_plot_data_matrix)[4]<- c("CD16+ Mono")

df <- melt(bubble_plot_data_matrix)
colnames(df) <- c("from", "to", "weight")
sct_to_ligands <- df
##Visualization ----
# png(file.path("figures","fig5_panel_c_circos_prioritized_ligands_from_all_cts.png"),width=15,height = 15,units = "cm",res=400)
circos.par(start.degree = 90, track.margin = c(0, 0))

chordDiagram(
  x = df,
  transparency = 0.5,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight  = -0.04,
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  preAllocateTracks = 1)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.9) # Adjust 'cex' for smaller text
}, bg.border = NA)
circos.clear()
# dev.off()

# write.csv(df,file.path("figures","fig5_panel_c_circos_prioritized_ligands_from_all_cts.csv"))


#PANEL D ----
##relies on prior sections ----
normalized_feature_statistics
ct_lr_markers
#called decipher_scores_by_cluster_df before
decipher_scores_by_cluster_bound
##data wrangling ----
#first we enrich decipher results with the information we will need for downstream visualization
##enrich ----
#used to be called decipher_scores_by_cluster_df_enriched
decipher_scores_by_cluster_bound_enriched <- decipher_scores_by_cluster_bound %>%
  select(interaction,receiver_cluster,decipher_score) %>%
  # Ensure complete combinations of interaction and receiver_cluster.
  tidyr::complete(interaction,receiver_cluster) %>%
  # Replace missing decipher scores with 0 and add sender_cluster column with value "mixed".
  mutate(decipher_score =tidyr::replace_na(decipher_score,0),
         sender_cluster = "mixed") %>%
  # Left join with ligand-receptor interaction data.
  left_join(select(L_set,ligand,receptor,interaction),by = "interaction") %>%
  # Add differential expression data for ligands.
  left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","ligand"="gene")) %>%
  rename(ligand.diff.expr = avg_log2FC) %>%
  # Add differential expression data for receptors.
  left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","receptor"="gene")) %>%
  rename(receptor.diff.expr = avg_log2FC) %>%
  # Replace missing values in differential expression columns with 0 and add condition column with value "case".
  mutate(ligand.diff.expr = replaceNAw0(ligand.diff.expr),
         receptor.diff.expr = replaceNAw0(receptor.diff.expr),
         condition = "case") %>%
  # Left join with normalized feature statistics for ligands.
  left_join(select(normalized_feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","ligand"="feature","condition"))%>%
  rename(ligand.frac = frac.normalized.counts.features.condition) %>%
  # Left join with normalized feature statistics for receptors.
  left_join(select(normalized_feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","receptor"="feature","condition"))%>%
  rename(receptor.frac = frac.normalized.counts.features.condition)

## top interactions ----
#used to be called decipher_for_overlap/merged_data and was comprised of decipher_bound and decipher_pre_processed
decipher_scores_by_cluster_bound_clean <- decipher_scores_by_cluster_bound %>%
  mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
  rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score) %>%
  select(interaction,ligand,receptor,receiver,prioritization_score) %>%
  arrange(prioritization_score)


#Merge the results from the three methods
# selected_cts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
# selected_rcts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
# selected_scts <- c("B_cell","Monocyte","HSC","CD4_T","CD8_Tem","NK_cell_1")

#now
#used to be called decipher_top_interactions_all_rcts
decipher_top_interactions_cluster_c8 <- decipher_scores_by_cluster_bound_clean %>%
  #filter(receiver %in% "Tumour_Cells") %>%
  mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
  #group_by(decipher_score_sign) %>%
  arrange(desc(abs(prioritization_score))) %>%
  select(interaction) %>%
  distinct() %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  left_join(decipher_scores_by_cluster_bound)

top_interactions <- decipher_top_interactions_cluster_c8 %>%
  select(interaction) %>%
  distinct() %>%
  unlist(use.names=FALSE)

#which(top_interactions == "IFNG-IFNGR2")

## visualization ----
base_data <- decipher_scores_by_cluster_bound_enriched %>%
  filter(interaction %in% top_interactions) %>%
  mutate(size = 1)

plot_limits_ligand <- list(max = max(base_data$ligand.diff.expr),min = min(base_data$ligand.diff.expr))
plot_limits_receptor <- list(max = max(base_data$receptor.diff.expr),min = min(base_data$receptor.diff.expr))
plot_limits_decipher <- list(max = max(base_data$decipher_score),min = min(base_data$decipher_score))

base_data$stroke <- 0.5

base_data <- base_data %>%
  mutate(stroke_ligand = if_else(ligand.frac > 0.05,0.5,NA)) %>%
  mutate(size_ligand = if_else(ligand.frac > 0.05,ligand.frac,NA)) %>%
  mutate(stroke_receptor = if_else(receptor.frac > 0.05,0.5,NA),
         size_receptor = if_else(receptor.frac > 0.05,receptor.frac,NA)) %>%
  mutate(size = if_else(abs(decipher_score) > 0.1,1,NA))%>%
  mutate(receiver_cluster=if_else(receiver_cluster == "Tumour_Cells","C8",receiver_cluster))


base_data$receiver_cluster <- convert_text_patterns(base_data$receiver_cluster)

ligand_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster" ,
  color.var = "ligand.diff.expr",
  size.var = "size_ligand",
  stroke.var = "stroke_ligand",
  plot.position = "left",
  col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
  plot.title = "Ligand",
  x_lab= "SCT",
  y_lab = "Interaction")

decipher_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster",
  color.var = "decipher_score",
  size.var = "size",
  stroke.var = "stroke",
  plot.position = "middle",
  col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
  plot.title = "Decipher score",
  x_lab= "RCT",
  y_lab = "")

receptor_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster",
  color.var = "receptor.diff.expr",
  size.var = "size_receptor",
  stroke.var = "stroke_receptor",
  plot.position = "middle",
  col.min.val=plot_limits_receptor$min,col.max.val=plot_limits_receptor$max,
  plot.title = "Receptor",
  x_lab= "RCT",
  y_lab = "")

#png(file.path("figures","fig5_panel_d_decipher_plot_prioritized_by_c8.csv.png"),width = 21,height = 11,units = "cm",res = 600)
ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot
#dev.off()

#write.csv(base_data,file.path("figures","fig5_panel_d_decipher_plot_prioritized_by_c8.csv"))

#PANEL E ----
##relies on prior sections ----
decipher_scores_by_regulon_and_cluster
regulon_signature
top_interactions
##data wrangling ----
to_plot <- decipher_scores_by_regulon_and_cluster$Tumour_Cells %>%
  filter(regulon %in% regulon_signature,
         interaction %in% top_interactions)


this_matrix <-reshape2::acast(to_plot,interaction~regulon,value.var = "decipher_score",fill = 0)

##visualization ----
#png(file.path("figures","fig5_panel_e_LR-TF connection.png"),width = 15, height = 9, units = "cm",res=600)
heatmap.2(
  this_matrix,
  trace="none",
  col = "bluered",
  breaks = 100,
  cexRow = 0.7,
  cexCol=0.7,
  scale = "none",
  key.title = "LR-TF Decipher score",Rowv = FALSE,
  dendrogram = "none",
  margins = c(6,6),
  key = FALSE,
  keysize = 0.3,
  Colv=TRUE)
#dev.off()
write.csv(this_matrix,file.path("figures","fig5_panel_e_LR-TF connection.csv"))

#PANEL F ----
## data wrangling ----
slamf7_high <- slamf7_high[,c("Gene","log2FoldChange.SF","log2FoldChange.PB")]


main_signature_synovial_fluid <- slamf7_high %>%
  dplyr::rename(avg_log2FC=log2FoldChange.SF,gene=Gene) %>%
  mutate(cell_type = "ref_syn_fluid")

main_signature_pbmc <- slamf7_high %>%
  dplyr::rename(avg_log2FC=log2FoldChange.PB,gene=Gene) %>%
  mutate(cell_type = "ref_pbmc")

cd14_monocytes_filtered <- de_markers_by_cluster$CD14_plus_monocytes %>%
  filter(gene %in% slamf7_high$Gene) %>%
  mutate(cell_type = 'CD14_plus_monocytes')

cd16_monocytes_filtered <- de_markers_by_cluster$CD16_plus_monocytes %>%
  filter(gene %in% slamf7_high$Gene)%>%
  mutate(cell_type = 'CD16_plus_monocytes')

c8_filtered <- de_markers_by_cluster$Tumour_Cells %>%
  filter(gene %in% slamf7_high$Gene)%>%
  mutate(cell_type = 'C8')

NK_filtered <- de_markers_by_cluster$NK %>%
  filter(gene %in% slamf7_high$Gene)%>%
  mutate(cell_type = 'NK')

CD8_filtered <- de_markers_by_cluster$CD8_T %>%
  filter(gene %in% slamf7_high$Gene)%>%
  mutate(cell_type = 'CD8 T')

pDC_filtered <- de_markers_by_cluster$pDC %>%
  filter(gene %in% slamf7_high$Gene)%>%
  mutate(cell_type = 'pDC')

combined <- bind_rows(pDC_filtered,NK_filtered,CD8_filtered,main_signature_synovial_fluid,main_signature_pbmc,cd14_monocytes_filtered, cd16_monocytes_filtered, c8_filtered)
new_matrix <- reshape2::acast(combined,gene~cell_type,value.var = "avg_log2FC")
colnames(new_matrix) <- convert_text_patterns(colnames(new_matrix))

new_matrix <- new_matrix[,c(8,7,1,2,3,4,5,6)]

##visualization ----
# Define the color palette from light red to dark red
red_palette <- colorRampPalette(c("lightcoral","red", "darkred"))(256)
divergent_palette <- colorRampPalette(c("lightcoral", "darkred", "purple"))(22)

png(file.path("figures","fig5_panel_f_slamf7_signature.png"),height = 12,width=10,units="cm",res=400)
heatmap.2(
  new_matrix,
  trace = "none",
  col = divergent_palette,
  density.info = "none",
  scale = "none",
  cexRow = 0.7,
  cexCol = 0.7,
  margin = c(5,5),
  key.title = "log2FC",
  dendrogram = "none",
  Rowv = FALSE,
  Colv= FALSE,
  keysize = 2)
dev.off()
write.csv(new_matrix,file.path("figures","fig5_panel_f_slamf7_signature.csv"))

#Supplementary


#other code ----


ligands_to_rct <- df
colnames(circos_data_matrix)
circos_data_matrix_subset <- circos_data_matrix[,c("B","CD4_T","CD8_T","Naive_CD8_T","NK")]

png(file.path("figures","circos_adaptive_cells.png"),width=15,height = 15,units = "cm",res=400)
circos.par(start.degree = 90, track.margin = c(0, 0))
chordDiagram(circos_data_matrix_subset,transparency = 0.5,
             annotationTrack = "grid",
             preAllocateTracks = 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6) # Adjust 'cex' for smaller text
}, bg.border = NA)
circos.clear()
dev.off()

decipher_scores_matrix_for_heatmap <- reshape2::acast(decipher_scores_by_cluster_bound_filtered,interaction~receiver_cluster,value.var = "decipher_score",fill=0)

png("all_immune_cells.png",width = 6,height = 15,units="cm",res=600)
gplots::heatmap.2(
  decipher_scores_matrix_for_heatmap,
  trace="none",
  cexRow=0.3,
  cexCol=0.3,
  col="bluered",
  breaks = 20,
  dendrogram = "row",
  key = FALSE,
  keysize = 0.3,
  margins = c(5,5))
dev.off()






#other code ----
df <- regulon_scores$Tumour_Cells
df <- df %>%
  filter(name %in% c("STAT2","IRF7","IRF2","IRF8","STAT1","FOSL2","ATF3","JUND","JUNB","FOS","FOSB","IRF1","STAT3"))
df <- df[order(df$deltaPagoda), ]

# Create the plot
p <- ggplot(df, aes(x = class, y = reorder(name, deltaPagoda), fill = deltaPagoda)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white",
                       midpoint = 0, limit = c(min(df$deltaPagoda), max(df$deltaPagoda)),
                       name = expression(Delta * " TF")) +
  theme_minimal() +
  labs(y = "", x = "C8") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

p <- ggplot(df, aes(x = class, y = reorder(name, deltaPagoda), fill = deltaPagoda)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white",
                       midpoint = 0, limit = c(min(df$deltaPagoda), max(df$deltaPagoda)),
                       name = expression(Delta * " TF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size=5),
        axis.text.y = element_blank())+
  labs(y = "", x = "C8")+
  theme(
    legend.text = element_text(size = 8),          # Smaller text in the legend
    legend.title = element_text(size = 7),         # Smaller legend title
    legend.key.size = unit(0.4, "lines"),          # Smaller keys in the legend
    legend.spacing = unit(0.5, "lines"),           # Tighter spacing
    legend.margin = margin(5, 5, 5, 5),             # Adjust margins around the legend
    legend.position = "right"
  )+
  coord_flip()

print(p)

regulon_deltas_by_cluster$Tumour_Cells

#ok so let's do this as a p-value and compare to other cell types


#jeez I need to pick a value here
#used to be 5
selected_interactions <- decipher_scores_by_regulon_and_cluster$Tumour_Cells %>%
  filter(regulon %in% regulon_signature,
         abs(decipher_score) > 12) %>%
  pull(interaction) %>%
  unique()

#top_interactions used to be selected_interactions, top_interactions are selected downstream so there is a break here
#another one names(table(regulon_grns_by_cluster_filtered$source)[ind]), which used to be regulon_signature
top_interactions <- selected_interactions
to_plot <- decipher_scores_by_regulon_and_cluster$Tumour_Cells %>%
  filter(regulon %in% names(table(regulon_grns_by_cluster_filtered$source)[ind]),
         interaction %in% top_interactions)

to_plot <- decipher_scores_by_regulon_and_cluster$Tumour_Cells %>%
  filter(regulon %in% regulon_signature,
         interaction %in% top_interactions)


this_matrix <-reshape2::acast(to_plot,interaction~regulon,value.var = "decipher_score",fill = 0)
png("LR-TF connection.png",width = 15, height = 9, units = "cm",res=600)
heatmap.2(
  this_matrix,
  trace="none",
  col = "bluered",
  breaks = 100,
  cexRow = 0.7,
  cexCol=0.7,
  scale = "none",
  key.title = "LR-TF Decipher score",Rowv = FALSE,
  dendrogram = "none",
  margins = c(6,6),
  key = FALSE,
  keysize = 0.3,
  Colv=TRUE)
dev.off()
to_plot_summary <- to_plot %>%
  filter(abs(decipher_score) > 1) %>%
  group_by(interaction) %>%
  summarize(weight = sum(decipher_score)) %>%
  ungroup()

print(to_plot_summary%>%
        arrange(weight),n=54)





 # Print the plot
png(file.path(output_data_filepath,"regulon_delta_findings.png"),height = 10,width = 5
    , units = "cm",res=500)
print(p)
dev.off()


case_cluster <- "Tumour_Cells"
control_cluster <- "CD14_plus_monocytes"

cells <- pseudobulk_seurat@meta.data %>%
  filter((cluster == case_cluster & condition == "case") | (cluster == control_cluster & condition == "case")) %>%
  pull(cell)

pseudobulk_seurat_case_cluster <- subset(pseudobulk_seurat,cells = cells)
pseudobulk_seurat_case_cluster@meta.data <- pseudobulk_seurat_case_cluster@meta.data %>%
  mutate(condition = if_else(cluster == control_cluster,"control","case"))

flag.normalize.non.log <- FALSE
param_max_n_cells <- 600
param_min_n_cells <- 100

if(flag.normalize.non.log){
  pseudobulk_seurat_case_cluster <- NormalizeData(pseudobulk_seurat_case_cluster,normalization.method = "RC",scale.factor=100000)
}

base_n_cells <- min(table(pseudobulk_seurat_case_cluster$condition))
if(base_n_cells < param_min_n_cells){
  next
}
if(base_n_cells > param_max_n_cells){
  base_n_cells <- param_max_n_cells
}

Idents(pseudobulk_seurat_case_cluster) <- pseudobulk_seurat_case_cluster@meta.data$condition
pseudobulk_seurat_case_cluster_downsampled <- subset(pseudobulk_seurat_case_cluster, downsample = base_n_cells)



data_case_cluster_downsampled <- pseudobulk_seurat_case_cluster_downsampled@assays$RNA@data


#dplyr doesn't work with sparse matrix
#data_case_cluster_downsampled_receptors <- data_case_cluster_downsampled[which(rownames(data_case_cluster_downsampled) %in% selected_receptors),]
data_case_cluster_downsampled_receptors <- data_case_cluster_downsampled[which(rownames(data_case_cluster_downsampled) %in% unique(L_set_relevant_features$receptor)),]
##PAGODA -----
#silence this function
regulon_scores_case_cluster <- getRegulonScores(
  seuratObject = pseudobulk_seurat_case_cluster,
  grn_df = regulon_case_cluster_capped)

##PAGODA DELTA ----
regulon_deltas_case_cluster <- getRegulonDeltas(
  regulon_scores_case_cluster,
  pseudobulk_seurat_case_cluster$condition)

regulon_deltas_case_cluster <- regulon_deltas_case_cluster %>%
  mutate(class = ifelse(stringr::str_detect(name,"sample"),"random","real"))

random.density <- density(subset(regulon_deltas_case_cluster, class == "random")$deltaPagoda, n = 2^10)
upper_threshold_random <- quantile(random.density,probs = 0.975,normalize = FALSE)
lower_threshold_random <- quantile(random.density,probs = 0.025,normalize = FALSE)

significant_regulon_deltas_case_cluster <-  regulon_deltas_case_cluster %>%
  filter(class == "real") %>%
  filter(deltaPagoda > upper_threshold_random | deltaPagoda < lower_threshold_random)%>%
  arrange(deltaPagoda)


# Assuming your data frame is named 'df'
# Replace 'df' with the actual name of your data frame

# Your data frame, sorted by 'deltaPagoda'
df <- significant_regulon_deltas_case_cluster
df <- df %>%
  filter(name %in% c("STAT2","IRF7","IRF2","IRF8","STAT1","FOSL2","ATF3","JUND","JUNB","FOS","FOSB","IRF1","STAT3"))
df <- df[order(df$deltaPagoda), ]

# Create the plot
p <- ggplot(df, aes(x = class, y = reorder(name, deltaPagoda), fill = deltaPagoda)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white",
                       midpoint = 0, limit = c(min(df$deltaPagoda), max(df$deltaPagoda)),
                       name = expression(Delta * " TF")) +
  theme_minimal() +
  labs(y = "", x = "C8") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

# Print the plot
png(file.path(output_data_filepath,"regulon_delta_findings.png"),height = 10,width = 5
    , units = "cm",res=500)
print(p)
dev.off()

#heatmap of ligands ----
ligands_expressed <- rownames(circos_data_matrix_subset)[abs(rowSums(circos_data_matrix_subset))>0]
#ligands_expressed <- rownames(circos_data_matrix_subset)

#ok time to normalize feature statisctics
feature_statistics$normalized.counts <- feature_statistics$sum.counts/feature_statistics$n.cell
#total counts feature condition
feature_statistics <- feature_statistics %>%
  group_by(condition,feature) %>%
  mutate(total.normalized.counts = sum(normalized.counts)) %>%
  ungroup() %>%
  mutate(frac.normalized.counts.features.condition = normalized.counts/total.normalized.counts)

feature_statistics %>%
  select(cluster,condition,feature,sum.counts,n.cell,normalized.counts,total.normalized.counts,frac.normalized.counts.features.condition)

plot(feature_statistics$sum.counts,feature_statistics$normalized.counts)
plot(x = feature_statistics$sum.counts,y = feature_statistics$frac.normalized.counts.features.condition)
plot(x = feature_statistics$frac.counts.features.condition,
     y = feature_statistics$frac.normalized.counts.features.condition,
     xlab = "before",
     ylab = "new")
dev.off()

to_plot <- feature_statistics %>% filter(feature %in% ligands_expressed & condition == "case") %>%
  select(cluster,feature,frac.normalized.counts.features.condition) %>%
  mutate(frac.normalized.counts.features.condition = if_else(is.na(frac.normalized.counts.features.condition), 0, frac.normalized.counts.features.condition))


to_plot_heatmap <- acast(to_plot,cluster~feature,value.var = "frac.normalized.counts.features.condition")
png("mechanistic_insight_ligand_expression.png",height = 15,width=15,units="cm",res=400)
heatmap.2(to_plot_heatmap,trace="none",cexRow = 0.6,cexCol=0.6,margins = c(7,7))
dev.off()


ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)

bubble_plot_data <- feature_statistics %>%
  dplyr::left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("cluster"="cluster","feature"="gene")) %>%
  mutate(avg_log2FC = replaceNAw0(avg_log2FC)) %>%
  select(cluster,condition,feature,frac.normalized.counts.features.condition,avg_log2FC) %>%
  filter(feature %in% ligands_expressed & condition == "case")

x_lab <- "CT"
y_lab <- "Ligand"
col.min.val <- -1*max(abs(bubble_plot_data$avg_log2FC))
col.max.val <- max(abs(bubble_plot_data$avg_log2FC))
plot.title <- "test"

bubble_plot_data$cluster <- clean_names(bubble_plot_data$cluster)

# Convert 'feature' to a factor and order it (if not already done)
bubble_plot_data$feature <- factor(bubble_plot_data$feature, levels = sort(unique(bubble_plot_data$feature)))

# Reverse the factor levels so A is at the top
bubble_plot_data$feature <- forcats::fct_rev(bubble_plot_data$feature)

data_for_circos <- bubble_plot_data %>%
  filter(frac.normalized.counts.features.condition > 0.1)

bubble_plot_data_matrix <- reshape2::acast(data_for_circos,cluster~feature,value.var = "frac.normalized.counts.features.condition",fill=0,fun.aggregate = sum)
rownames(bubble_plot_data_matrix)[2]<- c("C8")
rownames(bubble_plot_data_matrix)[3]<- c("CD14+ Mono")
rownames(bubble_plot_data_matrix)[4]<- c("CD16+ Mono")



df <- melt(bubble_plot_data_matrix)
colnames(df) <- c("from", "to", "weight")
ligands_to_rct
sct_to_ligands <- df
#CIRCOS plot
png(file.path("figures","circos_innate_cells_SCT.png"),width=15,height = 15,units = "cm",res=400)
circos.par(start.degree = 90, track.margin = c(0, 0))
# chordDiagram(bubble_plot_data_matrix,transparency = 0.5,
#              annotationTrack = "grid",
#              preAllocateTracks = 1)
chordDiagram(
  x = df,
  transparency = 0.5,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight  = -0.04,
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  preAllocateTracks = 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.9) # Adjust 'cex' for smaller text
}, bg.border = NA)
circos.clear()
dev.off()


this.plot <- ggplot(bubble_plot_data,aes(y=feature, x=cluster,fill = avg_log2FC)) +
  geom_point(aes(size = frac.normalized.counts.features.condition), shape = 21) +
  labs(x = x_lab, y = y_lab) +
  theme_bw()+
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    limits = c(col.min.val,col.max.val)
  )+ggtitle(label = "")+ guides(size = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        plot.margin = ggplot2::margin(t = 10, r = 2, b = 0, l = 2, unit = "pt"),
        plot.title = element_text(hjust = 0.5))

png("test.png",height = 20,width = 12, units="cm",res=400)
print(this.plot)
dev.off()

ligand_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster" ,
  color.var = "ligand.diff.expr",
  size.var = "size_ligand",
  stroke.var = "stroke_ligand",
  plot.position = "left",
  col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
  plot.title = "Ligand",
  x_lab= "SCT",
  y_lab = "Interaction")


base_data <- base_data %>%
  mutate(stroke_ligand = if_else(ligand.frac > 0.05,0.5,NA)) %>%
  mutate(size_ligand = if_else(ligand.frac > 0.05,ligand.frac,NA)) %>%
  mutate(stroke_receptor = if_else(receptor.frac > 0.05,0.5,NA),
         size_receptor = if_else(receptor.frac > 0.05,receptor.frac,NA)) %>%
  mutate(size = if_else(abs(decipher_score) > 1,1,NA))%>%
  mutate(receiver_cluster=if_else(receiver_cluster == "Tumour_Cells","C8",receiver_cluster))


ligand_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster" ,
  color.var = "ligand.diff.expr",
  size.var = "size_ligand",
  stroke.var = "stroke_ligand",
  plot.position = "left",
  col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
  plot.title = "Ligand",
  x_lab= "SCT",
  y_lab = "Interaction")

decipher_scores_bound <- bind_rows(decipher_scores_by_cluster)
decipher_scores_bound <- decipher_scores_bound %>%
  mutate(decipher_score = sign(decipher_score)*log10(abs(decipher_score)+1)) %>%
  filter(receiver_cluster %in% c("CD14_plus_monocytes","CD16_plus_monocytes","cDC2","pDC","Tumour_Cells"),
         abs(decipher_score) > 0.3)

decipher_scores_matrix_for_heatmap <- reshape2::acast(decipher_scores_bound,interaction~receiver_cluster,value.var = "decipher_score",fill=0)

png("innate_immune_cells.png",width = 5,height = 12,units="cm",res=600)
gplots::heatmap.2(
  decipher_scores_matrix_for_heatmap,
  trace="none",
  cexRow=0.3,
  cexCol=0.3,
  col="bluered",
  breaks = 20,
  dendrogram = "row",
  key = FALSE,
  keysize = 3,
  margins = c(5,5))
dev.off()

#LR-TF plot
#so which ligand receptors are associated to TFs of Interest


#network
regulon_subset_of_interest <- c("STAT2","IRF7","IRF2","IRF8","STAT1","FOSL2","ATF3","JUND","JUNB","FOS","FOSB","IRF1","STAT3")
lr_subset_of_interest <- c("CCL2-CCR1","C1QB-CD33","SPN-SIGLEC1","CCL2-CCR2")


df <- regulon_scores$Tumour_Cells
df <- df %>%
  filter(name %in% regulon_subset_of_interest)
df <- df[order(df$deltaPagoda), ]


decipher_scores_by_regulon_and_cluster_this_cluster <- decipher_scores_by_regulon_and_cluster$Tumour_Cells
data_for_network <- decipher_scores_by_regulon_and_cluster_this_cluster %>% filter(
  regulon %in%regulon_subset_of_interest
)

#here i want full positive, full negative and mixed signals, so to do so, I'm going to set a random threshold
hist(decipher_scores_by_regulon_and_cluster_this_cluster$imp.perm,breaks = 100)
#pick top 100 edges

data_for_network <- decipher_scores_by_regulon_and_cluster_this_cluster %>%
  filter(regulon %in%regulon_subset_of_interest) %>%
  slice_max(order_by = imp.perm,n=30)

# data_for_network <- decipher_scores_by_regulon_and_cluster_this_cluster %>% filter(
#   regulon %in%regulon_subset_of_interest,
#   interaction %in% lr_subset_of_interest
# )


data_for_network <- data_for_network %>%
  mutate(edge_width = sign(spearman.cor)*imp.perm/2,
         edge_color = if_else(edge_width > 0, "red","blue"))

# Assuming data_for_network and spearman.cor, imp.perm are already defined
data_for_network <- data_for_network %>%
  mutate(
    edge_width = sign(spearman.cor) * imp.perm / 2
  )

# Find the range of edge_width
min_width <- min(data_for_network$edge_width)
max_width <- max(data_for_network$edge_width)
max_width_abs <- max(c(abs(min_width),abs(max_width)))
# Create 10 intervals

breaks <- seq(-1*max_width_abs, max_width_abs, length.out = 11)

# Use cut to create a factor variable
data_for_network <- data_for_network %>%
  mutate(
    edge_color_factor = cut(edge_width, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  )

# Map factor levels to colors (gradient from blue to red)
color_palette <- colorRampPalette(c("blue","white", "red"))(10)
data_for_network$edge_color <- color_palette[data_for_network$edge_color_factor]

# View the result
head(data_for_network)


this_graph <- igraph::graph_from_data_frame(data_for_network[,c("interaction","regulon","edge_width")])

type_ind <- if_else(vertex_attr(this_graph)$name %in% data_for_network$interaction,0,1)

this_graph <- set_vertex_attr(this_graph,"type", value = type_ind)
#this_graph <- set_vertex_attr(this_graph,"vertex_value", value = c(7.82,7.16,3.41,3.40,3.55,5.76,6.04,6.41,6.52,8.36,9.34,10.35,10.37))

V(this_graph)$name <- gsub("-", "\n", V(this_graph)$name)

png("signalling_network.png",width = 20,height = 20,units="cm",res=400)
plot(this_graph,
     layout = layout_as_bipartite(this_graph,hgap = 1,vgap = 3),
     edge.width = abs(data_for_network$edge_width),
     edge.color = data_for_network$edge_color,
     vertex.label.cex = 0.5,
     vertex.label.font = 2,
     vertex.frame.color = NA,
     vertex.color = "gray",
     vertex.size = 18,
     vertex.stroke = NA,
     edge.arrow.size = 0.3,
     margin = c(0,0,0,0))
dev.off()

list.files(output_data_filepath)

regulon_grns_by_cluster_C8 <- regulon_grns_by_cluster$Tumour_Cells
regulon_grns_by_cluster_C8_selected <- regulon_grns_by_cluster_C8 %>% filter(source %in% regulon_subset_of_interest)
interactions_to_consider <- unique(data_for_network$interaction)
genes_to_consider <- unique(unlist(strsplit(interactions_to_consider, split = "-")))

feedback_loops <- regulon_grns_by_cluster_C8_selected %>%
  filter(target %in% genes_to_consider)

#ok now remove redundancy
data_for_network <- data_for_network %>%
  mutate(edge_width = sign(spearman.cor)*imp.perm/2,
         edge_color = if_else(edge_width > 0, "red","blue")) %>%
  filter(!(regulon %in% feedback_loops$source & (ligand %in% feedback_loops$target | receptor %in% feedback_loops$target)))

#CCC map ----

#ok time to normalize feature statisctics
feature_statistics$normalized.counts <- feature_statistics$sum.counts/feature_statistics$n.cell
#total counts feature condition
feature_statistics <- feature_statistics %>%
  group_by(condition,feature) %>%
  mutate(total.normalized.counts = sum(normalized.counts)) %>%
  ungroup() %>%
  mutate(frac.normalized.counts.features.condition = normalized.counts/total.normalized.counts)

#data pre-processing ----

#first we enrich decipher results with the information we will need for downstream visualization
#Here we enrich the Decipher plots with the information we require
ct_lr_markers <- getLigandReceptorDiffExprMarkersByCt(lr_markers_by_cluster)
decipher_scores_by_cluster_df <- bind_rows(decipher_scores_by_cluster)

### fill gaps and enrich
decipher_scores_by_cluster_df_enriched <- decipher_scores_by_cluster_df %>%
  select(interaction,receiver_cluster,decipher_score) %>%
  tidyr::complete(interaction,receiver_cluster) %>%
  mutate(decipher_score =tidyr::replace_na(decipher_score,0),
         sender_cluster = "mixed") %>%
  left_join(select(L_set,ligand,receptor,interaction),by = "interaction") %>%
  #add differential expression data
  left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","ligand"="gene")) %>%
  rename(ligand.diff.expr = avg_log2FC) %>%
  left_join(ct_lr_markers[,c("cluster","gene","avg_log2FC")],by = c("receiver_cluster"="cluster","receptor"="gene")) %>%
  rename(receptor.diff.expr = avg_log2FC) %>%
  mutate(ligand.diff.expr = replaceNAw0(ligand.diff.expr),
         receptor.diff.expr = replaceNAw0(receptor.diff.expr),
         condition = "case") %>%
  left_join(select(feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","ligand"="feature","condition"))%>%
  rename(ligand.frac = frac.normalized.counts.features.condition) %>%
  left_join(select(feature_statistics,cluster,feature,condition,frac.normalized.counts.features.condition), by = c("receiver_cluster"="cluster","receptor"="feature","condition"))%>%
  rename(receptor.frac = frac.normalized.counts.features.condition)

#this section is designed to prepare the results from the three methods for downstream visualization


decipher_bound <- bind_rows(decipher_scores_by_cluster)
decipher_pre_processed <- decipher_bound %>%
  mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
  rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score)



#here we calculate overlap, as for comparison purposes we require a matching set of interactions
decipher_for_overlap <- decipher_pre_processed %>%
  select(interaction,ligand,receptor,receiver,prioritization_score) %>%
  arrange(prioritization_score)

#Merge the results from the three methods
merged_data <- decipher_for_overlap
##pick the top n interactions and get all values associated to those interactions
#$merged_data %>% select(receiver) %>% distinct()
selected_cts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
selected_rcts <- c("B_cell","Monocyte","CD4_T","NK_cell_1")
selected_scts <- c("B_cell","Monocyte","HSC","CD4_T","CD8_Tem","NK_cell_1")

#now
decipher_top_interactions_all_rcts <- merged_data %>%
  filter(receiver %in% "Tumour_Cells") %>%
  mutate(decipher_score_sign=if_else(prioritization_score >=0,"positive","negative")) %>%
  #group_by(decipher_score_sign) %>%
  arrange(desc(abs(prioritization_score))) %>%
  select(interaction) %>%
  distinct() %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  left_join(decipher_scores_by_cluster_df)

top_interactions <- decipher_top_interactions_all_rcts %>%
  select(interaction) %>%
  distinct() %>%
  unlist(use.names=FALSE)

which(top_interactions == "IFNG-IFNGR2")
#top_interactions <- c("SERPING1-LRP1","FGF2-CD44","CCL2-CCR1","C1QB-CD33")

#main-analysis -----
##DECIPHER ----
### full data ----
base_data <- decipher_scores_by_cluster_df_enriched %>%
  filter(interaction %in% top_interactions) %>%
  mutate(size = 1)

plot_limits_ligand <- list(max = max(base_data$ligand.diff.expr),min = min(base_data$ligand.diff.expr))
plot_limits_receptor <- list(max = max(base_data$receptor.diff.expr),min = min(base_data$receptor.diff.expr))
plot_limits_decipher <- list(max = max(base_data$decipher_score),min = min(base_data$decipher_score))

base_data$stroke <- 0.5

base_data <- base_data %>%
  mutate(stroke_ligand = if_else(ligand.frac > 0.05,0.5,NA)) %>%
  mutate(size_ligand = if_else(ligand.frac > 0.05,ligand.frac,NA)) %>%
  mutate(stroke_receptor = if_else(receptor.frac > 0.05,0.5,NA),
         size_receptor = if_else(receptor.frac > 0.05,receptor.frac,NA)) %>%
  mutate(size = if_else(abs(decipher_score) > 0.1,1,NA))%>%
  mutate(receiver_cluster=if_else(receiver_cluster == "Tumour_Cells","C8",receiver_cluster))


base_data$receiver_cluster <- convert_text_patterns(base_data$receiver_cluster)

ligand_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster" ,
  color.var = "ligand.diff.expr",
  size.var = "size_ligand",
  stroke.var = "stroke_ligand",
  plot.position = "left",
  col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
  plot.title = "Ligand",
  x_lab= "SCT",
  y_lab = "Interaction")

decipher_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster",
  color.var = "decipher_score",
  size.var = "size",
  stroke.var = "stroke",
  plot.position = "middle",
  col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
  plot.title = "Decipher score",
  x_lab= "RCT",
  y_lab = "")

receptor_bubble_plot <- plotBubble(
  df = base_data,
  x_var = "receiver_cluster",
  color.var = "receptor.diff.expr",
  size.var = "size_receptor",
  stroke.var = "stroke_receptor",
  plot.position = "middle",
  col.min.val=plot_limits_receptor$min,col.max.val=plot_limits_receptor$max,
  plot.title = "Receptor",
  x_lab= "RCT",
  y_lab = "")

png("pfizer_mrna_v4_ccc_map_v3.png",width = 21,height = 11,units = "cm",res = 600)
ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot
dev.off()

### subset data ----
ligand_bubble_plot <- plotBubble(
  df = base_data %>% filter(receiver_cluster %in% selected_scts),
  x_var = "receiver_cluster" ,
  color.var = "ligand.diff.expr",
  size.var = "size_ligand",
  stroke.var = "stroke_ligand",
  plot.position = "left",
  col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
  plot.title = "Ligand",
  x_lab= "SCT",
  y_lab = "Interaction")

decipher_bubble_plot <- plotBubble(
  df = base_data %>% filter(receiver_cluster %in% selected_rcts),
  x_var = "receiver_cluster",
  color.var = "decipher_score",
  size.var = "size",
  stroke.var = "stroke",
  plot.position = "middle",
  col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
  plot.title = "Decipher score",
  x_lab= "RCT",
  y_lab = "")

receptor_bubble_plot <- plotBubble(
  df = base_data %>% filter(receiver_cluster %in% selected_rcts),
  x_var = "receiver_cluster",
  color.var = "receptor.diff.expr",
  size.var = "size_receptor",
  stroke.var = "stroke_receptor",
  plot.position = "middle",
  col.min.val=plot_limits_receptor$min,col.max.val=plot_limits_receptor$max,
  plot.title = "Receptor",
  x_lab= "RCT",
  y_lab = "")

png("decipher_ccc_map_subset.png",width = 24,height = 10,units = "cm",res = 600)
ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot
dev.off()

nichenet_plots_full <- list()
nichenet_plots_subset <- list()
isFirst <- TRUE
isFirstSubset <- TRUE
for(this_rct in names(nichenet_prior_table_all_clusters)){
  base_data <- nichenet_prior_table_all_clusters[[this_rct]] %>%
    mutate(interaction = paste(ligand,receptor,sep="-")) %>%
    filter(interaction %in% top_interactions) %>%
    mutate(receiver_cluster = receiver)

  nichenet_plots[[this_rct]] <- make_mushroom_plot_v3(
    base_data,
    rct_label = this_rct,
    first_flag = isFirst,
    show_rankings=TRUE)+
    theme(axis.text.x = element_text(size= 0.1),
          legend.position = "none")

  isFirst <- FALSE

  if(this_rct %in% selected_rcts){
    nichenet_plots_subset[[this_rct]] <- make_mushroom_plot_v3(
      base_data %>% filter(sender %in% selected_scts),
      rct_label = this_rct,
      first_flag = isFirstSubset,
      show_rankings=TRUE)+
      theme(axis.text.x = element_text(size= 0.1),
            legend.position = "none")
    isFirstSubset <- FALSE

    plot_limits_ligand <- list(max = max(base_data$scaled_lfc_ligand),min = min(base_data$scaled_lfc_ligand))
    plot_limits_receptor <- list(max = max(base_data$scaled_lfc_receptor),min = min(base_data$scaled_lfc_receptor))
    plot_limits_decipher <- list(max = 1,min = 0)

    base_data$stroke <- 0.5

    base_data <- base_data %>%
      mutate(stroke_ligand = if_else(scaled_avg_exprs_ligand > 0.01,0.5,NA)) %>%
      mutate(size_ligand = if_else(scaled_avg_exprs_ligand > 0.01,scaled_avg_exprs_ligand,NA)) %>%
      mutate(stroke_receptor = if_else(scaled_avg_exprs_receptor > 0.01,0.5,NA),
             size_receptor = if_else(scaled_avg_exprs_receptor > 0.01,scaled_avg_exprs_receptor,NA),
             size = 1)

    ligand_bubble_plot <- plotBubble(
      df = base_data %>% filter(sender %in% selected_scts),
      x_var = "sender",
      color.var = "scaled_lfc_ligand",
      size.var = "size_ligand",
      stroke.var = "stroke_ligand",
      plot.position = "left",
      col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
      plot.title = "Ligand",
      x_lab= "SCT",
      y_lab = this_rct)

    decipher_bubble_plot <- plotBubble(
      df = base_data %>% filter(sender %in% selected_scts),
      x_var = "sender",
      color.var = "prioritization_score",
      size.var = "size",
      stroke.var = "stroke",
      plot.position = "middle",
      col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
      plot.title = "NicheNet score",
      x_lab= "SCT",
      y_lab = "")

    receptor_bubble_plot <- plotBubble(
      df = base_data %>% filter(receiver %in% selected_rcts),
      x_var = "receiver",
      color.var = "scaled_lfc_receptor",
      size.var = "size_receptor",
      stroke.var = "stroke_receptor",
      plot.position = "middle",
      col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
      plot.title = "Receptor",
      x_lab= "RCT",
      y_lab = "")

    png(paste("nichenet_in_decipher_format_",this_rct,".png",sep=""),width = 24,height = 10,units = "cm",res = 600)
    print(ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot)
    dev.off()
  }
}

png("nichenet_results_full.png",width = 40,height = 60,units= "cm",res= 600)
patchwork::wrap_plots(nichenet_plots,nrow = length(nichenet_plots))
dev.off()

png("nichenet_results_subset.png",width = 10,height = 20,units= "cm",res= 600)
patchwork::wrap_plots(nichenet_plots_subset,nrow = length(nichenet_plots_subset))
dev.off()

#liana
top_interactions_for_liana <- stringr::str_replace(top_interactions,"-","^")
write.csv(top_interactions_for_liana,"top_interactions_for_liana_plot.csv")
write.csv(selected_rcts,"selected_rcts_for_liana_plot.csv")
write.csv(selected_scts,"selected_scts_for_liana_plot.csv")

#liana ----
for(this_rct in selected_rcts){
  base_data <- liana_results %>%
    filter(target == this_rct) %>%
    filter(source %in% selected_scts) %>%
    mutate(interaction = paste(ligand,receptor,sep="-")) %>%
    filter(interaction %in% top_interactions)

  plot_limits_ligand <- list(max =1,min = 0)
  plot_limits_receptor <- list(max = 1,min = 0)
  plot_limits_decipher <- list(max = 5,min = -5)

  base_data$stroke <- 0.5

  base_data <- base_data %>%
    mutate(stroke_ligand = if_else(ligand_expr > 0.01,0.5,NA)) %>%
    mutate(size_ligand = if_else(ligand_expr > 0.01,ligand_expr,NA)) %>%
    mutate(stroke_receptor = if_else(receptor_expr > 0.01,0.5,NA),
           size_receptor = if_else(receptor_expr > 0.01,receptor_expr,NA),
           size = 1)

  ligand_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "source",
    color.var = "ligand_padj",
    size.var = "ligand_expr",
    stroke.var = "stroke_ligand",
    plot.position = "left",
    col.min.val=plot_limits_ligand$min,col.max.val=plot_limits_ligand$max,
    plot.title = "Ligand",
    x_lab= "SCT",
    y_lab = this_rct)

  decipher_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "source",
    color.var = "interaction_stat",
    size.var = "size",
    stroke.var = "stroke",
    plot.position = "middle",
    col.min.val=plot_limits_decipher$min,col.max.val=plot_limits_decipher$max,
    plot.title = "LIANA+ score",
    x_lab= "SCT",
    y_lab = "")

  receptor_bubble_plot <- plotBubble(
    df = base_data,
    x_var = "target",
    color.var = "receptor_padj",
    size.var = "receptor_expr",
    stroke.var = "stroke_receptor",
    plot.position = "middle",
    col.min.val=plot_limits_receptor$min,col.max.val=plot_limits_receptor$max,
    plot.title = "Receptor",
    x_lab= "RCT",
    y_lab = "")

  png(paste("liana_subset_in_decipher_format_",this_rct,".png",sep=""),width = 24,height = 10,units = "cm",res = 600)
  print(ligand_bubble_plot+decipher_bubble_plot+receptor_bubble_plot)
  dev.off()
}


ligands_to_rct
sct_to_ligands

colnames(ligands_to_rct) <- colnames(sct_to_ligands)
ligands_to_rct_test <- ligands_to_rct
sct_to_ligands_test <- sct_to_ligands

ligands_to_rct_test$to <- paste("RCT",ligands_to_rct_test$to,sep="_")
sct_to_ligands_test$from <- paste("SCT",sct_to_ligands_test$from,sep="_")
colnames(ligands_to_rct_test) <- c("Source","Target","Value")
colnames(sct_to_ligands_test) <- c("Source","Target","Value")



links <- bind_rows(ligands_to_rct_test,sct_to_ligands_test)
colnames(links) <- c("Source","Target","Value")

# Create a list of all unique nodes
links <- links %>%
  filter(Value != 0)
nodes <- data.frame(name = unique(c(links$Source,links$Target)))

node_indices <- data.frame(
  index = c(0:(length(nodes$name)-1)),
  name = nodes$name
)


links$Source <- node_indices$index[match(links$Source,node_indices$name)]
links$Target <- node_indices$index[match(links$Target,node_indices$name)]




# Ensure all sources and targets are numeric and there are no NAs
links$Source <- as.numeric(links$Source)
links$Target <- as.numeric(links$Target)
links$Value <- as.numeric(links$Value)

# Remove rows with NAs or zeros in the value
links <- links[!is.na(links$Value) & links$Value != 0, ]

# Check for duplicated nodes
anyDuplicated(nodes$name)
# Plot the Sankey diagram
sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "Source", Target = "Target", Value = "Value", NodeID = "name")



# Assuming 'sankey' is your sankeyNetwork object
saveNetwork(sankey, file = "sankeyPlot.html")




sankey_data <- data.frame(
  SCT = c("SCT1", "SCT1", "SCT2", "SCT2", "SCT3"),
  Ligand = c("Ligand1", "Ligand2", "Ligand2", "Ligand3", "Ligand3"),
  Value_SCT_Ligand = c(10, 20, 5, 15, 25),
  RCT = c("RCT1", "RCT2", "RCT1", "RCT3", "RCT2"),
  Value_Ligand_RCT = c(30, 25, 15, 20, 10)
)

# Create a list of all unique nodes
nodes <- data.frame(name = unique(c(sankey_data$SCT, sankey_data$Ligand, sankey_data$RCT)))

# Create two sets of links for the different flows and bind them into one data frame
links_sct_ligand <- sankey_data %>%
  mutate(Source = match(SCT, nodes$name) - 1,
         Target = match(Ligand, nodes$name) - 1,
         Value = Value_SCT_Ligand)

links_ligand_rct <- sankey_data %>%
  mutate(Source = match(Ligand, nodes$name) - 1,
         Target = match(RCT, nodes$name) - 1,
         Value = Value_Ligand_RCT)

links <- bind_rows(links_sct_ligand, links_ligand_rct) %>%
  select(Source, Target, Value)

# Plot the Sankey diagram
sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "Source", Target = "Target", Value = "Value", NodeID = "name")
sankey



#OLD code (may have some useful stuff) ----

#PANEL B ----
##data wrangling ----
#this used to be called decipher_scores_bound
decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)
decipher_scores_by_cluster_bound_filtered <- decipher_scores_by_cluster_bound %>%
  mutate(decipher_score = sign(decipher_score)*log10(abs(decipher_score)+1)) %>%
  filter(abs(decipher_score) > 0.4)

decipher_scores_matrix_for_heatmap <- reshape2::acast(decipher_scores_by_cluster_bound_filtered,interaction~receiver_cluster,value.var = "decipher_score",fill=0)

circos_data  <- decipher_scores_by_cluster_bound_filtered %>%
  select(interaction,receiver_cluster,decipher_score) %>%
  dplyr::rename(from=interaction, to=receiver_cluster, value=decipher_score)
circos_data$to <- clean_names(circos_data$to)

# Replace 'desired_order' with your actual order of cell types
unique(circos_data$to)
desired_order <- c("B", "CD4_T", "CD8_T","Naive_CD8_T", "NK","pDC","cDC2","CD14+BDCA1+PD-L1+cells","CD14+monocytes","CD16+monocytes","Platelets") # Example order

# Reorder 'receiver_cluster' based on 'desired_order'
circos_data$to <- factor(circos_data$to, levels = desired_order)

## visualization ----
circos.par(start.degree = 90, track.margin = c(0, 0))
# Generate the chord diagram
chordDiagram(circos_data, transparency = 0.5,
             annotationTrack = "grid",
             preAllocateTracks = 1)

# Customize the sector names (labels) orientation and make text smaller
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1],
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6) # Adjust 'cex' for smaller text
}, bg.border = NA)

circos.clear()
dev.off()


##this used to be for a volcano plot ----
random_values <- regulon_deltas_by_cluster$Tumour_Cells %>%
  filter(class == "random") %>%
  pull(deltaPagoda)

real_features <- regulon_deltas_by_cluster$Tumour_Cells %>%
  filter(class == "real")
real_features

real_features$p_value <- calculate_p_value(base_values = random_values, real_features$deltaPagoda)

real_features$log_10 <- -1*log(real_features$p_value,base=10)

ggplot(real_features,aes(x = deltaPagoda,y=log_10))+
  geom_point()

abs_max_tf_delta <- max(abs(real_features$deltaPagoda))

p <- ggplot(real_features, aes(x = deltaPagoda, y = log_10)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey") +
  geom_point(shape = 16, alpha = 0.5) +
  scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
  scale_x_continuous(limits = c(-1*abs_max_tf_delta,abs_max_tf_delta), breaks = seq(-6,6,2)) +
  scale_y_continuous(limits = c(0,15)) +
  xlab("delta TF activity") +
  ylab("-log10(P)") +
  theme_classic(9) +
  theme(legend.position = "none")

p <- p + geom_text(data = subset(real_features, log_10 > -log(0.01,base=10)), aes(label = name),
                   vjust = "inward", hjust = "inward", check_overlap = TRUE,size=2)

# Print the plot
png(file.path(output_data_filepath,"volcano_plot_C8.png"),height = 10,width = 10
    , units = "cm",res=500)
print(p)
dev.off()

list.files(output_data_filepath)

list.files(output_data_filepath)


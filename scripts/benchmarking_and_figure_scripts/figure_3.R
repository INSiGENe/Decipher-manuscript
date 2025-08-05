#Load Manuscript package ----
library(devtools)
load_all()

#Required libraries ----
library(ggplot2)
library(dplyr)
library(gplots)
library(reshape2)
library(igraph)
library(scales)
library(Seurat)
library(stringr)
library(patchwork)


#Parameters ----
##SELECT parameters ----
set.seed(1)
dataset_path <- "results/covid"

pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- file.path("reference_data")
output_data_filepath <- file.path(dataset_path,"data")
figures_folder <- "figures_04_08_2025"

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
slamf7_high <- read.csv("reference_data/SLAMF7_high_signature.csv")


#PANEL A - volcano plot ----
## Data Wrangling ----
regulons_scores_by_clusters_c8 <- regulons_scores_by_clusters$CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells
condition_match <- match(colnames(regulons_scores_by_clusters_c8),names(decipher_seurat_lr$condition))
group_vector <- decipher_seurat_lr$condition[condition_match]

regulon_deltas_c8 <- regulon_deltas_by_cluster$CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells %>%
  filter(class == "real")
regulon_deltas_c8


group_vector[group_vector == "control"] <- 0
group_vector[group_vector == "case"] <- 1
group_factor <- factor(c(group_vector), levels = c(0, 1), labels = c("control", "case"))

diff_regulon_scores_p_values <- do_t_test_by_feature_by_grouping_factor(regulons_scores_by_clusters_c8,group_factor)

regulon_deltas_c8$p_value <- diff_regulon_scores_p_values[regulon_deltas_c8$name]
regulon_deltas_c8$log_10 <-  -1*log(regulon_deltas_c8$p_value,base=10)

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

ggsave(file.path(figures_folder,"figure_3a.png"), plot = p, width = 4, height = 7, dpi = 300)

write.csv(
  regulon_deltas_c8,
  file = file.path(figures_folder, "figure_3a.csv"),
  row.names = TRUE
)

#PANEL B  ----
plotDecipherPrioritizedMap("results/covid",top_n=6,priority_receiver_cells = "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells",dataset_name="figure_3b", width=21,height=11)
#plot is results/covid/figures

#PANEL C ----
##data wrangling ----
decipher_scores_by_cluster_bound <- bind_rows(decipher_scores_by_cluster)
decipher_scores_by_cluster_bound_filtered <- decipher_scores_by_cluster_bound %>%
  mutate(decipher_score = sign(decipher_score)*log10(abs(decipher_score)+1)) %>%
  filter(abs(decipher_score) > 0.4)

decipher_scores_by_cluster_bound_clean <- decipher_scores_by_cluster_bound %>%
  mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
  rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score) %>%
  select(interaction,ligand,receptor,receiver,prioritization_score) %>%
  arrange(prioritization_score)


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

to_plot <- decipher_scores_by_regulon_and_cluster$CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells %>%
  filter(regulon %in% regulon_signature,
         interaction %in% top_interactions)


this_matrix <-reshape2::acast(to_plot,interaction~regulon,value.var = "decipher_score",fill = 0)

##visualization ----
png(file.path(figures_folder,"figure_3c.png"),width = 15, height = 9, units = "cm",res=600)
heatmap.2(
  this_matrix,
  trace="none",
  col = "bluered",
  breaks = 100,
  cexRow = 0.7,
  cexCol=0.5,
  scale = "none",
  key.title = "LR-TF Decipher score",Rowv = FALSE,
  dendrogram = "none",
  margins = c(6,6),
  key = FALSE,
  keysize = 0.3,
  Colv=TRUE)
dev.off()
write.csv(this_matrix,file.path(figures_folder,"figure_3c.csv"),row.names=TRUE)

#PANEL D ----
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

c8_filtered <- de_markers_by_cluster$CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells %>%
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

new_matrix <- new_matrix[,c(7,6,1,2,3,4,5)]

##visualization ----
# Define the color palette from light red to dark red
red_palette <- colorRampPalette(c("lightcoral","red", "darkred"))(256)
divergent_palette <- colorRampPalette(c("lightcoral", "darkred", "purple"))(22)

# 2. compute breaks exactly matching your palette length + 1
brks <- seq(
  min(new_matrix, na.rm = TRUE),
  max(new_matrix, na.rm = TRUE),
  length.out = length(divergent_palette) + 1
)

png(file.path(figures_folder,"figure_3d.png"),height = 12,width=10,units="cm",res=400)
heatmap.2(
  new_matrix,
  trace = "none",
  col = divergent_palette,
  breaks       = brks,  
  symkey       = FALSE,
  density.info = "none",
  scale = "none",
  cexRow = 0.7,
  cexCol = 0.7,
  margin = c(5,5),
  key.title = "log2FC",
  key.xlab     = NA,
  dendrogram = "none",
  Rowv = FALSE,
  Colv= FALSE,
  keysize = 2)
dev.off()
write.csv(new_matrix,file.path(figures_folder,"figure_3d.csv"))

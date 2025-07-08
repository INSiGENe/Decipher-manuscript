
library(devtools)
load_all()

#global settings ----
set.seed(123)

#parameters ----
dataset_name <- "lupus"
dataset_path <- "manuscript_analysis/lupus"
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
meta_path <- "manuscript_analysis/data_for_meta_comparisons"
output_figures_filepath <- file.path(dataset_path,"figures")

reference_filepath <- "reference_data"
nichenet_reference_filepath <- "reference_data/nichenet"

decipher_filepath <- file.path(dataset_path,"data")
nichenet_filepath <- file.path(dataset_path,"nichenet/data")
connectome_filepath <- file.path(dataset_path,"connectome/data")
liana_filepath <- file.path(dataset_path,"liana/data")
natmi_filepath <- file.path(dataset_path,"natmi/data")
cytosig_filepath <- file.path(dataset_path,"cytosig/0_outputs")

dir.create(meta_path,recursive=TRUE)

#FIGURE 2 ----
##load data ----
liana_results <- read.csv(file.path(liana_filepath,"liana_p_interaction_results.csv"),header=TRUE,row.names=1)
nichenet_results <- readRDS(file.path(nichenet_filepath,"nichenet_results.rds"))
nichenet_prior_table_all_clusters <- readRDS(file.path(nichenet_filepath,"prior_table_all_clusters.rds"))
decipher_results <- readRDS(file.path(decipher_filepath,"decipher_scores_by_cluster.rds"))
connectome_results <- readRDS(file.path(connectome_filepath,"connectome_results.rds"))
natmi_results_all <- read.csv(file.path(natmi_filepath,"diff/Delta_edges_lrc2p/All_edges_mean.csv"))

##first data pre-processing ----
natmi_results_pre_processed <- preProcessNATMI(natmi_results_all)
connectome_results_pre_processed <- preProcessConnectome(connectome_results)
liana_pre_processed <- preProcessLIANA(liana_results)

#pre-process Decipher
decipher_bound <- bind_rows(decipher_results)
decipher_pre_processed <- decipher_bound %>%
  #add an interaction column
  mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
  #rename columns of interest
  rename(receiver=receiver_cluster,sender=sender_cluster,prioritization_score=decipher_score)
#scale score for comparison
decipher_pre_processed <- scale_prioritization_score(decipher_pre_processed,"prioritization_score")

#pre-process NicheNet
nichenet_bound <- bind_rows(nichenet_prior_table_all_clusters)
nichenet_pre_processed <- nichenet_bound %>%
  #add an interaction column
  mutate(interaction = paste(ligand, receptor, sep = "-"))
#scale score for comparison
nichenet_pre_processed <- scale_prioritization_score(nichenet_pre_processed,"prioritization_score")

#aggregate into a larger object
results_to_compare_full <- list(
  "NicheNet"=nichenet_pre_processed %>% select(sender,receiver,interaction,prioritization_score,scaled_score),
  "Decipher"=decipher_pre_processed %>% select(sender,receiver,interaction,prioritization_score,scaled_score),
  "LIANA+"=liana_pre_processed %>% select(sender,receiver,interaction,prioritization_score,scaled_score),
  "NATMI"=natmi_results_pre_processed %>% select(sender,receiver,interaction,prioritization_score,scaled_score),
  "Connectome"=connectome_results_pre_processed %>% select(sender,receiver,interaction,prioritization_score,scaled_score)
)
data_frame_names <- names(results_to_compare_full)

# Adding the new column to each data frame in the list
named_data_frames <- Map(add_name_column, results_to_compare_full, data_frame_names)

# Binding all the data frames together
combined_data_frame <- do.call(rbind, named_data_frames)

# Rename the DataFrameName column to method
combined_data_frame <- combined_data_frame %>%
  dplyr::rename(method=DataFrameName)

##data pre-processing for correlation ----
nichenet_for_correlation <- prepareDataForCorrelationAnalysis(nichenet_pre_processed)
liana_for_correlation <- prepareDataForCorrelationAnalysis(liana_pre_processed)
connectome_results_for_correlation <- prepareDataForCorrelationAnalysis(connectome_results_pre_processed)
natmi_results_for_correlation <- prepareDataForCorrelationAnalysis(natmi_results_pre_processed)
#we don't need the same pre-processing for Decipher as our scores are already aggregated
decipher_for_correlation <- decipher_pre_processed %>%
  select(interaction,receiver,prioritization_score)  %>%
  arrange(prioritization_score)

#main analysis -----
results_to_compare <- list(
  "NicheNet"=nichenet_for_correlation,
  "Decipher"=decipher_for_correlation,
  "LIANA+"=liana_for_correlation,
  "NATMI"=natmi_results_for_correlation,
  "Connectome"=connectome_results_for_correlation
)

## Correlation and Search Space ----
interaction_results_correlation_search_space <- getInteractionCorrelationAndSearchSpaceBetweenMethods(results_to_compare)

spearman_matrix <- interaction_results_correlation_search_space$spearman
k_matrix <- interaction_results_correlation_search_space$k_matrix

#save in results
saveRDS(spearman_matrix,file.path(output_figures_filepath,paste(dataset_name,"spearman_matrix.rds")))
saveRDS(k_matrix,file.path(output_figures_filepath,paste(dataset_name,"k_matrix.rds")))

#but also save for further aggregation
saveRDS(spearman_matrix,file.path(meta_path,paste(dataset_name,"spearman_matrix.rds")))
saveRDS(k_matrix,file.path(meta_path,paste(dataset_name,"k_matrix.rds")))

#plot spearman results
plotInteractionCorrelation(spearman_matrix,dataset_name,output_figures_filepath)
plotSearchSpace(k_matrix,dataset_name,output_figures_filepath)

## Overlap ----
### summary table entry ---
overlap_table_row <-getoverlapTable(results_to_compare)
overlap_table_file <- paste(dataset_name,"overlap_table_row.rds")
saveRDS(overlap_table_row,file.path(meta_path,overlap_table_file))

### upset plots -----
##data pre-processing ----
n_top <- 100
liana_set <- getSet(results_to_compare$`LIANA+`,n_top)
nichenet_set <-  getSet(results_to_compare$NicheNet,n_top)
decipher_set <- getSet(results_to_compare$Decipher,n_top)
connectome_set <- getSet(results_to_compare$Connectome,n_top)
natmi_set <- getSet(results_to_compare$NATMI,n_top)

list_input <- list(`LIANA+` = liana_set,
                   NicheNet = nichenet_set,
                   Decipher = decipher_set,
                   Connectome = connectome_set,
                   NATMI = natmi_set)

plotUpsetPlot(list_input,dataset_name,output_figures_filepath)

### summary boxplot entry---
meta_overlap <- getInteractionOverlap(list_input)
overlap_boxplot_file <- paste(dataset_name,"overlap_box_plot.rds")
saveRDS(meta_overlap,file.path(meta_path,overlap_boxplot_file))

##accuracy AUC plot ----
###parameters ----


###load data ----
#LR information
L.set <- getForrestLRDatabase(file.path(reference_filepath,"connectomedb_forrest_lrc2p.csv"))
L.set <- L.set %>%
  mutate(interaction = paste(ligand,receptor,sep="-"),
         lr = interaction) %>%
  unique()

#Cytosig results
z_score_folder <- file.path(cytosig_filepath,"z_score/")
p_value_folder <- file.path(cytosig_filepath,"p_value/")
z_score_files <- list.files(z_score_folder)
p_value_files <- list.files(p_value_folder)

seurat_object_oi <- readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))
mapping_table <- read.csv(file.path(reference_filepath,"cytosig_mapping_table_ligands_genes.csv"),header=TRUE)

###data pre-processing ----
cytosig_significance <- summarizeZScores(z_score_files,z_score_folder,mapping_table)

liana_results_for_comparison <- prepareLianaForCytosigComparison(liana_results)
nichenet_results_for_comparison <- generateComparisonObjectFromNicheNet(nichenet_results)
names(nichenet_results_for_comparison) <- names(nichenet_results)
decipher_results_for_comparison <- lapply(decipher_results, "renameDecipherScore")
connectome_results_for_comparison <- prepareConnectomeForCytosigComparison(connectome_results)
natmi_results_for_comparison <- prepareNatmiForCytosigComparison(natmi_results_all)

###main analysis -----
results_to_compare <- list(
  "NicheNet"=nichenet_results_for_comparison,
  "Decipher"=decipher_results_for_comparison,
  "LIANA+"=liana_results_for_comparison,
  "NATMI"=natmi_results_for_comparison,
  "Connectome"=connectome_results_for_comparison
)

predictions_and_responses <- getPredictionsResponsesForMethods(
  results_to_compare,
  cytosig_significance,
  L.set = L.set,
  seurat_object_oi,
  output_figures_filepath
)

all_predictions_across_methods <- predictions_and_responses$predictions
all_responses_across_methods <- predictions_and_responses$responses

AUC_scores <- plotROCAndExtractAUC(all_predictions_across_methods,all_responses_across_methods,output_figures_filepath)
saveRDS(AUC_scores,file.path(meta_path,paste(dataset_name,"auc_scores.rds")))

###cytosig plots ----
colnames(cytosig_significance) <- convert_text_patterns(colnames(cytosig_significance))
if("CD14+BDCA1+PD-L1+cells" %in% colnames(cytosig_significance)){
  colnames(cytosig_significance)[which(colnames(cytosig_significance) == "CD14+BDCA1+PD-L1+cells" )] <- "C8"
}

#TODO: here if a median z score is greater than absolute (15) then the cell will be greyed out
plotCytosigSignificanceMatrix(cytosig_significance,output_figures_filepath)

# FIGURE 3 ----
#data ----
decipher_scores_by_cluster <- readRDS(file.path(decipher_filepath,"decipher_scores_by_cluster.rds"))
decipher_scores_by_regulon_and_cluster <- readRDS(file.path(decipher_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
regulon_deltas_by_cluster <- readRDS(file.path(decipher_filepath,"regulon_deltas_by_cluster.rds"))
significant_regulons_by_cluster <- readRDS(file.path(decipher_filepath,"significant_regulons_by_cluster.rds"))
seurat_object_oi = readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))
nichenet_results <- readRDS(file.path(nichenet_filepath,"nichenet_results.rds"))
nichenet_prior_table_all_clusters <- readRDS(file.path(nichenet_filepath,"prior_table_all_clusters.rds"))

#reference data----
lr_network = readRDS(file.path(nichenet_reference_filepath,"lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(file.path(nichenet_reference_filepath,"ligand_target_matrix_nsga2r_final.rds"))
weighted_networks = readRDS(file.path(nichenet_reference_filepath,"weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix = readRDS(file.path(nichenet_reference_filepath,"ligand_tf_matrix_nsga2r_final.rds"))


#check for results from older version of Decipher
regulon_deltas_by_cluster <- checkAndAdaptDecipherVersion(regulon_deltas_by_cluster)
#FIGURE 3 ----
## TF activity heatmap ----
top_regulons_and_cts <- getTopRegulonsAndCts(regulon_deltas_by_cluster)
top_regulons_and_cts_delta_matrix <- pull_top_regulons_cts_delta_matrix(regulon_deltas_by_cluster,top_regulons_and_cts)

max_vectors <- apply(top_regulons_and_cts_delta_matrix,MARGIN=1,FUN = "get_abs_max")
#for ERP we had to lower this
top_30_regulons <- max_vectors[order(max_vectors,decreasing=TRUE)][1:30]
#top_30_regulons <- top_30_regulons[-which(is.na(top_30_regulons))]
top_30_regulons_delta_matrix <- top_regulons_and_cts_delta_matrix[names(top_30_regulons),]

# Perform hierarchical clustering
distance_matrix <- dist(t(top_30_regulons_delta_matrix)) # default is euclidean distance
hclust_object_rows <- hclust(distance_matrix) # default is complete linkage
distance_matrix <- dist(top_30_regulons_delta_matrix) # default is euclidean distance
hclust_object_cols <- hclust(distance_matrix) # default is complete linkage

#https://www.tandfonline.com/doi/full/10.1080/21541264.2023.2294623?scroll=top&needAccess=true
tf_colors <- assign_tf_family_colors_to_each_regulon(top_30_regulons_delta_matrix)

# Your heatmap.2 code with modifications
tf_deltas_heatmap <- paste(dataset_name,"delta_heatmap.png")

top_30_regulons_delta_matrix_shortened_labels <- top_30_regulons_delta_matrix
colnames(top_30_regulons_delta_matrix_shortened_labels) <- convert_text_patterns(colnames(top_30_regulons_delta_matrix))

# Generate a sequence of values from low to high to create fixed color breaks
# The number of breaks will determine the resolution of the color scale
color_breaks <- seq(-20, 20, length.out=20)

#mapping for pFizer BioNTech mRNA vaccine
#colnames(top_30_regulons_delta_matrix_shortened_labels)[c(2,3,11)] <- c("CD14+ Mono","CD16+ Mono","C8")
png(file.path(output_figures_filepath,tf_deltas_heatmap),width = 12,height = 8, units = "cm",res=400)
heatmap.2(t(top_30_regulons_delta_matrix_shortened_labels),
          Rowv=as.dendrogram(hclust_object_rows),
          Colv=as.dendrogram(hclust_object_cols),
          trace="none",
          col=colorpanel(19, "orchid4", "white", "chartreuse4"),
          breaks=color_breaks,
          cexRow=0.7,
          cexCol=0.7,
          margins=c(5,5),
          dendrogram="none",
          key=FALSE,
          keysize=0.3,
          ColSideColors=tf_colors,
          xlab = "",
          ylab = "")

dev.off()
saveRDS(top_30_regulons_delta_matrix_shortened_labels,file = file.path(meta_path,paste(dataset_name,"top_30_regulons_delta_matrix_shortened_labels.rds")))

##spearman matrix for all TFS ----
#Figure 3B
spearman_matrix <- GetSpearmanLRTF(decipher_scores_by_regulon_and_cluster)[["spearman_matrix"]]
colnames(spearman_matrix) <- convert_text_patterns(colnames(spearman_matrix))
rownames(spearman_matrix) <- convert_text_patterns(rownames(spearman_matrix))

distance_matrix <- dist(as.matrix(spearman_matrix)) # default is euclidean distance
hclust_object_rows <- hclust(distance_matrix) # default is complete linkage
spearman_matrix <- spearman_matrix[hclust_object_rows$order,hclust_object_rows$order]
plotSpearmanHeatmap(spearman_matrix,labels=NULL,output_filepath = output_figures_filepath,manuscript = FALSE)
saveRDS(spearman_matrix,file = file.path(meta_path,paste(dataset_name,"spearman_lr_tf.rds")))
#for dataset_name used in mansucript, preceed with figure_3_panel_a1_
write.csv(spearman_matrix,file.path(output_figures_filepath,"spearman_regulon_scores.csv"))

##spearman matrix  for overlapping TFs ----
#top 10 shared TFs
top_regulons <- getTopOverlappingRegulons(regulon_deltas_by_cluster,10)
spearman_matrix <- GetSpearmanLRTF(decipher_scores_by_regulon_and_cluster,top_regulons)[["spearman_matrix"]]
colnames(spearman_matrix) <- convert_text_patterns(colnames(spearman_matrix))
rownames(spearman_matrix) <- convert_text_patterns(rownames(spearman_matrix))
label_matrix <- GetSpearmanLRTF(decipher_scores_by_regulon_and_cluster,top_regulons)[["label_matrix"]]
colnames(label_matrix) <- convert_text_patterns(colnames(label_matrix))
rownames(label_matrix) <- convert_text_patterns(rownames(label_matrix))
distance_matrix <- dist(as.matrix(spearman_matrix)) # default is euclidean distance
hclust_object_rows <- hclust(distance_matrix) # default is complete linkage
spearman_matrix <- spearman_matrix[hclust_object_rows$order,hclust_object_rows$order]
plotSpearmanHeatmap(spearman_matrix,labels=label_matrix,output_filepath = output_figures_filepath,manuscript = FALSE)
saveRDS(spearman_matrix,file = file.path(meta_path,paste(dataset_name,"spearman_lr_top_overlapping_tf.rds")))
#for dataset_name used in mansucript, preceed with figure_3_panel_a1_
write.csv(spearman_matrix,file.path(output_figures_filepath,"spearman_regulon_scores_top_overlapping_tf.csv"))
print("finished generating images")

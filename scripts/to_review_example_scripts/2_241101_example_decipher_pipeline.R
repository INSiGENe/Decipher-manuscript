#main Decipher analysis pipeline
#V1
#load decipher package -----
library(devtools)
load_all()

#global options ----
selected_random_seed = 123
set.seed(selected_random_seed)

#global directories
analysis_path <- "older_stim_vs_media"
pre_processing_path <- file.path(analysis_path,"pre_processing")
seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))
this_comparison_path <- "older_stim_vs_media"
#Parameters: dataset ----
#targeted run
sender_cts <- c("ISG_EM_CD4","Th2","EM_Sen_CD4","ISG_Naive_CD4","HSP_CD4","EM_CD4","ISG_CM_CD4","CM_CD4","Naive_CD4","EM_CCR4_CD4")
receiver_cts <- c("Age_B","Mem_B","Naive_B")
#user parameters
case_condition <- "older_stim"
control_condition <- "older_media"
min_cells_per_cluster_condition <- 100
species <-  "human"
#select which comparison to run
condition_name <- "condition"
k_parameter = 1
flag.normalize.non.log <- FALSE
#local directories ----
pre_processing_path <- file.path(this_comparison_path,"pre_processing")
reference_filepath <- "reference_data"
output_filepath <- this_comparison_path
output_data_filepath <- file.path(output_filepath,"data")
output_figures_filepath <- file.path(output_filepath,"figures")
output_importances_filepath <- file.path(output_filepath,"importances")
#directory set up----
dir.create(pre_processing_path)
dir.create(file.path(pre_processing_path,"h5ad_by_cluster"))
dir.create(output_data_filepath,recursive=TRUE)
dir.create(output_figures_filepath,recursive=TRUE)
dir.create(output_importances_filepath,recursive=TRUE)

#load dataset ----
seurat_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

#data pre-processing ----
#WARNING: if not using the sample dataset, please run the following two functions
seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi,condition_name,case_condition,control_condition)


#load reference data ----
if(species == "human"){
  L.set <- loadLSet(reference_filepath,species)

} else {
  L.set.mouse <- readRDS("reference_data/mouse_lr_pair.rds")

  L.set <- L.set.mouse %>%
    # Rename the columns
    rename(
      ligand = ligand_gene_symbol,
      receptor = receptor_gene_symbol
    ) %>%
    # Create the 'lr' and 'interaction' columns by combining ligand and receptor
    mutate(
      lr = paste(ligand, receptor, sep = "-"),
      interaction = lr  # 'interaction' is the same as 'lr'
    ) %>%
    # Select the desired columns in the correct order
    select(ligand, receptor, lr, interaction)
}


cytosig_ligands <- loadCytosigLigands(reference_filepath,species)

##############
##QC ----
##############
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
plotQC_CpC(CpC_data,outputPath=output_figures_filepath)

#PARAM: select the minimum number of cells per cluster + condition
#clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,minCpc = 100)
#seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed=NULL]
plotQC_UpC(seuratObject = seurat_oi,outputPath = output_figures_filepath,id = "_sc")

min_meta_cells_parameter = 100
##############
##Meta cells ----
##############
decipher_seurat <- metaCellModule(
  seurat_object = seurat_oi[, which(seurat_oi$cluster %in% receiver_cts), seed=NULL],
  min_meta_cells = min_meta_cells_parameter,
  k = k_parameter
)

plotQC_UpC(seuratObject = decipher_seurat,outputPath = output_figures_filepath,id = "_meta")

CpC_data_meta <- generateQCDataByClusterAndCondition(decipher_seurat,max(stringr::str_length(unique(decipher_seurat$cluster))))
plotQC_CpC(CpC_data_meta,outputPath=output_figures_filepath,id = "_meta")


parameter_record <- data.frame(
  "k" = k_parameter,
  "min_meta_cells" = min_meta_cells_parameter
)
write.csv(parameter_record,file.path(output_data_filepath,"parameter_record.csv"))

saveRDS(decipher_seurat,file.path(output_data_filepath,"pseudobulk_seurat.rds"))

##############
#data pre-processing: main analysis ----
##############

seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% sender_cts), seed=NULL]

decipher_seurat_lr <- decipher_seurat[unique(c(L.set$ligand,L.set$receptor)),, seed=NULL]

feature_statistics_receptor <- getFeatureStatistics(
  features=unique(c(L.set$receptor)),
  seuratObj=decipher_seurat)


feature_statistics_ligand <- getFeatureStatistics(
  features=unique(c(L.set$ligand)),
  seuratObj=seurat_oi)

feature_statistics <- bind_rows(feature_statistics_ligand,feature_statistics_receptor)

expressed_ligands <- getFilteredLigands(
  seurat_oi,
  L.set,
  param_min_ligand_expr_in_cluster = 0.1)

expressed_receptors_all_clusters <- getExpressedReceptorsForEachCluster(
  decipher_seurat,
  L.set)

L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(
  L.set,
  expressed_ligands,
  expressed_receptors_all_clusters)

regulon_grns_by_cluster <- getRegulonsAllClusters(
  output_filepath,
  decipher_seurat)

capped_regulons_all_clusters <- capRegulonsAllClusters(
  regulon_grns_by_cluster,
  decipher_seurat,
  flag.normalize.non.log)

#used to be called regulon_scores_this_cluster within the loop
#does not pass identical test but I did a spot check and it looked identical, likely number formatting (despite both being doubles)
regulon_scores_by_cluster <- getRegulonScoresAllClusters(
  capped_regulons_all_clusters,
  decipher_seurat)

regulon_deltas_by_cluster <- getRegulonDeltasAllClusters(
  regulon_scores_by_cluster,
  decipher_seurat)

significant_regulons_by_cluster <- getSignificantRegulonsAllClusters(
  regulon_deltas_by_cluster)

significant_regulon_markers_by_cluster <- getDifferentiallyExpressedTargetsForRegulonsAllClusters(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_grns_by_cluster,
  flag.normalize.non.log,
  random.seed=selected_random_seed)

#used to be called interaction_potentials_matrix_this_cluster
#careful with this one, though it looks fine, just double check a few times
interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClusters(
  decipher_seurat,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log,
  sc_seurat = seurat_oi)

interaction_deltas_by_cluster <- calculateInteractionDeltasAllClusters(
  interaction_potential_by_clusters,
  decipher_seurat_lr)

filtered_interaction_potentials_matrix_all_clusters <- filterIntPotByDeltas(
  interaction_potential_by_clusters,
  interaction_deltas_by_cluster)


#careful with this one, not identical but I believe this is due to number encoding
#used to be called interaction_potentials_matrix_clusters
interaction_potentials_matrix_clusters_all_clusters <-
  getInteractionPotentialMatrixForRepresentativeInteractionsAllClusters(
    decipher_seurat,
    L_set_relevant_features_all_clusters,
    filtered_interaction_potentials_matrix_all_clusters,
    cytosig_ligands = cytosig_ligands,
    flag.normalize.non.log)

#so the seed has to be set before randomForest for reproducibility, is this ok SEB?
decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClusters(
  decipher_seurat,
  significant_regulon_deltas_all_clusters = significant_regulons_by_cluster,
  regulon_scores_by_cluster,
  interaction_potentials_matrix_clusters_all_clusters,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log)

lr_markers_by_cluster <- FindLRMarkersAllClusters(
  decipher_seurat,
  decipher_scores_by_regulon_and_cluster,
  flag.normalize.non.log,
  random.seed = selected_random_seed
)

l_markers_by_cluster_sc <- FindLMarkersAllClustersSC(
  seurat_oi,
  expressed_ligands,
  flag.normalize.non.log,
  random.seed = selected_random_seed
)

de_markers_by_cluster <- FindMarkersAllClusters(
  decipher_seurat,
  flag.normalize.non.log,
  random.seed= selected_random_seed
)

#DECIPHER analysis-----
decipher_scores_by_regulon_and_cluster <- lapply(
  decipher_scores_by_regulon_and_cluster,
  FUN = "listOfDFsRenameColumn",
  original_name = "weighted.spearman.cont",
  new_name = "decipher_score")

decipher_scores_by_cluster <- lapply(
  decipher_scores_by_regulon_and_cluster,
  FUN = "calculateScoresByCluster")

decipher_scores_by_cluster <- addListNameToDFElements(
  decipher_scores_by_cluster,
  "receiver_cluster")

#save DECIPHER ----
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(regulon_scores_by_cluster,file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
saveRDS(interaction_potential_by_clusters,file.path(output_data_filepath,"interaction_potential_by_clusters.rds"))
saveRDS(regulon_deltas_by_cluster,file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
saveRDS(significant_regulons_by_cluster,file.path(output_data_filepath,"significant_regulons_by_cluster.rds"))
saveRDS(significant_regulon_markers_by_cluster,file.path(output_data_filepath,"significant_regulon_markers_by_cluster.rds"))
saveRDS(interaction_deltas_by_cluster,file.path(output_data_filepath,"interaction_deltas_by_cluster.rds"))
saveRDS(regulon_grns_by_cluster,file.path(output_data_filepath,"regulon_grns_by_cluster.rds"))
saveRDS(lr_markers_by_cluster,file.path(output_data_filepath,"lr_markers_by_cluster.rds"))
saveRDS(l_markers_by_cluster_sc,file.path(output_data_filepath,"l_markers_by_cluster_sc.rds"))
saveRDS(de_markers_by_cluster,file.path(output_data_filepath,"de_markers_by_cluster.rds"))
saveRDS(feature_statistics,file.path(output_data_filepath,"feature_statistics.rds"))
saveRDS(feature_statistics_ligand,file.path(output_data_filepath,"feature_statistics_ligand.rds"))
saveRDS(feature_statistics_receptor,file.path(output_data_filepath,"feature_statistics_receptor.rds"))
saveRDS(decipher_seurat_lr,file.path(output_data_filepath,"decipher_seurat_lr.rds"))
saveRDS(L.set,file.path(output_data_filepath,"L_set.rds"))
saveRDS(decipher_scores_by_cluster,file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
saveRDS(interaction_potentials_matrix_clusters_all_clusters,file.path(output_data_filepath,"interaction_potentials_matrix_clusters_all_clusters.rds"))
saveRDS(expressed_receptors_all_clusters, file.path(output_data_filepath, "expressed_receptors_all_clusters.rds"))
saveRDS(capped_regulons_all_clusters, file.path(output_data_filepath, "capped_regulons_all_clusters.rds"))
saveRDS(L_set_relevant_features_all_clusters, file.path(output_data_filepath, "L_set_relevant_features_all_clusters.rds"))

#plot results ----

plotDecipherPrioritizedMap(
  this_comparison_path,
  top_n = 6,
  sc_feature_statistics = TRUE,
  selected_receiver_cells = NULL,
  primary_ct = NULL,
  dataset_name = "",
  split_by_direction = FALSE,
  direction = "pos")

plotDecipherPrioritizedMap(
  this_comparison_path,
  top_n = 12,
  sc_feature_statistics = TRUE,
  selected_receiver_cells = NULL,
  primary_ct = NULL,
  dataset_name = "",
  split_by_direction = TRUE,
  direction = "neg")

plotDecipherPrioritizedMap(
  this_comparison_path,
  top_n = 12,
  sc_feature_statistics = TRUE,
  selected_receiver_cells = NULL,
  primary_ct = NULL,
  dataset_name = "",
  split_by_direction = TRUE,
  direction = "pos")

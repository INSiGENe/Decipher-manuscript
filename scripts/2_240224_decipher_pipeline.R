#main Decipher analysis pipeline
#load decipher package -----
library(devtools)
load_all()

#global options ----
set.seed(123)

#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
species <-  "human"
case_condition = "stim"
control_condition = "ctrl"

#Parameters: directories ----
dataset_path <- "sample_analysis"
dir.create(dataset_path)
pre_processing_path <- file.path(dataset_path,"pre_processing")
reference_filepath <- "reference_data"
output_filepath <- dataset_path
output_data_filepath <- file.path(output_filepath,"data")
output_figures_filepath <- file.path(output_filepath,"figures")
output_importances_filepath <- file.path(output_filepath,"importances")
#directory set up----
dir.create(pre_processing_path)
dir.create(file.path(pre_processing_path,"h5ad_by_cluster"))
dir.create(output_data_filepath,recursive=TRUE)
dir.create(output_figures_filepath,recursive=TRUE)
dir.create(output_importances_filepath,recursive=TRUE)

#Parameters: analysis ----
flag.normalize.non.log <- FALSE

#create sample dataset ----
#including seurat object and h5ad objects
seurat_oi <- generateSampleSeuratFromExperimentHub(min_cells_per_cluster_condition,case_condition,control_condition)
##save outputs for Decipher analysis
saveRDS(seurat_oi,file.path("sample_analysis/pre_processing","seurat_object_oi.rds"))
#in addition, we need to create python-compatible h5ad objects for the CO pipeline, here, I've opted against it
#as they are not necessary for this script
# for(this_cluster in unique(kang.seurat$cluster)){
#   seurat_object_oi_this_cluster <- subset(kang.seurat,subset = cluster == this_cluster)
#   sce.object = as.SingleCellExperiment(seurat_object_oi_this_cluster)
#   sce.object@assays@data[["logcounts"]] <- NULL
#   writeH5AD(sce.object, file.path("pre_processing/h5ad_by_cluster",paste(this_cluster,".h5ad",sep="")),X_name = "counts")
# }

#load reference data ----
L.set <- loadLSet(reference_filepath,species)
enrichr_database <- loadEnrichrDatabase(reference_filepath,species)
cytosig_ligands <- loadCytosigLigands(reference_filepath,species)

#data pre-processing ----
#moved this functions to generateSampleSeuratFromExperimentHub() but need alternative when user actually starts with seurat object
#seurat_oi$orig.condition <- seurat_oi$condition
#map conditions to case and control because the code internally has 'case' and 'control'references
#seurat_oi <- mapConditionsInSeurat(seurat_oi,"condition",case_condition,control_condition)

##############
##QC ----
##############
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
#plotQC_CpC(CpC_data,outputPath=output_figures_filepath)

#PARAM: select the minimum number of cells per cluster + condition
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,minCpc = 100)
seurat_oi <- subset(seurat_oi,subset = cluster %in% clusters_passing_CpC_filter)

##############
##Meta cells ----
##############
decipher_seurat <- metaCellModule(
  seurat_object = seurat_oi,
  min_meta_cells = 100
)

saveRDS(decipher_seurat,file.path(output_data_filepath,"pseudobulk_seurat.rds"))

##############
#data pre-processing: main analysis ----
##############
decipher_seurat_lr <- subset(decipher_seurat,features = unique(c(L.set$ligand,L.set$receptor)))

feature_statistics <- getFeatureStatistics(
  features=unique(c(L.set$ligand,L.set$receptor)),
  seuratObj=decipher_seurat)

expressed_ligands <- getFilteredLigands(
  decipher_seurat,
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
  flag.normalize.non.log)

#used to be called interaction_potentials_matrix_this_cluster
#careful with this one, though it looks fine, just double check a few times
interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClusters(
  decipher_seurat,decipher_seurat_this_cluster,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log)

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
    cytosig_ligands,
    flag.normalize.non.log)

#so the seed has to be set before randomForest for reproducibility, is this ok SEB?
decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClusters(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_scores_by_cluster,
  interaction_potentials_matrix_clusters_all_clusters,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log)

lr_markers_by_cluster <- FindLRMarkersAllClusters(
  decipher_seurat,
  decipher_scores_by_regulon_and_cluster,
  flag.normalize.non.log
)

de_markers_by_cluster <- FindMarkersAllClusters(
  decipher_seurat,
  flag.normalize.non.log
)

#this function takes a while so would be best to add a progress bar for the user
# enrichr_results_by_cluster <- enrichResultsAllClusters(
#   de_markers_by_cluster,
#   significant_regulons_by_cluster,
#   regulon_grns_by_cluster,
#   enrichr_database)

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
saveRDS(de_markers_by_cluster,file.path(output_data_filepath,"de_markers_by_cluster.rds"))
saveRDS(feature_statistics,file.path(output_data_filepath,"feature_statistics.rds"))
saveRDS(decipher_seurat_lr,file.path(output_data_filepath,"decipher_seurat_lr.rds"))
saveRDS(L.set,file.path(output_data_filepath,"L_set.rds"))
saveRDS(decipher_scores_by_cluster,file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
saveRDS(interaction_potentials_matrix_clusters_all_clusters,file.path(output_data_filepath,"interaction_potentials_matrix_clusters_all_clusters.rds"))
saveRDS(expressed_receptors_all_clusters, file.path(output_data_filepath, "expressed_receptors_all_clusters.rds"))
saveRDS(capped_regulons_all_clusters, file.path(output_data_filepath, "capped_regulons_all_clusters.rds"))
saveRDS(L_set_relevant_features_all_clusters, file.path(output_data_filepath, "L_set_relevant_features_all_clusters.rds"))

#saveRDS(enrichr_results_by_cluster,file.path(output_data_filepath,"enrichr_results_by_cluster.rds"))


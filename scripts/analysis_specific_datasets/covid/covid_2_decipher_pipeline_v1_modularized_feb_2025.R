
library(devtools)
load_all()

#global options ----
#Set this seed to NULL if you don't need reproducible results
selected_random_seed = 123
set.seed(selected_random_seed)

#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
species <-  "human"
#for sample dataset, dataset_path is "sample_analysis", condition_name is "condition", case_condition is "stim" and control_condition is "ctrl"
#if using another dataset, ensure the following four variables are updated and data pre-processing lines below (marked WARNING) are uncommented
dataset_path <- "results/covid"
condition_name <- "condition"
case_condition <- 22
control_condition <- 0
k_parameter = 2

#Parameters: directories ----
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

#load dataset ----
#replace dataset_path above if using another dataset
seurat_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

#data pre-processing ----
#WARNING: if not using the sample dataset, please run the following two functions
seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi,condition_name,case_condition,control_condition)


#load reference data ----
L.set <- loadLSet(reference_filepath,species)
enrichr_database <- loadEnrichrDatabase(reference_filepath,species)
cytosig_ligands <- loadCytosigLigands(reference_filepath,species)

##############
##QC ----
##############
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
plotQC_CpC(CpC_data,outputPath=output_figures_filepath)

#PARAM: select the minimum number of cells per cluster + condition
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,minCpc = 100)
#here he retain the C8 cluster as it is of interest to us
clusters_passing_CpC_filter<- c(clusters_passing_CpC_filter,"CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells")

seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed=NULL]
plotQC_UpC(seuratObject = seurat_oi,outputPath = output_figures_filepath,id = "_sc")

min_meta_cells_parameter = 100


##############
##Pairings ----
##############

paramPairings <- data.frame(
  case = clusters_passing_CpC_filter,
  control = clusters_passing_CpC_filter
)
paramPairings <- paramPairings %>%
  mutate(control = if_else(case == "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells","CD14_plus_monocytes",control))

min_counts <- seurat_oi@meta.data %>%
  group_by(cluster) %>%
  count(condition) %>%
  ungroup()


groups <- createGroupsFromPairings(paramPairings)
paramPairings_min_n <- calculateMinimumN(groups,min_counts,paramPairings)

##############
##Meta cells ----
##############

seurat_oi@meta.data$barcode <- rownames(seurat_oi@meta.data)
MetaCellMatrices <- generateMetaCellMatricesWPairings(
  seuratObj = seurat_oi,
  paramMaxScCells = 1200*(k_parameter+1),
  paramK = k_parameter,
  paramPairings = paramPairings_min_n)

seurat_pseudo_bulk <- generatePseudoBulkSeurat(
  pseudobulkList = MetaCellMatrices)

decipher_seurat <- Seurat::NormalizeData(seurat_pseudo_bulk,normalization.method="RC",scale.factor = 100000)
rm(MetaCellMatrices)
rm(seurat_pseudo_bulk)
DefaultAssay(decipher_seurat) <- "RNA"

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
decipher_seurat_lr <- decipher_seurat[unique(c(L.set$ligand,L.set$receptor)),, seed=NULL]

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
regulon_scores_by_cluster <- getRegulonScoresAllClustersWParamPairings(
  capped_regulons_all_clusters,
  decipher_seurat,
  paramPairings = paramPairings)

regulon_deltas_by_cluster <- getRegulonDeltasAllClustersWParamPairings(
  regulon_scores_by_cluster,
  decipher_seurat,
  paramPairings)

significant_regulons_by_cluster <- getSignificantRegulonsAllClusters(
  regulon_deltas_by_cluster)

significant_regulon_markers_by_cluster <- getDifferentiallyExpressedTargetsForRegulonsAllClustersWParamPairings(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_grns_by_cluster,
  flag.normalize.non.log,
  paramPairings,
  random.seed = selected_random_seed)

#used to be called interaction_potentials_matrix_this_cluster
#careful with this one, though it looks fine, just double check a few times
interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClustersWParamPairings(
  decipher_seurat,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log,
  paramPairings)

interaction_deltas_by_cluster <- calculateInteractionDeltasAllClusters(
  interaction_potential_by_clusters,
  decipher_seurat_lr)

filtered_interaction_potentials_matrix_all_clusters <- filterIntPotByDeltas(
  interaction_potential_by_clusters,
  interaction_deltas_by_cluster)

#careful with this one, not identical but I believe this is due to number encoding
#used to be called interaction_potentials_matrix_clusters
interaction_potentials_matrix_clusters_all_clusters <-
  getInteractionPotentialMatrixForRepresentativeInteractionsAllClustersWParamPairings(
    decipher_seurat,
    L_set_relevant_features_all_clusters,
    filtered_interaction_potentials_matrix_all_clusters,
    cytosig_ligands,
    flag.normalize.non.log,
    paramPairings)

#so the seed has to be set before randomForest for reproducibility, is this ok SEB?
decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClustersWParamPairings(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_scores_by_cluster,
  interaction_potentials_matrix_clusters_all_clusters,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log,
  paramPairings)

#missing these with paramPairings
lr_markers_by_cluster <- FindLRMarkersAllClustersWParamPairings(
  decipher_seurat,
  decipher_scores_by_regulon_and_cluster,
  flag.normalize.non.log,
  paramPairings,
  random.seed = selected_random_seed
)

de_markers_by_cluster <- FindMarkersAllClustersWParamPairings(
  decipher_seurat,
  flag.normalize.non.log,
  paramPairings,
  random.seed = selected_random_seed
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

#plot results ----
plotDecipherPrioritizedMap(dataset_path,top_n=6,dataset_name="covid")
plotDecipherPrioritizedMap_v2(
  dataset_path = dataset_path,
  receiver_cell_type = "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells",
  output_filename = "decipher_plot_prioritized_c8_log",
  log_transform = TRUE,
  slice_n = 11)
#saveRDS(enrichr_results_by_cluster,file.path(output_data_filepath,"enrichr_results_by_cluster.rds"))


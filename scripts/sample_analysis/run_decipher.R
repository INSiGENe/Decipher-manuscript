renv::restore()
library(devtools)
load_all()
devtools::install_github("immunogenomics/presto")
library(presto)
#global options ----
i = 123
set.seed(i)
selected_random_seed <- i

#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
species <-  "human"
#for sample dataset condition_name is "condition", case_condition is "stim" and control_condition is "ctrl"
condition_name <- "condition"
case_condition = "stim"
control_condition = "ctrl"
k_parameter = 2
min_meta_cells_parameter = 100



create_project_dirs <- function(dataset_path) {
  dirs <- c(
    dataset = dataset_path,
    pre_processing = file.path(dataset_path, "pre_processing"),
    co_input = file.path(dataset_path,'pre_processing',"h5ad_by_cluster"),
    data = file.path(dataset_path, "data"),
    figures = file.path(dataset_path, "validity", "figures"),
    importances = file.path(dataset_path, "validity", "importances")
  )

  dir.create(dataset_path, recursive = TRUE, showWarnings = FALSE)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  invisible(dirs)   # returns the paths (named) invisibly
}

#Parameters: directories ----
reference_filepath <- "reference_data"
dataset_path <- "results/sample_analysis"
paths <- create_project_dirs("results/sample_analysis")
#not sure what this is about: dir.create(file.path(paths['pre_processing'],"validity/h5ad_by_cluster"))

#Parameters: analysis ----
flag.normalize.non.log <- FALSE

# #create sample dataset ----
# #including seurat object and h5ad objects
seurat_oi <- readRDS(file.path(paths["data"],"seurat_object_oi.rds"))

#load reference data ----
L.set <- loadLSet(reference_filepath,species)
cytosig_ligands <- loadCytosigLigands(reference_filepath,species)

#data pre-processing ----
#moved this functions to generateSampleSeuratFromExperimentHub() but need alternative when user actually starts with seurat object
seurat_oi$orig.condition <- seurat_oi[["condition"]]
#map conditions to case and control because the code internally has 'case' and 'control'references

##############
##QC ----
##############
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
plotQC_CpC(CpC_data,outputPath=paths['figures'])

#PARAM: select the minimum number of cells per cluster + condition
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,minCpc = 100)
seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed=NULL]
plotQC_UpC(seuratObject = seurat_oi, outputPath = paths['figures'],id = "_sc")
##############
##Meta cells ----
##############
decipher_seurat <- metaCellModule(
  seurat_object = seurat_oi,
  min_meta_cells = min_meta_cells_parameter,
  k = k_parameter
)

plotQC_UpC(seuratObject = decipher_seurat,outputPath = paths['figures'],id = "_meta")

CpC_data_meta <- generateQCDataByClusterAndCondition(decipher_seurat,max(stringr::str_length(unique(decipher_seurat$cluster))))
plotQC_CpC(CpC_data_meta,outputPath=paths['figures'],id = "_meta")


parameter_record <- data.frame(
  "k" = k_parameter,
  "min_meta_cells" = min_meta_cells_parameter
)
write.csv(parameter_record,file.path(paste0(paths['data'],"/",i,"_parameter_record.csv")))

saveRDS(decipher_seurat,file.path(paste0(paths['data'],"/",i,"_pseudobulk_seurat.rds")))

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
  dataset_path,
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

saveRDS(decipher_scores_by_cluster,file.path(paths['data'],paste0("Run_",i , "_decipher_scores_by_cluster.rds")))
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(paths['data'],paste0("Run_",i , "_decipher_scores_by_regulon_and_cluster.rds")))



dataset_path <- file.path(paths['data'],"for_plotting")
paths['figures'] <- file.path(dataset_path,"figures")
paths['data'] <- file.path(dataset_path,"data")
dir.create(paths['data'],recursive = TRUE)
dir.create(paths['figures'],recursive=TRUE)
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(paths['data'],"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(regulon_scores_by_cluster,file.path(paths['data'],"regulon_scores_by_cluster.rds"))
saveRDS(interaction_potential_by_clusters,file.path(paths['data'],"interaction_potential_by_clusters.rds"))
saveRDS(regulon_deltas_by_cluster,file.path(paths['data'],"regulon_deltas_by_cluster.rds"))
saveRDS(significant_regulons_by_cluster,file.path(paths['data'],"significant_regulons_by_cluster.rds"))
saveRDS(significant_regulon_markers_by_cluster,file.path(paths['data'],"significant_regulon_markers_by_cluster.rds"))
saveRDS(interaction_deltas_by_cluster,file.path(paths['data'],"interaction_deltas_by_cluster.rds"))
saveRDS(regulon_grns_by_cluster,file.path(paths['data'],"regulon_grns_by_cluster.rds"))
saveRDS(lr_markers_by_cluster,file.path(paths['data'],"lr_markers_by_cluster.rds"))
saveRDS(de_markers_by_cluster,file.path(paths['data'],"de_markers_by_cluster.rds"))
saveRDS(feature_statistics,file.path(paths['data'],"feature_statistics.rds"))
saveRDS(decipher_seurat_lr,file.path(paths['data'],"decipher_seurat_lr.rds"))
saveRDS(L.set,file.path(paths['data'],"L_set.rds"))
saveRDS(decipher_scores_by_cluster,file.path(paths['data'],"decipher_scores_by_cluster.rds"))
saveRDS(interaction_potentials_matrix_clusters_all_clusters,file.path(paths['data'],"interaction_potentials_matrix_clusters_all_clusters.rds"))
saveRDS(expressed_receptors_all_clusters, file.path(paths['data'], "expressed_receptors_all_clusters.rds"))
saveRDS(capped_regulons_all_clusters, file.path(paths['data'], "capped_regulons_all_clusters.rds"))
saveRDS(L_set_relevant_features_all_clusters, file.path(paths['data'], "L_set_relevant_features_all_clusters.rds"))


#plot for a particular seeda
plotDecipherPrioritizedMap(dataset_path,top_n=4,dataset_name="sample_1", abs_decipher_plot_limit = 20,width=21,height=9)
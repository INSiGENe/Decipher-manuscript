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

##output objects initialize ----
ligand_scores_result <- list()
decipher_scores_by_regulon_and_cluster <- list()
regulon_scores_by_cluster <- list()
interaction_potential_by_clusters <- list()
regulon_deltas_by_cluster <- list()
significant_regulons_by_cluster <- list()
significant_regulon_markers_by_cluster <- list()
regulon_grns_by_cluster <- list()
lr_markers_by_cluster <- list()
de_markers_by_cluster <- list()
enrichr_results_by_cluster <- list()
interaction_deltas_by_cluster <- list()

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

#used to be called expressed_receptors within the loop
expressed_receptors_all_clusters <- getExpressedReceptorsForEachCluster(decipher_seurat,L.set)

#used to be called L_set_relevant_features within the loop
L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(L.set,expressed_ligands,expressed_receptors_all_clusters)

#used to be called regulon_this_cluster within the loop
regulons_all_clusters <- getRegulonsAllClusters(output_filepath,decipher_seurat)

#used to be called regulon_this_cluster_capped within the loop
capped_regulons_all_clusters <- capRegulonsAllClusters(regulons_all_clusters,decipher_seurat,flag.normalize.non.log)

#used to be called regulon_scores_this_cluster within the loop
#does not pass identical test but I did a spot check and it looked identical, likely number formatting (despite both being doubles)
regulon_scores_all_clusters <- getRegulonScoresAllClusters(capped_regulons_all_clusters,decipher_seurat)

#DECIPHER analysis-----
start_time <- Sys.time()
for(this_cluster in unique(decipher_seurat$cluster)[1]){

  #main object
  decipher_seurat_this_cluster <- subset(decipher_seurat,subset = cluster == this_cluster)
  #see how I can simplify this, I don't think I should use downsampled...
  SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition
  data_this_cluster <- decipher_seurat_this_cluster@assays$RNA@data

  if(flag.normalize.non.log){
    decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster,normalization.method = "RC",scale.factor=100000)
  }

  data_this_cluster_receptors <- data_this_cluster[which(rownames(data_this_cluster) %in% unique(L_set_relevant_features$receptor)),]

  ##PAGODA -----
  #TODO: silence this function
  regulon_scores_this_cluster <- getRegulonScores(
    seuratObject = decipher_seurat_this_cluster,
    grn_df = regulon_this_cluster_capped)

  ##PAGODA DELTA ----
  regulon_deltas_this_cluster <- getRegulonDeltas(
    regulon_scores_this_cluster,
    decipher_seurat_this_cluster$condition)

  significant_regulon_deltas_this_cluster <- getSignificantRegulons(regulon_deltas_this_cluster)

  #### find target genes for each top differentially expressed regulons and calculate diff expr. ----
  SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster$condition
  #wait but this needs to align to my GRN right?
  significant_regulon_markers_by_cluster[[this_cluster]] <- getDifferentiallyExpressedTargetsForRegulons(
    seuratObj = decipher_seurat_this_cluster,
    regulonNames = significant_regulon_deltas_this_cluster$name,
    logFcThreshold = 0.58,
    grnDf = regulon_this_cluster,
    targetCt = this_cluster
  )

  ## calculate Interaction Potential Matrix ----
  interaction_potentials_matrix_this_cluster <- getInteractionPotentialsMatrixThisCluster(
    seurat_obj = decipher_seurat,
    seurat_obj_this_cluster_ds = decipher_seurat_this_cluster,
    selected_lr_pairs = L_set_relevant_features
  )

  #this cannot be moved from this location
  interaction_deltas <- calculateInteractionDeltas(interaction_potentials_matrix_this_cluster,decipher_seurat_lr)
  interaction_deltas_by_cluster[[this_cluster]] <- interaction_deltas

  #subset interaction potential matrix to those interactions that have changed between conditions
  interaction_potentials_matrix_this_cluster <- interaction_potentials_matrix_this_cluster[rownames(interaction_deltas),]

  ## subset interaction_potential matrix by correlation clusters for cluster-based RF ----
  interaction_potentials_matrix_clusters <- getInteractionPotentialMatrixForRepresentativeInteractions(
    data_this_cluster_receptors,
    selected_lr_pairs = L_set_relevant_features,
    interaction_potentials_matrix_this_cluster,
    cytosig_ligands
  )

  ## run random forest on each regulon -----
  all_rf_results <- list()
  for(this.tf in significant_regulon_deltas_this_cluster$name){

    ind.this.tf <- which(significant_regulon_deltas_this_cluster$name == this.tf)

    val.this.tf <- significant_regulon_deltas_this_cluster$deltaPagoda[ind.this.tf]

    print(paste("calculating forest for",this.tf))

    tf.merged <- regulon_scores_this_cluster[this.tf,colnames(interaction_potentials_matrix_clusters)]

    rf <- randomForest::randomForest(
      x = t(interaction_potentials_matrix_clusters),
      y=tf.merged,
      ntree = 100,
      importance=T)

    imp.df <- extractDecipherResults(
      random_forest_results = rf,
      interaction_potentials_matrix_clusters,
      data_this_cluster_receptors,
      selected_lr_pairs = L_set_relevant_features,
      this.tf,
      val.this.tf
    )

    all_rf_results[[this.tf]] <- imp.df
  }

  #convert interaction_potential list into a matrix
  all_rf_results_matrix <- convertListOfMatricesToMatrix(all_rf_results)

  #stuff for visualization
  lr_markers_this_cluster <- FindMarkers(decipher_seurat_this_cluster,
                                         ident.1 = "case",
                                         ident.2 = "control",
                                         feature = unique(c(all_rf_results_matrix$ligand,all_rf_results_matrix$receptor)),
                                         logfc.threshold = 0,
                                         min.pct = 0,
                                         only.pos = FALSE)

  de_markers_this_cluster <- FindMarkers(decipher_seurat_this_cluster,
                                         ident.1 = "case",
                                         ident.2 = "control",
                                         logfc.threshold = 0.58,
                                         only.pos = FALSE)

  de_markers_this_cluster$gene <- rownames(de_markers_this_cluster)

  #regulon_results_df <- enrichResults(de_markers_this_cluster,significant_regulon_deltas_this_cluster,regulon_this_cluster,enrichr_database)

  #enrichr_results_by_cluster[[this_cluster]] <- regulon_results_df
  regulon_grns_by_cluster[[this_cluster]] <- regulon_this_cluster
  regulon_scores_by_cluster[[this_cluster]] <- regulon_scores_this_cluster
  regulon_deltas_by_cluster[[this_cluster]] <- regulon_deltas_this_cluster
  interaction_potential_by_clusters[[this_cluster]] <- interaction_potentials_matrix_this_cluster
  decipher_scores_by_regulon_and_cluster[[this_cluster]] <- all_rf_results_matrix
  lr_markers_by_cluster[[this_cluster]] <- lr_markers_this_cluster
  de_markers_by_cluster[[this_cluster]] <- de_markers_this_cluster
  significant_regulons_by_cluster[[this_cluster]] <- significant_regulon_deltas_this_cluster
  }
end_time <- Sys.time()
decipher_scores_by_regulon_and_cluster <- lapply(decipher_scores_by_regulon_and_cluster,FUN = "listOfDFsRenameColumn",original_name = "weighted.spearman.cont",new_name = "decipher_score")
decipher_scores_by_cluster <- lapply(decipher_scores_by_regulon_and_cluster,FUN = "calculateScoresByCluster")
decipher_scores_by_cluster <- addListNameToDFElements(decipher_scores_by_cluster,"receiver_cluster")

#save DECIPHER ----
saveRDS(regulon_scores_by_cluster,file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
saveRDS(regulon_grns_by_cluster,file.path(output_data_filepath,"regulon_grns_by_cluster.rds"))
saveRDS(regulon_deltas_by_cluster,file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
saveRDS(significant_regulons_by_cluster,file.path(output_data_filepath,"significant_regulons_by_cluster.rds"))
saveRDS(significant_regulon_markers_by_cluster,file.path(output_data_filepath,"significant_regulon_markers_by_cluster.rds"))
saveRDS(interaction_potential_by_clusters,file.path(output_data_filepath,"interaction_potential_by_clusters.rds"))
saveRDS(interaction_deltas_by_cluster,file.path(output_data_filepath,"interaction_deltas_by_cluster.rds"))
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(lr_markers_by_cluster,file.path(output_data_filepath,"lr_markers_by_cluster.rds"))
saveRDS(de_markers_by_cluster,file.path(output_data_filepath,"de_markers_by_cluster.rds"))
saveRDS(enrichr_results_by_cluster,file.path(output_data_filepath,"enrichr_results_by_cluster.rds"))
saveRDS(feature_statistics,file.path(output_data_filepath,"feature_statistics.rds"))
saveRDS(decipher_seurat_lr,file.path(output_data_filepath,"decipher_seurat_lr.rds"))
saveRDS(L.set,file.path(output_data_filepath,"L_set.rds"))
saveRDS(decipher_scores_by_regulon_and_cluster, file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(decipher_scores_by_cluster,file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))


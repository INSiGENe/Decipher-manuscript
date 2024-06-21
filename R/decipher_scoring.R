#' Extract and Analyze Decipher Results from Random Forest and Interaction Data
#'
#' This function processes the results from a random forest analysis, computing the importance
#' of features and their correlations with transcription factors (TFs). It merges these results
#' with interaction potential data to create a detailed summary of interactions, including ligands,
#' receptors, and their importance metrics. Additionally, it evaluates correlations using both
#' Spearman and Pearson methods and calculates contributions based on these correlations.
#'
#' @param random_forest_results A randomForest object containing the results from a random forest analysis.
#' @param interaction_potentials_matrix_clusters A matrix of interaction potentials between clusters.
#' @param data_this_cluster_downsampled_receptors A matrix or dataframe containing receptor data for a downsampled cluster.
#' @param selected_lr_pairs A dataframe specifying ligand-receptor pairs.
#' @param this_tf The specific transcription factor (TF) being analyzed.
#' @param val.this_tf The value associated with the transcription factor in the analysis.
#'
#' @return A dataframe containing a comprehensive summary of interactions, importance metrics,
#'         correlation scores, and additional computed values related to the transcription factor.
#'
#' @examples
#' # Assuming you have all the necessary data prepared:
#' decipher_results_df <- extractDecipherResults(
#'   random_forest_results = rf_results,
#'   interaction_potentials_matrix_clusters = interaction_matrix,
#'   data_this_cluster_downsampled_receptors = receptor_data,
#'   selected_lr_pairs = lr_pairs,
#'   this_tf = "TF_name",
#'   val.this_tf = 0.5
#' )
#'
#' @importFrom randomForest importance
#' @importFrom stats cor
#' @export
extractDecipherResults <- function(random_forest_results,interaction_potentials_matrix_clusters,data_this_cluster_downsampled_receptors,selected_lr_pairs,this.tf,val.this.tf,tf.merged){

  interaction_mapping_table <-  getInteractionMappingTable(
    receptorMatrix = data_this_cluster_downsampled_receptors,
    ligandSet = selected_lr_pairs
  )

  imp.perm.merged <- importance(random_forest_results,type=1, scale = F)

  spearman.cor <- cor(t(interaction_potentials_matrix_clusters),tf.merged,method = "spearman")
  pearson.cor <- cor(t(interaction_potentials_matrix_clusters),tf.merged,method = "pearson")

  imp <- importance(random_forest_results, scale = F)

  index_match_interaction_mapping_table <- match(rownames(imp),interaction_mapping_table$interaction)


  imp.df <- data.frame(
    interaction = interaction_mapping_table$interaction[index_match_interaction_mapping_table],
    ligand =  interaction_mapping_table$ligand[index_match_interaction_mapping_table],
    receptor =  interaction_mapping_table$receptor[index_match_interaction_mapping_table],
    imp.perm = imp[,1],
    perm.rank = length(imp[,1])-rank(imp[,1]),
    imp.gini = imp[,2],
    gini.rank = length(imp[,2])-rank(imp[,2]),
    gene = rownames(imp),
    regulon = this.tf,
    regulon.val = val.this.tf,
    pearson.cor =  pearson.cor,
    spearman.cor = spearman.cor,
    possible.spearman.cont = spearman.cor*val.this.tf,
    weighted.spearman.cont = imp[,1]*sign(spearman.cor)*val.this.tf
  )

  imp.df <- imp.df[order(imp.df$perm.rank,decreasing=FALSE),]

  return(imp.df)
}


#' Get Random Forest Weights for All Clusters
#'
#' This function calculates random forest weights for all clusters in a Seurat object using significant regulon deltas, regulon scores, interaction potentials, and relevant ligand-receptor features.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param significant_regulon_deltas_all_clusters A list of significant regulon deltas for each cluster.
#' @param regulon_scores_all_clusters A list of regulon scores for each cluster.
#' @param interaction_potentials_matrix_clusters_all_clusters A list of interaction potentials matrices for each cluster.
#' @param L_set_relevant_features_all_clusters A list of relevant ligand-receptor features for each cluster.
#' @param flag.normalize.non.log A logical flag indicating whether to normalize non-log-transformed data.
#'
#' @return A list where each element corresponds to a cluster and contains the random forest importance results for that cluster.
#'
#' @details The function iterates through each unique cluster in the Seurat object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates random forest importance scores for the transcription factors using the `randomForest` package. A progress bar is displayed to show the progress of the computation.
#'
#' @examples
#' \dontrun{
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' significant_regulon_deltas_all_clusters <- getSignificantRegulonsAllClusters(regulon_deltas_all_clusters)
#' regulon_scores_all_clusters <- getRegulonScoresAllClusters(capped_regulons_all_clusters, decipher_seurat)
#' interaction_potentials_matrix_clusters_all_clusters <- getInteractionPotentialMatrixForRepresentativeInteractionsAllClusters(
#'   decipher_seurat, L_set_relevant_features_all_clusters, filtered_interaction_potentials_matrix_all_clusters, cytosig_ligands, TRUE)
#' rf_results_all_clusters <- getRandomForestWeightsAllClusters(
#'   decipher_seurat, significant_regulon_deltas_all_clusters, regulon_scores_all_clusters, interaction_potentials_matrix_clusters_all_clusters, L_set_relevant_features_all_clusters, TRUE)
#' }
#'
#' @import progress
#' @importFrom randomForest randomForest
#' @export
getRandomForestWeightsAllClusters <- function(decipher_seurat, significant_regulon_deltas_all_clusters, regulon_scores_all_clusters, interaction_potentials_matrix_clusters_all_clusters, L_set_relevant_features_all_clusters, flag.normalize.non.log) {
  rf_results_all_clusters <- list()

  # Calculate total number of tasks
  total_tasks <- sum(sapply(significant_regulon_deltas_all_clusters, function(x) length(x$name)))

  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "  [:bar] :current/:total (:percent) :elapsedfull",
    total = total_tasks, clear = FALSE, width = 60
  )

  for(this_cluster in names(significant_regulon_deltas_all_clusters)){
    all_rf_results <- list()
    significant_regulon_deltas_this_cluster <- significant_regulon_deltas_all_clusters[[this_cluster]]
    regulon_scores_this_cluster <- regulon_scores_all_clusters[[this_cluster]]
    interaction_potentials_matrix_clusters <- interaction_potentials_matrix_clusters_all_clusters[[this_cluster]]

    # main object
    decipher_seurat_this_cluster <- subset(decipher_seurat, subset = cluster == this_cluster)
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition
    data_this_cluster <- decipher_seurat_this_cluster@assays$RNA@data

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    data_this_cluster_receptors <- data_this_cluster[which(rownames(data_this_cluster) %in% unique(L_set_relevant_features_all_clusters[[this_cluster]]$receptor)),]
    for(this.tf in significant_regulon_deltas_this_cluster$name){
      ind.this.tf <- which(significant_regulon_deltas_this_cluster$name == this.tf)
      val.this.tf <- significant_regulon_deltas_this_cluster$deltaPagoda[ind.this.tf]
      tf.merged <- regulon_scores_this_cluster[this.tf, colnames(interaction_potentials_matrix_clusters)]

      set.seed(123)
      rf <- randomForest::randomForest(
        x = t(interaction_potentials_matrix_clusters),
        y=tf.merged,
        ntree = 100,
        importance=T)

      imp.df <- extractDecipherResults(
        random_forest_results = rf,
        interaction_potentials_matrix_clusters,
        data_this_cluster_receptors,
        selected_lr_pairs = L_set_relevant_features_all_clusters[[this_cluster]],
        this.tf,
        val.this.tf,
        tf.merged
      )

      all_rf_results[[this.tf]] <- imp.df

      # Update progress bar
      pb$tick()
    }
    all_rf_results_matrix <- convertListOfMatricesToMatrix(all_rf_results)
    rf_results_all_clusters[[this_cluster]] <- all_rf_results_matrix

  }
  return(rf_results_all_clusters)
}

#' Find LR Marker Genes for Ligand-Receptor Pairs in All Clusters
#'
#' This function identifies marker genes for ligand-receptor pairs in each cluster of a Seurat object based on random forest results.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param rf_results_all_clusters A list of random forest results for each cluster.
#' @param flag.normalize.non.log A logical flag indicating whether to normalize non-log-transformed data.
#'
#' @return A list where each element corresponds to a cluster and contains the marker genes for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and identifies marker genes for ligand-receptor pairs using the `FindMarkers` function.
#'
#' @examples
#' \dontrun{
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' rf_results_all_clusters <- getRandomForestWeightsAllClusters(
#'   decipher_seurat, significant_regulon_deltas_all_clusters, regulon_scores_all_clusters,
#'   interaction_potentials_matrix_clusters_all_clusters, L_set_relevant_features_all_clusters, TRUE)
#' lr_markers_all_clusters <- FindLRMarkersAllClusters(decipher_seurat, rf_results_all_clusters, TRUE)
#' }
#'
#' @importFrom Seurat NormalizeData FindMarkers
#' @export
FindLRMarkersAllClusters <- function(decipher_seurat, rf_results_all_clusters, flag.normalize.non.log) {
  lr_markers_all_clusters <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){

    # main object
    decipher_seurat_this_cluster <- subset(decipher_seurat, subset = cluster == this_cluster)
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition
    data_this_cluster <- decipher_seurat_this_cluster@assays$RNA@data

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster$condition

    all_rf_results_matrix <- rf_results_all_clusters[[this_cluster]]
    lr_markers_all_clusters[[this_cluster]] <- FindMarkers(
      object = decipher_seurat_this_cluster,
      ident.1 = "case",
      ident.2 = "control",
      features = unique(c(all_rf_results_matrix$ligand, all_rf_results_matrix$receptor)),
      logfc.threshold = 0,
      min.pct = 0,
      only.pos = FALSE
    )
  }
  return(lr_markers_all_clusters)
}

#' Find Differentially Expressed Markers for All Clusters
#'
#' This function identifies differentially expressed markers for each cluster in a Seurat object.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param flag.normalize.non.log A logical flag indicating whether to normalize non-log-transformed data.
#'
#' @return A list where each element corresponds to a cluster and contains the differentially expressed markers for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and identifies differentially expressed markers using the `FindMarkers` function. The results for each cluster are stored in a list.
#'
#' @examples
#' \dontrun{
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' markers_all_clusters <- FindMarkersAllClusters(decipher_seurat, TRUE)
#' }
#'
#' @importFrom Seurat NormalizeData FindMarkers
#' @export
FindMarkersAllClusters <- function(decipher_seurat, flag.normalize.non.log) {
  markers_all_clusters <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){

    # main object
    decipher_seurat_this_cluster <- subset(decipher_seurat, subset = cluster == this_cluster)
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition
    data_this_cluster <- decipher_seurat_this_cluster@assays$RNA@data

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster$condition

    de_markers_this_cluster <- FindMarkers(
      object = decipher_seurat_this_cluster,
      ident.1 = "case",
      ident.2 = "control",
      logfc.threshold = 0.58,
      only.pos = FALSE
    )

    de_markers_this_cluster$gene <- rownames(de_markers_this_cluster)

    markers_all_clusters[[this_cluster]] <- de_markers_this_cluster
  }
  return(markers_all_clusters)
}




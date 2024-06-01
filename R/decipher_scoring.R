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
extractDecipherResults <- function(random_forest_results,interaction_potentials_matrix_clusters,data_this_cluster_downsampled_receptors,selected_lr_pairs,this.tf,val.this.tf){

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

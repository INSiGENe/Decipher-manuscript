
#' Generate Interaction Mapping Table for Ligands and Receptors
#'
#' This function creates a mapping table for ligands and receptors based on
#' a receptor expression matrix and a ligand-receptor set. For each receptor
#' in the receptor matrix, it identifies the associated ligands from the
#' ligand-receptor set and constructs a table with the interaction pairs.
#'
#' @param receptorMatrix A matrix with receptors as row names.
#' @param ligandSet A data frame with columns 'ligand' and 'receptor', representing
#' ligand-receptor pairs.
#' @return A data frame where each row represents an interaction pair, with columns
#' for the interaction name (concatenation of ligand and receptor names), ligand name,
#' and receptor name.
#' @export
#'
#' @examples
#' # receptorMatrix is a matrix with receptor genes as row names
#' # ligandSet is a data frame with ligand-receptor pairs
#' interactionTable <- getInteractionMappingTable(receptorMatrix, ligandSet)
getInteractionMappingTable <- function(receptorMatrix, ligandSet) {
  # Validate arguments
  # if (!is.matrix(receptorMatrix)) {
  #   stop("receptorMatrix must be a matrix")
  # }
  # if (!is.data.frame(ligandSet)) {
  #   stop("ligandSet must be a data frame")
  # }

  # Extract highly expressed receptors
  highlyExpressedReceptors <- rownames(receptorMatrix)

  # Initialize empty data frame for interaction mapping
  interactionMappingTable <- data.frame(interaction = character(), ligand = character(), receptor = character())

  # Iterate over each receptor
  for (receptor in highlyExpressedReceptors) {
    # Find associated ligands
    associatedLigands <- ligandSet$ligand[ligandSet$receptor == receptor]

    # Construct interaction names
    interactions <- paste(associatedLigands, receptor, sep = "-")

    # Create a temporary data frame for current receptor and its ligands
    tempMappingTable <- data.frame(interaction = interactions, ligand = associatedLigands, receptor = receptor)

    # Combine with the main mapping table
    interactionMappingTable <- rbind(interactionMappingTable, tempMappingTable)
  }

  interactionMappingTable
}


#' Calculate Interaction Potential Matrix
#'
#' This function computes an interaction potential matrix for ligand-receptor pairs
#' based on the receptor expression matrix, a condition vector, mean ligand values,
#' and a ligand-receptor set (LRSet). It calculates the interaction potential for each
#' ligand-receptor pair under different conditions (e.g., 'case' and 'control').
#'
#' @param receptorMatrix A matrix with receptors as row names and samples as columns.
#' @param conditionVector A vector indicating the condition (e.g., 'case' or 'control')
#' for each sample in `receptorMatrix`.
#' @param ligandMeans A data frame with mean values of ligands under different conditions.
#' @param LRSet A data frame representing ligand-receptor pairs and their interactions.
#' @return An interaction potential matrix where rows represent ligand-receptor interactions
#' and columns represent samples. Each element in the matrix is the calculated interaction
#' potential for a ligand-receptor pair under a specific condition.
#' @export
#'
#' @examples
#' # receptorMatrix is a matrix with receptor genes as row names
#' # conditionVector is a vector indicating the condition of each sample
#' # ligandMeans is a data frame with mean ligand values under different conditions
#' # LRSet is a data frame with ligand-receptor pairs
#' interactionMatrix <- calculateInteractionMatrix(receptorMatrix, conditionVector, ligandMeans, LRSet)
calculateInteractionMatrix <- function(receptorMatrix, conditionVector, ligandMeans, LRSet) {
  # Validate arguments
  # if (!is.matrix(receptorMatrix)) {
  #   stop("receptorMatrix must be a matrix")
  # }
  if (!is.vector(conditionVector) || length(conditionVector) != ncol(receptorMatrix)) {
    stop("conditionVector must be a vector of the same length as the number of columns in receptorMatrix")
  }
  if (!is.data.frame(ligandMeans)) {
    stop("ligandMeans must be a data frame")
  }
  if (!is.data.frame(LRSet)) {
    stop("LRSet must be a data frame")
  }

  # Initialize interaction potential matrix
  interactionPotentialMatrix <- matrix(0, nrow = nrow(LRSet), ncol = ncol(receptorMatrix),
                                       dimnames = list(LRSet$interaction, colnames(receptorMatrix)))

  # Identify case and control indices
  indexCase <- which(conditionVector == "case")
  indexControl <- which(conditionVector == "control")

  # Iterate through ligand-receptor pairs
  for (thisEntry in seq_len(nrow(LRSet))) {
    receptor <- LRSet$receptor[thisEntry]
    ligand <- LRSet$ligand[thisEntry]
    interaction <- LRSet$interaction[thisEntry]

    # Extract receptor vectors for case and control
    receptorVectorCase <- receptorMatrix[receptor, indexCase]
    receptorVectorControl <- receptorMatrix[receptor, indexControl]

    # Calculate mean ligand values
    meanLigandValueCase <- ligandMeans[ligand, "case"]
    meanLigandValueControl <- ligandMeans[ligand, "control"]

    # Calculate interaction potential
    interactionPotentialCase <- sqrt(receptorVectorCase * meanLigandValueCase)
    interactionPotentialControl <- sqrt(receptorVectorControl * meanLigandValueControl)

    # Assign values to interaction matrix
    interactionPotentialMatrix[interaction, c(indexCase, indexControl)] <- c(interactionPotentialCase, interactionPotentialControl)
  }

  interactionPotentialMatrix
}


#' Extract One-to-One Interactions from a Mapping Table
#'
#' This function identifies one-to-one interactions in a mapping table, where each
#' receptor is associated with exactly one ligand. It returns the interactions
#' that have a unique receptor-ligand pair.
#'
#' @param mappingTable A data frame representing a mapping table with at least
#' two columns: 'receptor' and 'interaction'.
#' @return A vector of interaction names that represent one-to-one ligand-receptor
#' pairs from the mapping table.
#' @export
#'
#' @examples
#' # mappingTable is a data frame with columns 'interaction' and 'receptor'
#' oneToOneInteractions <- getOneToOneInteractions(mappingTable)
getOneToOneInteractions <- function(mappingTable) {
  # Validate arguments
  if (!is.data.frame(mappingTable)) {
    stop("mappingTable must be a data frame")
  }

  # Count occurrences of each receptor
  receptorCounts <- table(mappingTable$receptor)

  # Identify receptors that appear exactly once
  oneToOneReceptors <- names(receptorCounts)[receptorCounts == 1]

  # Subset the mapping table to get one-to-one interactions
  oneToOneInteractions <- mappingTable$interaction[mappingTable$receptor %in% oneToOneReceptors]

  oneToOneInteractions
}


#' Extract Many-to-One Interactions from a Mapping Table
#'
#' This function identifies many-to-one interactions in a mapping table, where each
#' receptor is associated with multiple ligands. It returns the interactions
#' that have a receptor linked to more than one ligand.
#'
#' @param mappingTable A data frame representing a mapping table with at least
#' two columns: 'receptor' and 'interaction'.
#' @return A vector of interaction names that represent many-to-one ligand-receptor
#' pairs from the mapping table.
#' @export
#'
#' @examples
#' # mappingTable is a data frame with columns 'interaction' and 'receptor'
#' manyToOneInteractions <- getManyToOneInteractions(mappingTable)
getManyToOneInteractions <- function(mappingTable) {
  # Validate arguments
  if (!is.data.frame(mappingTable)) {
    stop("mappingTable must be a data frame")
  }

  # Count occurrences of each receptor
  receptorCounts <- table(mappingTable$receptor)

  # Identify receptors that are associated with more than one ligand
  manyToOneReceptors <- names(receptorCounts)[receptorCounts > 1]

  # Filter the mapping table for many-to-one interactions
  manyToOneInteractions <- mappingTable$interaction[mappingTable$receptor %in% manyToOneReceptors]

  manyToOneInteractions
}


#' Cluster Many-to-One Receptors Based on Interaction Potentials
#'
#' This function performs hierarchical clustering on many-to-one receptors based on
#' their interaction potentials. It computes a correlation matrix for the interaction
#' potentials and then applies hierarchical clustering to identify clusters of
#' receptors.
#'
#' @param interactionPotentialsMatrixMTO A matrix of interaction potentials for
#' many-to-one (MTO) receptor-ligand pairs.
#' @param interactionMappingTable A data frame representing a mapping table with
#' columns for receptor-ligand interactions.
#' @param pctMTOReceptors A numeric value between 0 and 1 indicating the percentage
#' of MTO receptors used to determine the number of clusters.
#' @param correlationMethod The method of correlation computation (default is "pearson").
#' @param clusteringMethod The method used for hierarchical clustering (default is "complete").
#' @return A vector indicating the cluster assignment for each receptor.
#' @export
#'
#' @examples
#' # interactionPotentialsMatrixMTO is a matrix of interaction potentials
#' # interactionMappingTable is a data frame with interaction mapping
#' # pctMTOReceptors is the percentage of receptors to consider for clustering
#' clusters <- getCorrelationClusters(interactionPotentialsMatrixMTO, interactionMappingTable, 0.5)
getCorrelationClusters <- function(interactionPotentialsMatrixMTO, interactionMappingTable, pctMTOReceptors, correlationMethod = "pearson", clusteringMethod = "complete") {
  # Validate arguments
  if (!is.matrix(interactionPotentialsMatrixMTO)) {
    stop("interactionPotentialsMatrixMTO must be a matrix")
  }
  if (!is.data.frame(interactionMappingTable)) {
    stop("interactionMappingTable must be a data frame")
  }

  # Identify MTO receptors
  mtoReceptors <- unique(interactionMappingTable$receptor[interactionMappingTable$interaction %in% rownames(interactionPotentialsMatrixMTO)])

  # Compute correlation matrix
  correlationMatrix <- cor(t(interactionPotentialsMatrixMTO), method = correlationMethod)
  correlationMatrix[is.na(correlationMatrix)] <- 0

  # Perform hierarchical clustering
  hclustResults <- hclust(as.dist(1 - abs(correlationMatrix)), method = clusteringMethod)

  # Determine the number of clusters
  numberOfClusters <- floor(pctMTOReceptors * length(mtoReceptors))

  # Cut the dendrogram into clusters
  clusters <- cutree(hclustResults, k = numberOfClusters)

  clusters
}


#' Select Representative Interactions for Many-to-One Clusters
#'
#' This function selects a representative interaction for each cluster of many-to-one
#' (MTO) receptor-ligand pairs. If provided, it prioritizes certain ligands. The function
#' also issues a warning if a cluster contains multiple unique receptors.
#'
#' @param mtoInteractionsClusters A vector representing the cluster assignments for each
#' many-to-one interaction.
#' @param interactionMappingTable A data frame representing a mapping table with columns
#' 'interaction', 'receptor', and 'ligand'.
#' @param prioritizedBenchmarkingLigands An optional vector of ligand names to prioritize
#' when selecting representative interactions. If not NULL, interactions with these ligands
#' are preferred.
#' @return A data frame with two columns: 'cluster', indicating the cluster number, and
#' 'interaction', the name of the representative interaction for that cluster.
#' @export
#'
#' @examples
#' # mtoInteractionsClusters is a vector of cluster assignments
#' # interactionMappingTable is a data frame with interaction mappings
#' # prioritizedBenchmarkingLigands is a vector of prioritized ligand names
#' representativeInteractions <- getRepresentativeInteractionsForMTOClusters(mtoInteractionsClusters, interactionMappingTable, prioritizedBenchmarkingLigands)
getRepresentativeInteractionsForMTOClusters <- function(mtoInteractionsClusters, interactionMappingTable, prioritizedBenchmarkingLigands = NULL) {
  # Validate arguments
  if (!is.vector(mtoInteractionsClusters)) {
    stop("mtoInteractionsClusters must be a vector")
  }
  if (!is.data.frame(interactionMappingTable)) {
    stop("interactionMappingTable must be a data frame")
  }

  # Initialize data frame for representative interactions
  representativeInteractions <- data.frame(cluster = integer(), interaction = character())

  # Iterate over each cluster
  for (cluster in unique(mtoInteractionsClusters)) {
    interactionsInCluster <- names(mtoInteractionsClusters)[mtoInteractionsClusters == cluster]

    # Get interactions from mapping table
    interactionIndices <- match(interactionsInCluster, interactionMappingTable$interaction)
    receptorsInCluster <- interactionMappingTable$receptor[interactionIndices]
    ligandsInCluster <- interactionMappingTable$ligand[interactionIndices]

    # Check if a single receptor is unique to this cluster
    if (length(unique(receptorsInCluster)) > 1) {
      warning("Cluster ", cluster, " has multiple receptors. Selecting one at random.")
    }

    # Prioritize benchmarking ligands if provided
    if (!is.null(prioritizedBenchmarkingLigands)) {
      prioritizedIndices <- which(ligandsInCluster %in% prioritizedBenchmarkingLigands)
      if (length(prioritizedIndices) > 0) {
        interactionsInCluster <- interactionsInCluster[prioritizedIndices]
      }
    }

    # Select a representative interaction
    representativeInteraction <- sample(interactionsInCluster, size = 1)

    # Add entry to the result
    representativeInteractions <- rbind(representativeInteractions, data.frame(cluster = cluster, interaction = representativeInteraction))
  }

  representativeInteractions
}


#' Calculate Interaction Deltas for Cluster-Specific Interaction Potentials
#'
#' This function creates a new Seurat object using a matrix of interaction potentials specific to a cluster.
#' It then performs differential expression analysis to calculate interaction deltas between conditions
#' specified in the metadata. This analysis helps identify significant interactions in terms of adjusted p-values.
#'
#' @param interaction_potentials_matrix_this_cluster A matrix of interaction potentials for a specific cluster.
#' @param decipher_seurat_lr A Seurat object from which metadata will be used to annotate the newly created Seurat object.
#'
#' @return A data frame containing the interaction deltas, including log-fold changes, p-values,
#'         adjusted p-values, and gene names, filtered by significance based on adjusted p-value.
#'
#' @examples
#' # Assuming 'interaction_matrix' is a matrix and 'decipher_seurat' is a Seurat object:
#' interaction_deltas <- calculateInteractionDeltas(
#'   interaction_potentials_matrix_this_cluster = interaction_matrix,
#'   decipher_seurat_lr = decipher_seurat
#' )
#'
#' @importFrom Seurat CreateSeuratObject FindMarkers
#' @importFrom dplyr filter
#' @export
calculateInteractionDeltas <- function(interaction_potentials_matrix_this_cluster,decipher_seurat_lr){
  new_seurat <- Seurat::CreateSeuratObject(counts = interaction_potentials_matrix_this_cluster,meta.data = decipher_seurat_lr@meta.data[colnames(interaction_potentials_matrix_this_cluster),])
  SeuratObject::Idents(new_seurat) <- new_seurat$condition
  # Perform differential expression analysis to find markers between specified conditions
  interaction_deltas <- FindMarkers(new_seurat,ident.1 = "case",logfc.threshold = 0.1)
  interaction_deltas <- interaction_deltas %>%
    dplyr::filter(p_val_adj < 0.01)
  interaction_deltas$name = rownames(interaction_deltas)
  return(interaction_deltas)
}


#' Compute Interaction Potentials Matrix for a Specific Cluster
#'
#' This function calculates the interaction potentials between ligands and receptors
#' within a specific cluster. It separates the case and control conditions in the
#' main Seurat object, computes the mean expression of selected ligands for both
#' conditions, and then uses these means along with the expression of receptors
#' in a downsampled cluster-specific Seurat object to compute an interaction
#' potentials matrix.
#'
#' @param seurat_obj A Seurat object containing the full dataset.
#' @param seurat_obj_this_cluster_ds A downsampled Seurat object for the specific cluster.
#' @param selected_lr_pairs A dataframe specifying ligand-receptor pairs to consider,
#'        which must have columns 'ligand' and 'receptor'.
#'
#' @return An interaction potentials matrix where rows represent receptors and
#'         columns correspond to interaction potentials derived from case-control comparisons
#'         for the ligands.
#'
#' @examples
#' # Assuming 'full_seurat' is the full Seurat object, 'cluster_seurat_ds' is the downsampled Seurat object for a cluster,
#' # and 'lr_pairs' is a dataframe of selected ligand-receptor pairs:
#' interaction_matrix <- getInteractionPotentialsMatrixThisCluster(
#'   seurat_obj = full_seurat,
#'   seurat_obj_this_cluster_ds = cluster_seurat_ds,
#'   selected_lr_pairs = lr_pairs
#' )
#'
#' @importFrom Matrix rowMeans
#' @export
getInteractionPotentialsMatrixThisCluster <- function(seurat_obj,seurat_obj_this_cluster_ds,selected_lr_pairs){

  # Extract indices for case condition and data for case and control
  ind_case <- which(seurat_obj$condition == "case")
  data_seurat_obj_case <- seurat_obj@assays$RNA@data[,ind_case]
  data_seurat_obj_control <- seurat_obj@assays$RNA@data[,-ind_case]

  # Extract receptor data from the downsampled Seurat object for this cluster
  data_seurat_obj_this_cluster_ds <- seurat_obj_this_cluster_ds@assays$RNA@data
  data_seurat_obj_this_cluster_ds_receptors <- data_seurat_obj_this_cluster_ds[which(rownames(data_seurat_obj_this_cluster_ds) %in% unique(selected_lr_pairs$receptor)),]

  # Compute mean expression levels of ligands in case and control
  unique_ligands <- unique(selected_lr_pairs$ligand)
  ligand_means <- data.frame(
    case = Matrix::rowMeans(data_seurat_obj_case[unique_ligands,]),
    control = Matrix::rowMeans(data_seurat_obj_control[unique_ligands,]),
    row.names = unique_ligands
  )

  # Calculate interaction matrix
  interaction_potentials_matrix_this_cluster <- calculateInteractionMatrix(
    receptorMatrix = data_seurat_obj_this_cluster_ds_receptors,
    conditionVector =  seurat_obj_this_cluster_ds$condition,
    ligandMeans = ligand_means,
    LRSet = selected_lr_pairs
  )

  # Remove rows with no interaction information
  ind_no_information <- which(rowSums(interaction_potentials_matrix_this_cluster) == 0)
  if(length(ind_no_information) > 0){
    interaction_potentials_matrix_this_cluster <- interaction_potentials_matrix_this_cluster[-ind_no_information,]
  }

  return(interaction_potentials_matrix_this_cluster)
}


#' Generate Interaction Potential Matrix for Representative Interactions
#'
#' This function analyzes interaction potentials for a specific cluster and categorizes
#' interactions between receptors and ligands into one-to-one and many-to-one types.
#' It then identifies representative interactions for the many-to-one category based on
#' clustering of interaction potentials. This helps in reducing complexity and focusing
#' on representative interaction dynamics within the cluster.
#'
#' @param data_this_cluster_downsampled_receptors Matrix of receptor expression data
#'        for a downsampled cluster.
#' @param selected_lr_pairs Data frame of selected ligand-receptor pairs.
#' @param interaction_potentials_matrix_this_cluster Matrix of interaction potentials for the cluster.
#' @param cytosig_ligands Vector of prioritized benchmarking ligands used to determine
#'        representative interactions in many-to-one clusters.
#'
#' @return A matrix that combines the interaction potentials for one-to-one interactions
#'         and representative interactions for many-to-one clusters.
#'
#' @examples
#' # Assuming the required matrices and data frames are predefined:
#' interaction_matrix <- getInteractionPotentialMatrixForRepresentativeInteractions(
#'   data_this_cluster_downsampled_receptors = receptor_data,
#'   selected_lr_pairs = lr_pairs,
#'   interaction_potentials_matrix_this_cluster = interaction_matrix,
#'   cytosig_ligands = benchmarking_ligands
#' )
#'
#' @export
getInteractionPotentialMatrixForRepresentativeInteractions <- function(data_this_cluster_downsampled_receptors,selected_lr_pairs,interaction_potentials_matrix_this_cluster,cytosig_ligands){

  interaction_mapping_table <-  getInteractionMappingTable(
    receptorMatrix = data_this_cluster_downsampled_receptors,
    ligandSet = selected_lr_pairs
  )

  #split matrix into interactions comprised of receptors with a unique ligand (one-to-one), and interactions of receptors with multiple ligands (many-to-one)
  one_to_one_interactions <- intersect(getOneToOneInteractions(interaction_mapping_table),rownames(interaction_potentials_matrix_this_cluster))
  many_to_one_interactions <- intersect(getManyToOneInteractions(interaction_mapping_table),rownames(interaction_potentials_matrix_this_cluster))

  interaction_potentials_matrix_this_cluster_oto <- interaction_potentials_matrix_this_cluster[one_to_one_interactions,]
  interaction_potentials_matrix_this_cluster_mto <- interaction_potentials_matrix_this_cluster[many_to_one_interactions,]


  ## correlation clusters for many-to-one interactions ----
  mto_interactions_clusters <- getCorrelationClusters(
    interactionPotentialsMatrixMTO = interaction_potentials_matrix_this_cluster_mto,
    interactionMappingTable = interaction_mapping_table,
    pctMTOReceptors = 1.15,
    correlationMethod = "spearman",
    clusteringMethod = "complete")

  ## representative interaction for each cluster ----
  representative_interactions_mto <- getRepresentativeInteractionsForMTOClusters(
    mtoInteractionsClusters = mto_interactions_clusters,
    interactionMappingTable = interaction_mapping_table,
    prioritizedBenchmarkingLigands = cytosig_ligands
  )

  ## cluster-based matrix from random forest -----
  interaction_potentials_matrix_this_cluster_mto_representative <- interaction_potentials_matrix_this_cluster[representative_interactions_mto$interaction,]
  interaction_potentials_matrix_clusters <- rbind(interaction_potentials_matrix_this_cluster_oto,interaction_potentials_matrix_this_cluster_mto_representative)

}


#' Get Interaction Potentials Matrix for All Clusters
#'
#' This function calculates the interaction potentials matrix for each cluster in a Seurat object.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param decipher_seurat_this_cluster A Seurat object subset for a specific cluster.
#' @param L_set_relevant_features_all_clusters A list of relevant ligand-receptor features for each cluster.
#' @param flag.normalize.non.log A logical flag indicating whether to normalize non-log-transformed data.
#'
#' @return A list where each element corresponds to a cluster and contains the interaction potentials matrix for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates the interaction potentials matrix for each cluster using the `getInteractionPotentialsMatrixThisCluster` function.
#'
#' @examples
#' \dontrun{
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(L.set, expressed_ligands, expressed_receptors_all_clusters)
#' interaction_potentials_matrix <- getInteractionPotentialsMatrixAllClusters(
#'   decipher_seurat, decipher_seurat_this_cluster, L_set_relevant_features_all_clusters, TRUE)
#' }
#'
#' @export
getInteractionPotentialsMatrixAllClusters <- function(decipher_seurat, decipher_seurat_this_cluster, L_set_relevant_features_all_clusters, flag.normalize.non.log) {
  interaction_potentials_matrix_all_clusters <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){
    # main object
    decipher_seurat_this_cluster <- subset(decipher_seurat, subset = cluster == this_cluster)
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster$condition

    interaction_potentials_matrix_this_cluster <- getInteractionPotentialsMatrixThisCluster(
      seurat_obj = decipher_seurat,
      seurat_obj_this_cluster_ds = decipher_seurat_this_cluster,
      selected_lr_pairs = L_set_relevant_features_all_clusters[[this_cluster]]
    )

    interaction_potentials_matrix_all_clusters[[this_cluster]] <- interaction_potentials_matrix_this_cluster
  }

  return(interaction_potentials_matrix_all_clusters)
}



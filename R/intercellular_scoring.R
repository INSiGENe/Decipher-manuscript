
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


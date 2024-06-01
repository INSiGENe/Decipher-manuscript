#' Generate Pseudobulk Matrices from Seurat Object
#'
#' Iterates over clusters and conditions in a Seurat object to generate pseudobulk matrices.
#' It handles minimum and maximum cell count constraints for pseudobulk generation.
#'
#' @param seuratObj A Seurat object.
#' @param paramMinMetaCells The minimum number of meta cells to generate.
#' @param paramMaxMetaCells The maximum number of meta cells to generate.
#' @param paramMaxScCells The maximum number of single cells to consider within each condition.
#' @param paramK Parameter 'k' used in calculatePseudoBulkCell function.
#' @return A list of pseudobulk matrices organized by cluster and condition.
#' @importFrom dplyr filter rename sample_n
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @export
generateMetaCellMatrices <- function(seuratObj, paramMinMetaCells = 100, paramMaxMetaCells = 600, paramMaxScCells, paramK) {
  B_matrices <- list()

  for (this_cluster in unique(seuratObj$cluster)) {
    cat("Calculating pseudobulk matrices for cluster:", this_cluster, "\n")

    seuratObjectCluster <- subset(seuratObj, subset = cluster == this_cluster)
    minCellCount <- min(table(seuratObjectCluster$condition))
    minCellCount <- floor(minCellCount/(paramK+1))

    if (minCellCount < paramMinMetaCells) {
      next
    }

    minCellCount <- min(minCellCount, paramMaxMetaCells)
    cat("Number of pseudobulk cells:", minCellCount, "\n")

    for (this_condition in unique(seuratObjectCluster$condition)) {
      cat("Calculating pseudobulk matrices for condition:", this_condition, "\n")

      conditionData <- seuratObjectCluster@meta.data %>%
        tibble::rownames_to_column(var = "cell") %>%
        dplyr::filter(condition == this_condition)

      if (nrow(conditionData) > paramMaxScCells) {
        cellsToKeep <- conditionData %>% dplyr::select(cell) %>% dplyr::sample_n(size = paramMaxScCells) %>% unlist(use.names = FALSE)
      } else {
        cellsToKeep <- conditionData$cell
      }

      seuratObjectCondition <- subset(seuratObjectCluster, cells = cellsToKeep)

      B_matrix <- calculatePseudoBulkCell(
        seuratObject = seuratObjectCondition,
        numNearestNeighbors = paramK,
        numMetaCells = minCellCount
      )

      B_matrices[[this_cluster]][[this_condition]] <- B_matrix
    }
  }

  return(B_matrices)
}


#' Calculate Pseudo Bulk Cell Data
#'
#' This function calculates pseudo bulk cell data from a Seurat object.
#' It first retrieves the RNA counts matrix, then computes a distance matrix,
#' and finally calculates the pseudo bulk matrix using the specified number
#' of nearest neighbors and meta cells.
#'
#' @param seuratObject A Seurat object containing single-cell expression data.
#' @param numNearestNeighbors Integer, the number of nearest neighbors to consider.
#' @param numMetaCells Integer, the number of meta cells to calculate.
#' @return A matrix representing pseudo bulk data calculated from the Seurat object.
#'
#' @examples
#' # seuratObject is a pre-existing Seurat object
#' pseudoBulkMatrix <- calculatePseudoBulkCell(seuratObject, 5, 10)
calculatePseudoBulkCell <- function(seuratObject, numNearestNeighbors, numMetaCells) {
  # Validate arguments
  if (!inherits(seuratObject, "Seurat")) {
    stop("seuratObject must be a Seurat object")
  }

  # Obtain RNA counts matrix
  rnaCountsMatrix <- Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "counts")
  rnaCountsMatrix <- as.matrix(rnaCountsMatrix)
  # Calculate distance matrix
  distanceMatrix <- calculateDistVarF(seuratObject = seuratObject)

  # Calculate pseudo bulk matrix
  pseudoBulkMatrix <- calculatePseudoBulkMatrix(
    rnaCountsMatrix = rnaCountsMatrix,
    distanceMatrix = distanceMatrix,
    numNearestNeighbors = numNearestNeighbors,
    numMetaCells = numMetaCells
  )

  # Return the result
  pseudoBulkMatrix
}

#' Calculate Distance Matrix from Variable Features of a Seurat Object
#'
#' This function calculates a distance matrix based on the variable features
#' of a Seurat object. It identifies variable features in the RNA assay of the
#' Seurat object and computes a distance matrix using these features.
#'
#' @param seuratObject A Seurat object containing single-cell RNA sequencing data.
#' @return A distance matrix computed from the variable features of the RNA assay
#' in the Seurat object.
#'
#' @importFrom Seurat FindVariableFeatures VariableFeatures GetAssayData
#'
#' @examples
#' # seuratObject is a pre-existing Seurat object
#' distanceMatrix <- calculateDistVarF(seuratObject)
calculateDistVarF <- function(seuratObject) {
  # Validate arguments
  if (!inherits(seuratObject, "Seurat")) {
    stop("seuratObject must be a Seurat object")
  }

  # Find and retrieve variable features
  seuratObject <- Seurat::FindVariableFeatures(seuratObject)
  variableFeatures <- Seurat::VariableFeatures(seuratObject, assay = "RNA")

  # Retrieve data matrix for RNA assay
  dataMatrix <-  Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "data")
  dataMatrix <- as.matrix(dataMatrix)
  # Calculate and return the distance matrix based on variable features
  CorDist(dataMatrix[variableFeatures, ])
}

#' Calculate Distance Matrix Using Correlation
#'
#' This function computes a distance matrix based on the correlation of the input matrix.
#' Optionally, it can first perform a principal component analysis (PCA) and use the
#' specified number of PCA features for the distance calculation.
#'
#' @param matrix A numeric matrix on which the distance calculation is based.
#' @param pcaFeats An optional integer specifying the number of principal components
#' to retain for distance calculation. If NULL, PCA is not performed.
#' @param corMethod A character string indicating the correlation method to be used.
#' Can be one of 'pearson', 'kendall', or 'spearman'.
#' @return A distance matrix, with each element representing the distance
#' calculated based on the correlation between columns of the input matrix.
#'
#' @importFrom stats as.dist
#'
#' @examples
#' dataMatrix <- matrix(rnorm(100), ncol = 10)
#' distanceMatrix <- CorDist(dataMatrix)
#' # With PCA transformation
#' pcaDistanceMatrix <- CorDist(dataMatrix, pcaFeats = 5)
CorDist <- function(matrix, pcaFeats = NULL, corMethod = 'spearman') {
  # Validate arguments
  if (!is.null(pcaFeats) && (!is.numeric(pcaFeats) || pcaFeats <= 0)) {
    stop("pcaFeats must be a positive integer")
  }

  # Perform PCA if pcaFeats is specified
  if (!is.null(pcaFeats)) {
    pcaResult <- stats::prcomp(t(matrix))
    matrix <- t(pcaResult$x[, 1:pcaFeats])
  }

  # Calculate and return distance matrix
  result <- as.dist(sqrt(1 - stats::cor(matrix, method = corMethod)), upper = TRUE, diag = TRUE)
  result <- as.matrix(result)
  return(result)
}


#' Calculate Pseudo Bulk Matrix
#'
#' This function creates a pseudo bulk matrix from a given RNA counts matrix
#' and a distance matrix. It samples a subset of columns (representing cells)
#' from the RNA counts matrix and calculates the sum of counts for each
#' cell and its nearest neighbors.
#'
#' @param rnaCountsMatrix A matrix containing RNA counts data, with rows
#' representing genes and columns representing cells.
#' @param distanceMatrix A precomputed distance matrix where each element
#' represents the distance between cells.
#' @param numNearestNeighbors Integer, the number of nearest neighbors to include
#' for each sampled cell.
#' @param numMetaCells Integer, the number of cells (columns from rnaCountsMatrix)
#' to sample for pseudo bulk calculation.
#' @return A matrix where each column represents a pseudo bulk sample, containing
#' the sum of counts for each gene across a cell and its nearest neighbors.
#'
#' @examples
#' # rnaCountsMatrix and distanceMatrix are pre-existing matrices
#' pseudoBulkMatrix <- calculatePseudoBulkMatrix(rnaCountsMatrix, distanceMatrix, 5, 10)
calculatePseudoBulkMatrix <- function(rnaCountsMatrix, distanceMatrix, numNearestNeighbors, numMetaCells) {
  # Validate arguments
  if (!is.matrix(rnaCountsMatrix)) {
    stop("rnaCountsMatrix must be a matrix")
  }
  if (!is.matrix(distanceMatrix)) {
    stop("distanceMatrix must be a matrix")
  }

  # Sampling column names
  subSamples <- sample(colnames(rnaCountsMatrix), numMetaCells)

  # Getting nearest neighbors
  knnMatrix <- getNearestNeighbors(distanceMatrix, numNearestNeighbors, subSamples)

  # Initialize imputed matrix
  imputedMatrix <- matrix(0, nrow = nrow(rnaCountsMatrix), ncol = numMetaCells)
  rownames(imputedMatrix) <- rownames(rnaCountsMatrix)
  colnames(imputedMatrix) <- subSamples

  # Impute matrix using nearest neighbors
  for (sample in subSamples) {
    neighbors <- knnMatrix[sample, ]
    sampleData <- rnaCountsMatrix[, c(sample, colnames(rnaCountsMatrix)[neighbors])]
    imputedMatrix[, sample] <- rowSums(sampleData)
  }

  # Return the imputed matrix
  imputedMatrix
}


#' Get Nearest Neighbors
#'
#' This function computes the nearest neighbors for a subset of samples
#' in a distance matrix. It identifies the closest cells (neighbors) based
#' on the provided distances, which is useful in various clustering and
#' cell similarity analyses.
#'
#' @param distanceMatrix A distance matrix where each element represents
#' the distance between pairs of samples (e.g., cells).
#' @param numNearestNeighbors Integer, specifying the number of nearest
#' neighbors to identify for each sample.
#' @param subSamples A vector of sample names (or indices) for which
#' the nearest neighbors should be calculated.
#' @return A matrix where each row corresponds to a sample in subSamples
#' and each column contains the indices of the nearest neighbors in the
#' order of proximity.
#'
#' @examples
#' # distanceMatrix is a pre-existing square matrix of distances
#' # subSamples is a vector of sample names or indices
#' nearestNeighbors <- getNearestNeighbors(distanceMatrix, 5, c("Sample1", "Sample2"))
getNearestNeighbors <- function(distanceMatrix, numNearestNeighbors, subSamples) {
  # Validate arguments
  if (!is.matrix(distanceMatrix)) {
    stop("distanceMatrix must be a matrix")
  }

  # Converting distance matrix to matrix format and adjusting diagonal
  distanceMatrix <- as.matrix(distanceMatrix)
  diag(distanceMatrix) <- max(distanceMatrix) + 1

  # Initialize neighbor matrix
  neighborMatrix <- matrix(0, nrow = length(subSamples), ncol = numNearestNeighbors)
  rownames(neighborMatrix) <- subSamples

  # Find nearest neighbors
  for (i in seq_along(subSamples)) {
    neighborMatrix[i, ] <- order(distanceMatrix[subSamples[i], ])[1:numNearestNeighbors]
  }

  # Return the neighbor matrix
  neighborMatrix
}

#' Generate PseudoBulk Seurat Object
#'
#' Joins pseudo bulk matrices and creates a Seurat object with the combined data
#' and corresponding metadata.
#'
#' @param pseudobulkList A list of matrices representing pseudo bulk data,
#'        named by clusters and conditions.
#'
#' @return A Seurat object created from the combined pseudo bulk data.
#' @importFrom Seurat CreateSeuratObject
#' @export
generatePseudoBulkSeurat <- function(pseudobulkList) {
  meta_data_list <- list()
  B_matrix_all <- NULL

  for (cluster in names(pseudobulkList)) {
    B_matrix_cluster <- pseudobulkList[[cluster]]
    for (condition in names(B_matrix_cluster)) {
      B_matrix_condition <- B_matrix_cluster[[condition]]
      meta_data_list[[paste(cluster, condition, sep = "_")]] <- data.frame(
        cell = colnames(B_matrix_condition),
        cluster = cluster,
        condition = condition
      )
      B_matrix_all <- cbind(B_matrix_all, B_matrix_condition)
    }
  }

  meta_data <- do.call("rbind", meta_data_list)
  rownames(meta_data) <- meta_data$cell
  pseudobulkSeurat <- CreateSeuratObject(counts = B_matrix_all, meta.data = meta_data, assay = "RNA")
  return(pseudobulkSeurat)
}



#' Generate Pseudobulk Matrices from Seurat Object
#'
#' Iterates over clusters and conditions in a Seurat object to generate pseudobulk matrices.
#' It handles minimum and maximum cell count constraints for pseudobulk generation.
#'
#' @param seuratObj A Seurat object.
#' @param paramMinMetaCells The minimum number of meta cells to generate.
#' @param paramMaxMetaCells The maximum number of meta cells to generate.
#' @param paramMaxScCells The maximum number of single cells to consider within each condition.
#' @param paramK Parameter 'k' used in calculatePseudoBulkCell function.
#' @param paramPairings pairings for case and control clusters
#' @return A list of pseudobulk matrices organized by cluster and condition.
#' @importFrom dplyr filter rename sample_n
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @export
generateMetaCellMatricesWPairings <- function(seuratObj, paramMinMetaCells = 100, paramMaxMetaCells = 600, paramMaxScCells, paramK,paramPairings) {
  B_matrices <- list()

  for (this_row in c(1:nrow(paramPairings))) {
    cat("Calculating pseudobulk matrices for cluster:", paramPairings$case[this_row], "\n")
    case_cluster <- paramPairings$case[this_row]
    control_cluster <- paramPairings$control[this_row]
    case_cells = seuratObj@meta.data %>% filter(cluster %in% case_cluster & condition == "case") %>% pull(barcode)
    control_cells = seuratObj@meta.data %>% filter(cluster %in% control_cluster & condition == "control") %>% pull(barcode)
    all_cells = c(case_cells,control_cells)

    seuratObjectCluster <- subset(seuratObj, cells = all_cells)
    minCellCount <- min(table(seuratObjectCluster$condition))
    minCellCount <- floor(minCellCount/(paramK+1))

    if (minCellCount < paramMinMetaCells) {
      next
    }

    minCellCount <- min(minCellCount, paramMaxMetaCells)
    cat("Number of pseudobulk cells:", minCellCount, "\n")

    for (this_condition in unique(seuratObjectCluster$condition)) {
      cat("Calculating pseudobulk matrices for condition:", this_condition, "\n")

      conditionData <- seuratObjectCluster@meta.data %>%
        tibble::rownames_to_column(var = "cell") %>%
        filter(condition == this_condition)

      if (nrow(conditionData) > paramMaxScCells) {
        cellsToKeep <- conditionData %>% select(cell) %>% sample_n(size = paramMaxScCells) %>% unlist(use.names = FALSE)
      } else {
        cellsToKeep <- conditionData$cell
      }

      seuratObjectCondition <- subset(seuratObjectCluster, cells = cellsToKeep)

      B_matrix <- calculatePseudoBulkCell(
        seuratObject = seuratObjectCondition,
        numNearestNeighbors = paramK,
        numMetaCells = minCellCount
      )

      B_matrices[[unique(seuratObjectCondition$cluster)]][[this_condition]] <- B_matrix
    }
  }

  return(B_matrices)
}


#' createGroupsFromPairings
#'
#' @param paramPairings a data frame of comparisons (case and control columns) whose values correspond to clusters
#'
#' @return a group from which to calculate minimum sizes
#' @export
#'
#' @examples
createGroupsFromPairings <- function(paramPairings) {
  # Append 'case' and 'control' to cluster names
  df <- transform(paramPairings,
                  case = paste(case, "case"),
                  control = paste(control, "control"))

  # Create a graph from the data frame
  g <- graph_from_data_frame(df)

  # Find connected components
  components <- components(g)

  # Extract each component and remove duplicates
  groups <- lapply(components$membership, function(x) {
    names(which(components$membership == x))
  })

  # Return the unique groups
  return(groups)
}

#' calculateMinimumN
#'
#' @param groups groups from createGroupsFromPairings
#' @param min_counts the minimum number of counts accepted per cluster/condition
#' @param paramPairings a data frame of comparisons (case and control columns) whose values correspond to clusters
#'
#' @return minimum size of each comparison
#' @export
#'
#' @examples
calculateMinimumN <- function(groups, min_counts, paramPairings) {
  paramPairings$min_n <- 0
  # Iterate over groups to calculate and assign minimum 'n'
  for(this_group in groups) {
    # Extract cluster names from group strings
    # Extract the first word when the string contains 'case'
    case_clusters <- sapply(this_group, function(x) {
      if (grepl(" case", x)) {
        return(strsplit(x, " ")[[1]][1])
      }
      return(NA)
    })

    # Filter out NA values if you only want the words extracted from strings containing 'case'
    case_clusters <- case_clusters[!is.na(case_clusters)]

    # Extract the first word when the string contains 'control'
    control_clusters <- sapply(this_group, function(x) {
      if (grepl(" control", x)) {
        return(strsplit(x, " ")[[1]][1])
      }
      return(NA)
    })

    # Filter out NA values if you only want the words extracted from strings containing 'case'
    control_clusters <- control_clusters[!is.na(control_clusters)]

    # Calculate minimum 'n' for the group
    this_min_n <- min_counts %>%
      filter((cluster %in% case_clusters & condition == "case") |
               (cluster %in% control_clusters & condition == "control")) %>%
      summarise(min_n = min(n)) %>%
      pull(min_n)

    # Update the minimum 'n' in paramPairings
    paramPairings$min_n <- ifelse(paramPairings$case %in% case_clusters &
                                    paramPairings$control %in% control_clusters,
                                  this_min_n,
                                  paramPairings$min_n)
  }

  # Return the updated paramPairings
  return(paramPairings)
}


#' Remove Clusters with Insufficient Cells Per Condition
#'
#' This function filters out clusters from a Seurat object where any condition within a cluster
#' has fewer than a specified number of cells, N. It processes the metadata to identify clusters
#' that fail to meet this criterion across their conditions, and then subsets the Seurat object
#' to exclude these clusters.
#'
#' @param seurat_object A Seurat object containing cell metadata with columns `cluster` and `condition`.
#' @param N The minimum number of cells required per condition within each cluster.
#'
#' @return A Seurat object with only clusters that have all conditions meeting the specified cell count threshold.
#'
#' @examples
#' # Assuming 's' is a Seurat object and we require at least 100 cells per condition
#' s_filtered <- KeepClustersWithMtNCellsPerCondition(s, 100)
#'
#' @export
KeepClustersWithMtNCellsPerCondition <- function(seurat_object,N){
  # Convert row names to a column for easier manipulation
  metadata <- seurat_object@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  # Identify clusters that have any condition with fewer than N cells
  insufficient_clusters <- metadata %>%
    group_by(cluster, condition) %>%
    summarize(cell_count = n(), .groups = 'drop') %>%
    group_by(cluster) %>%
    filter(any(cell_count < N)) %>%
    pull(cluster) %>%
    unique()
  # Filter out cells that belong to insufficient clusters
  cells_to_keep <- metadata %>%
    filter(!cluster %in% insufficient_clusters) %>%
    pull(cell)

  seurat_object@meta.data$cell <- NULL

  seurat_object <- subset(seurat_object,cells = cells_to_keep)

  return(seurat_object)
}


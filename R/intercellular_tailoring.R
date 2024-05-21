
#' Retrieve and Format Ligand-Receptor Data from a Forrest Database
#'
#' This function reads a CSV file containing a ligand-receptor database (typically
#' from the Forrest database) and reformats it to a specific structure. It selects
#' ligand and receptor gene symbols from the database and creates a new column combining
#' these two symbols.
#'
#' @param forrest.database.filepath A string specifying the path to the CSV file
#'        containing the Forrest ligand-receptor database.
#'
#' @return A data frame with three columns: 'ligand', 'receptor', and 'lr'.
#'         The 'lr' column is a concatenation of the 'ligand' and 'receptor' columns,
#'         separated by a hyphen.
#'
#' @importFrom utils read.csv
#' @examples
#' # Example usage:
#' filepath <- "path/to/forrest/database.csv"
#' lr_data <- getForrestLRDatabase(filepath)
#'
#' @export
getForrestLRDatabase <- function(forrest.database.filepath) {
  # Validate the input filepath
  if (!file.exists(forrest.database.filepath)) {
    stop("The specified file does not exist.")
  }

  # Read the CSV file and create a simplified data frame
  L.set <- read.csv(forrest.database.filepath, header = TRUE)
  L.set <- L.set[c("Ligand.gene.symbol", "Receptor.gene.symbol")]

  # Rename columns for clarity
  names(L.set) <- c("ligand", "receptor")

  # Create a combined ligand-receptor column
  L.set$lr <- paste(L.set$ligand, L.set$receptor, sep = "-")

  return(L.set)
}



#TODO: getFeatureStatistics here sample size may skew the contribution by each cell-type, this needs to be addressed

#' Calculate Feature Statistics in Seurat Object
#'
#' This function calculates various statistics for specified features across different clusters
#' and conditions in a Seurat object. Statistics include the sum of counts, binary counts,
#' number of cells, fraction of cells with counts, and fractions of total counts for each feature.
#'
#' @param features A vector of feature (gene) names for which to calculate statistics.
#' @param seuratObj A Seurat object containing single-cell RNA sequencing data.
#'
#' @return A data frame with rows representing features and columns representing different statistics:
#'         cluster, condition, feature, sum of counts, binary counts, number of cells, fraction of cells
#'         with counts, total counts for the feature in the condition, and fraction of counts for the
#'         feature in the condition.
#'
#' @importFrom Matrix rowSums
#' @importFrom dplyr bind_rows
#' @importFrom SeuratObject FetchData GetAssayData
#'
#' @export
#'
#' @examples
#' # Assuming 'seurat' is a valid Seurat object and 'features' is a vector of gene names:
#' stats <- getFeatureStatistics(features, seurat)
getFeatureStatistics <- function(features, seuratObj) {
  # Ensure features are unique and present in the Seurat object
  validFeatures <- unique(features[features %in% rownames(seuratObj)])

  # Initialize an empty list for results
  featureStats <- list()

  # Extract metadata and RNA data
  metaData <- SeuratObject::FetchData(seuratObj, vars = c("cluster", "condition"))
  rnaData <- SeuratObject::GetAssayData(seuratObj, assay = "RNA", slot = "data")

  # Iterate through conditions
  for (condition in unique(metaData$condition)) {
    indCondition <- which(metaData$condition == condition)
    totalCountsCondition <- Matrix::rowSums(rnaData[validFeatures, indCondition])

    # Iterate through clusters within each condition
    for (cluster in unique(metaData$cluster[indCondition])) {
      indCluster <- which(metaData$cluster == cluster & metaData$condition == condition)
      sumCounts <- Matrix::rowSums(rnaData[validFeatures, indCluster])
      binaryCounts <- Matrix::rowSums(rnaData[validFeatures, indCluster] != 0)
      nCells <- length(indCluster)
      fracCellsWCounts <- binaryCounts / nCells
      fracCountsFeaturesCondition <- sumCounts / totalCountsCondition[validFeatures]

      stats <-  data.frame(cluster = cluster,
                           condition = condition,
                           feature = validFeatures,
                           sum.counts = sumCounts,
                           binary.counts = binaryCounts,
                           n.cell = nCells,
                           frac.cells.w.counts = fracCellsWCounts,
                           total.counts.feature.condition = totalCountsCondition[validFeatures],
                           frac.counts.features.condition = fracCountsFeaturesCondition)

      featureStats[[paste(cluster, condition, sep = "_")]] <- stats
    }
  }

  # Combine results and handle NA values
  featureStatsDf <- bind_rows(featureStats)
  featureStatsDf$frac.counts.features.condition[is.na(featureStatsDf$frac.counts.features.condition)] <- 0
  featureStatsDf$frac.cells.w.counts[is.na(featureStatsDf$frac.cells.w.counts)] <- 0
  rownames(featureStatsDf) <- NULL
  return(featureStatsDf)
}

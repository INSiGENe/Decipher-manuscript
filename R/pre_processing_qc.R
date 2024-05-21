#' Generate QC Data by Cluster and Condition
#'
#' Aggregates quality control data by cluster and condition from a Seurat object.
#' This function creates a summary data frame that counts the number of cells
#' in each combination of cluster and condition. It is useful for quality control
#' and exploratory data analysis in single-cell studies.
#'
#' @param seuratObj A Seurat object containing single-cell data.
#' @param maxClusterNameLength The maximum length of cluster names to consider.
#'                             Names longer than this length will be truncated.
#'
#' @return A data frame with counts of cells in each combination of cluster and condition.
#'         The data frame has columns for 'cluster', 'condition', and 'n', where 'n'
#'         represents the count of cells.
#'
#' @importFrom dplyr count
#' @importFrom stringr str_sub
#'
#' @export
#'
#' @examples
#' # Assuming 'seurat' is a valid Seurat object with proper metadata:
#' qcData <- generateQCDataByClusterAndCondition(seurat, maxClusterNameLength = 10)
generateQCDataByClusterAndCondition <- function(seuratObj, maxClusterNameLength) {
  # Check if seuratObj is a Seurat object
  if (!inherits(seuratObj, "Seurat")) {
    stop("seuratObj must be a Seurat object")
  }

  # Check if maxClusterNameLength is a positive integer
  if (!is.numeric(maxClusterNameLength) || maxClusterNameLength <= 0) {
    stop("maxClusterNameLength must be a positive integer")
  }

  # Extract base data from Seurat object
  metaData <- seuratObj@meta.data
  baseData <- dplyr::tibble(
    cell = colnames(seuratObj),
    cluster = stringr::str_sub(metaData$cluster, start = 1, end = maxClusterNameLength),
    condition = metaData$condition
  )

  # Aggregate QC plot data
  qcPlotData <- baseData %>%
    dplyr::count(cluster, condition)

  return(qcPlotData)
}


#' Plot Quality Control - Cells per Cluster
#'
#' Creates a scatter plot showing the log2-transformed count of cells per cluster,
#' distinguished by condition. Horizontal lines are drawn at specified cell counts
#' (100, 300, 500). The plot is saved as a PNG file.
#'
#' @param qc_plot_data_cells_per_cluster A data frame with QC data.
#'        Expected columns: 'cluster', 'condition', and 'n' (cell counts).
#' @param outputPath Directory path to save the plot file.
#'        Default is "data/figures".
#'
#' @import ggplot2
#' @importFrom grDevices png
#' @importFrom tidyr spread gather
#' @importFrom dplyr mutate case_when if_else
#'
#' @export
#'
#' @examples
#' # Assuming 'qcData' is a valid data frame with the necessary columns:
#' plotQC_CpC(qcData)
plotQC_CpC <- function(qc_plot_data_cells_per_cluster,outputPath = "data/figures"){
  fileName <- paste("cells_per_cluster",".png",sep="")
  filePath <- file.path(outputPath, fileName)

  df <- qc_plot_data_cells_per_cluster

  # Spread the data so we have 'case' and 'control' as separate columns
  df_spread <- df %>%
    spread(key = condition, value = n)


  # Create a column for color based on the value of n
  df_spread <- df_spread %>%
    mutate(
      color_case = case_when(
        case < 100 ~ "tomato3",
        case >= 100 & case < 300 ~ "tomato1",
        case >= 300 & case < 500 ~ "goldenrod2",
        case >= 500 ~ "forestgreen",
        TRUE ~ "white"  # Assuming you want to have a default color as white
      ),
      color_control = case_when(
        control < 100 ~ "tomato3",
        control >= 100 & control < 300 ~ "tomato1",
        control >= 300 & control < 500 ~ "goldenrod2",
        control >= 500 ~ "forestgreen",
        TRUE ~ "white"
      )
    )

  # Melt the data back for ggplot2
  df_melted <- df_spread %>%
    gather(key = "condition", value = "n", -cluster, -color_case, -color_control)

  # Create a new color column based on condition
  df_melted <- df_melted %>%
    mutate(
      color = dplyr::if_else(condition == "case", color_case, color_control)
    )

  # Plot
  p <- ggplot(df_melted, aes(x = condition, y = reorder(cluster, n), fill = color)) +
    geom_tile(color = "black") +
    scale_fill_identity() +
    geom_text(aes(label = n), color = "black", size = 5) +
    labs(x = "Condition", y = "Cluster", fill = "Cell Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"  # Remove legend if not needed
    )+
    ggtitle("Cells per cluster/condition")

  grDevices::png(filePath,width=12,height=16,units="cm",res=400)
  print(p)
  dev.off()
}


#' Get Clusters Passing Cells per Cluster (CpC) Filter
#'
#' This function filters clusters based on a minimum cells per cluster (CpC) threshold.
#' It returns the clusters that have exactly 2 occurrences after applying the threshold.
#'
#' @param cpcData A data frame with CpC data, expected to contain columns 'n' and 'cluster'.
#' @param minCpc The minimum threshold for CpC to filter the clusters.
#'
#' @return A vector of clusters that pass the CpC filter with exactly 2 occurrences.
#' @importFrom dplyr filter count
#'
#' @export
#'
#' @examples
#' # Assuming 'cpcData' is a valid data frame with the necessary columns:
#' filteredClusters <- getClustersPassingCpCFilter(cpcData, 50)
getClustersPassingCpCFilter <- function(cpcData, minCpc = 100) {
  cpcData %>%
    dplyr::filter(n > minCpc) %>%
    dplyr::count(cluster) %>%
    dplyr::filter(n == 2) %>%
    dplyr::pull(cluster)
}

#' Plot Quality Control - UMIs per Cluster
#'
#' Creates a boxplot showing the distribution of Unique Molecular Identifiers (UMIs) per cluster,
#' distinguished by condition. Horizontal lines are drawn at specified UMI counts for reference.
#'
#' @param seuratObject A Seurat object containing single-cell RNA sequencing data.
#' @param outputPath Directory path to save the plot file.
#'        Default is "data/figures".
#'
#' @details
#' This function calculates the sum of RNA counts (UMIs) for each cell in the Seurat object and then
#' creates a boxplot for these UMIs per cluster. Horizontal lines are added at 5000, 7000, and 10000
#' UMIs for visual reference.
#'
#' The plot is saved as a PNG file named "UMIs_per_cluster.png" in the specified output directory.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline theme element_text ggtitle ylim
#' @importFrom grDevices png
#' @importFrom Seurat GetAssayData
#' @importFrom stringr str_sub
#' @importFrom grDevices dev.off
#' @importFrom Matrix colSums
#'
#' @export
#'
#' @examples
#' # Assuming 'seurat' is a valid Seurat object:
#' plotQC_UpC(seurat)
plotQC_UpC <- function(seuratObject,outputPath = "data/figures"){
  rnaCountsMatrix <- Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "counts")
  UpC <- colSums(rnaCountsMatrix)
  UpCdf <- data.frame(
    cell = colnames(seuratObject),
    cluster = stringr::str_sub(seuratObject@meta.data$cluster, start = 1, end = 20),
    UpC = UpC,
    condition = seuratObject@meta.data$condition
  )

  p <- ggplot(UpCdf,aes(y = UpC,x=cluster,colour=condition))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_hline(aes(yintercept = 5000),colour = "firebrick1")+
    geom_hline(aes(yintercept = 7000), colour = "darkgoldenrod1")+
    geom_hline(aes(yintercept = 10000), colour = "chartreuse4")+
    ggtitle(label = "UMIs per Cluster")+
    ylim(0,100000)

  fileName <- paste("UMIs_per_cluster",".png",sep="")
  filePath <- file.path(outputPath, fileName)
  png(filePath,width=30,height=20,units="cm",res=400)
  print(p)
  dev.off()
}

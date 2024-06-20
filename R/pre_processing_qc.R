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
#' Creates a tile plot showing the counts of cells per cluster,
#' separated by condition. Colours are drawn at specified cell counts
#' (100, 300, 500), red indicates low cell counts, yellow moderate, green high. The plot is saved as a PNG file.
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

  #grDevices::png(filePath,width=12,height=16,units="cm",res=400)
  #print(p)
  #dev.off()
  return(p)
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

#' Filter Ligands Based on Expression Threshold in Clusters
#'
#' This function extracts ligands that are expressed above a specified threshold in one or more
#' clusters from a Seurat object. It integrates information about ligands and receptors from
#' a given dataset, computes feature statistics for these molecules, and filters out those ligands
#' that meet the expression criteria.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data.
#' @param L.set A dataframe containing columns `ligand` and `receptor`, which lists
#'        the ligands and their corresponding receptors.
#' @param param_min_ligand_expr_in_cluster The minimum fraction of cells (expressed as a decimal)
#'        within any cluster that must express the ligand for it to be considered expressed.
#'
#' @return A vector of ligands that are expressed above the specified threshold across clusters.
#'
#' @examples
#' # Assuming 'decipher_seurat' is a Seurat object and 'L.set' is a dataframe
#' # with ligand and receptor information:
#' expressed_ligands <- getFilteredLigands(
#'   decipher_seurat = seurat_object,
#'   L.set = ligand_receptor_set,
#'   param_min_ligand_expr_in_cluster = 0.1
#' )
#'
#' @importFrom dplyr filter pull
#' @export
getFilteredLigands <- function(decipher_seurat,L.set,param_min_ligand_expr_in_cluster){

  receptors <- unique(L.set$receptor)
  ligands <- unique(L.set$ligand)
  all_ligand_receptors <- unique(c(ligands,receptors))

  feature_statistics <- getFeatureStatistics(
    features=all_ligand_receptors,
    seuratObj=decipher_seurat)

  expressed_ligands <- feature_statistics %>%
    dplyr::filter(feature %in% ligands & frac.cells.w.counts > param_min_ligand_expr_in_cluster) %>%
    dplyr::pull(feature) %>%
    unique()

  return(expressed_ligands)
}


#' Get Expressed Receptors for Each Cluster
#'
#' This function identifies and retrieves expressed receptors for each cluster in a Seurat object.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param L.set A set of ligands used for filtering receptors.
#'
#' @return A list where each element corresponds to a cluster and contains the expressed receptors for that cluster.
#'
#' @details The function iterates through each unique cluster in the Seurat object and applies the `getFilteredReceptorsForCluster` function to identify receptors expressed in that cluster. The parameter `param_min_receptor_expr_in_cluster` is set to 0.1 to filter receptors based on their expression level.
#'
#' @examples
#' \dontrun{
#' seurat_object <- CreateSeuratObject(counts = your_counts_matrix)
#' L.set <- c("Ligand1", "Ligand2", "Ligand3")
#' expressed_receptors <- getExpressedReceptorsForEachCluster(seurat_object, L.set)
#' }
#'
#' @export
getExpressedReceptorsForEachCluster <- function(seurat_object, L.set) {

  expressed_receptors_all_clusters <- list()

  for(this_cluster in unique(seurat_object$cluster)){

    expressed_receptors_this_cluster <- getFilteredReceptorsForCluster(
      seurat_object,
      L.set,
      param_min_receptor_expr_in_cluster = 0.1,
      this_cluster)

    expressed_receptors_all_clusters[[this_cluster]] <- expressed_receptors_this_cluster

  }
  return(expressed_receptors_all_clusters)

}



#' Filter Receptors Based on Expression Threshold within a Specific Cluster
#'
#' This function filters out receptors based on their expression levels within a specific cluster
#' of a Seurat object. It uses a dataset that contains information about both ligands and receptors,
#' calculates feature statistics for these molecules across the dataset, and returns receptors
#' that are expressed above a specified threshold in a specified cluster.
#'
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data.
#' @param L.set A dataframe containing columns `ligand` and `receptor`, which lists
#'        the ligands and their corresponding receptors.
#' @param param_min_receptor_expr_in_cluster The minimum fraction of cells (expressed as a decimal)
#'        within the specified cluster that must express the receptor for it to be considered expressed.
#' @param this_cluster The cluster within the Seurat object for which receptors are being filtered.
#'
#' @return A vector of receptors that are expressed above the specified threshold within the given cluster.
#'
#' @examples
#' # Assuming 'decipher_seurat' is a Seurat object and 'L.set' is a dataframe
#' # with ligand and receptor information:
#' filtered_receptors <- getFilteredReceptorsForCluster(
#'   decipher_seurat = seurat_object,
#'   L.set = ligand_receptor_set,
#'   param_min_receptor_expr_in_cluster = 0.1,
#'   this_cluster = "Cluster_1"
#' )
#'
#' @importFrom dplyr filter pull
#' @export
getFilteredReceptorsForCluster <- function(decipher_seurat,L.set,param_min_receptor_expr_in_cluster,this_cluster){

  receptors <- unique(L.set$receptor)
  ligands <- unique(L.set$ligand)
  all_ligand_receptors <- unique(c(ligands,receptors))

  feature_statistics <- getFeatureStatistics(
    features=all_ligand_receptors,
    seuratObj=decipher_seurat)

  expressed_receptor <- feature_statistics %>%
    dplyr::filter(cluster == this_cluster) %>%
    dplyr::filter(frac.cells.w.counts > param_min_receptor_expr_in_cluster) %>%
    dplyr::pull(feature) %>%
    unique()

  return(expressed_receptor)

}


#' Get Relevant Features for Each Cluster
#'
#' This function identifies and retrieves relevant ligand-receptor features for each cluster based on expressed ligands and receptors.
#'
#' @param L.set A data frame containing ligand-receptor pairs.
#' @param expressed_ligands A vector of expressed ligands.
#' @param expressed_receptors_all_clusters A list of expressed receptors for each cluster, typically obtained from `getExpressedReceptorsForEachCluster`.
#'
#' @return A list where each element corresponds to a cluster and contains the relevant ligand-receptor features for that cluster.
#'
#' @details The function iterates through each cluster in the `expressed_receptors_all_clusters` list, filters the `L.set` data frame to include only those ligand-receptor pairs where the receptor is expressed in the cluster and the ligand is in the `expressed_ligands` vector. The filtered results are then stored in a list, with each element corresponding to a cluster.
#'
#' @examples
#' \dontrun{
#' L.set <- data.frame(ligand = c("Ligand1", "Ligand2"), receptor = c("Receptor1", "Receptor2"))
#' expressed_ligands <- c("Ligand1", "Ligand3")
#' expressed_receptors_all_clusters <- list(Cluster1 = c("Receptor1"), Cluster2 = c("Receptor2"))
#' relevant_features <- getRelevantFeaturesForEachCluster(L.set, expressed_ligands, expressed_receptors_all_clusters)
#' }
#'
#' @export
getRelevantFeaturesForEachCluster <- function(L.set, expressed_ligands, expressed_receptors_all_clusters){

  L_set_relevant_features_all_clusters <- list()

  for(this_cluster in names(expressed_receptors_all_clusters)){
    expressed_receptors_this_cluster <- expressed_receptors_all_clusters[[this_cluster]]

    L_set_relevant_features <- L.set %>%
      filter(receptor %in% expressed_receptors_this_cluster & ligand %in% expressed_ligands)

    L_set_relevant_features_all_clusters[[this_cluster]] <- L_set_relevant_features
  }

  return(L_set_relevant_features_all_clusters)
}


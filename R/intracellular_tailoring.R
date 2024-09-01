
#' Trim Gene Regulatory Network (GRN) Data Frame
#'
#' This function trims a GRN data frame based on specified criteria.
#' It filters the GRN by p-value, limits the number of edges based on their absolute coefficients,
#' and removes sources with fewer targets than a specified minimum.
#'
#' @param grn_df A data frame representing the gene regulatory network, expected to have
#' columns for p-values and source-target relationships.
#' @param pValue A numeric threshold for p-values; edges with p-values above this threshold will be discarded.
#' @param topEdges The maximum number of edges to retain, based on the highest absolute coefficients.
#' If NULL, no limit is applied.
#' @param minTargets The minimum number of targets required per source; sources with fewer targets than this
#' will be removed from the network.
#' @return A trimmed GRN data frame.
#' @importFrom dplyr filter slice_max arrange count pull
#' @export
#'
#' @examples
#' # grn_df is a data frame representing a gene regulatory network
#' trimmedGRN <- trimGRN(grn_df, 0.05, 100, 3)
#'
#' @note The function assumes that 'coef_abs' is a column in `grn_df` used for sorting edges
#' in `applyTopEdgesFilter`. If your data frame has a different column for edge coefficients,
#' you'll need to modify the function accordingly.
trimGRN <- function(grn_df, pValue, topEdges, minTargets) {
  # Validate input
  if (!is.data.frame(grn_df)) {
    stop("grn_df must be a data frame")
  }

  grn_df %>%
    filter(p < pValue) %>%
    applyTopEdgesFilter(topEdges) %>%
    removeSourcesWithFewTargets(minTargets)
}
applyTopEdgesFilter <- function(df, topEdges) {
  if (!is.null(topEdges)) {
    df %>%
      arrange(desc(coef_abs)) %>%
      dplyr::slice_max(n = min(nrow(df),topEdges),order_by = coef_abs)
  } else {
    df
  }
}
removeSourcesWithFewTargets <- function(df, minTargets) {
  sourcesToRemove <- df %>%
    count(source) %>%
    filter(n < minTargets) %>%
    pull(source)

  if (length(sourcesToRemove) > 0) {
    df %>%
      filter(!source %in% sourcesToRemove)
  } else {
    df
  }
}


#TODO: this assumes that the dataframe is arranged in the right order (descending coef_abs I assume)
#make it such that this doesn't have to be an assumption

#' Cap Regulon by Top N Targets
#'
#' This function limits the number of targets for each source in a gene regulatory network (GRN)
#' to the top N entries.
#'
#' @param this_grn A data frame representing the gene regulatory network, expected to contain
#'        a column 'source'.
#' @param n_top The maximum number of top entries to retain for each source.
#'        Default is 100.
#'
#' @return A data frame similar to `this_grn` but with each source limited to a maximum of `n_top` entries.
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' # Assuming 'grn' is a gene regulatory network data frame:
#' cappedRegulon <- capRegulon(grn, n_top = 100)
capRegulon <- function(this_grn, n_top = 100) {
  sources <- unique(this_grn$source)

  cappedEntries <- lapply(sources, function(source) {
    this_grn_source <- this_grn[this_grn$source == source,]
    head(this_grn_source, n_top)
  })

  result <- dplyr::bind_rows(cappedEntries)
  return(result)
}


#' Cap Regulon by Top N Targets
#'
#' This function limits the number of targets for each source in a gene regulatory network (GRN)
#' to the top N entries.
#'
#' @param this_grn A data frame representing the gene regulatory network, expected to contain
#'        a column 'source'.
#' @param n_top The maximum number of top entries to retain for each source.
#'        Default is 100.
#'
#' @return A data frame similar to `this_grn` but with each source limited to a maximum of `n_top` entries.
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' # Assuming 'grn' is a gene regulatory network data frame:
#' cappedRegulon <- capRegulon(grn, n_top = 100)
capRegulon_2 <- function(this_grn, n_top = 100) {
  sources <- unique(this_grn$source)

  cappedEntries <- lapply(sources, function(source) {
    this_grn_source <- this_grn[this_grn$source == source,]
    this_grn_source %>%
      arrange(-coef_abs) %>%
      top_n(n_top,wt=coef_abs)
    #head(this_grn_source, n_top)
  })

  result <- dplyr::bind_rows(cappedEntries)
  return(result)
}

#' Retrieve and Optionally Trim Gene Regulatory Network (GRN) for a Specific Cluster
#'
#' This function loads a gene regulatory network (GRN) from a CSV file specific to a given cluster.
#' It optionally trims the GRN based on thresholds specified for p-values, the number of top edges,
#' and the minimum number of targets, using the CellOracle toolkit if enabled.
#'
#' @param output_filepath The base directory where the CellOracle data and GRN files are stored.
#' @param this_cluster The name of the cluster for which the GRN is to be retrieved.
#' @param flag.co.grn Logical flag indicating whether to trim the GRN using CellOracle criteria.
#'        If TRUE, the GRN is trimmed according to specified parameters.
#'
#' @return A data frame representing the GRN for the specified cluster. If no GRN file is found,
#'         a warning is issued and the function returns `NULL`.
#'
#' @examples
#' # Assuming the necessary files are in 'path/to/output', and you're interested in 'Cluster_1':
#' regulon_data <- getRegulon(
#'   output_filepath = "path/to/output",
#'   this_cluster = "Cluster_1",
#'   flag.co.grn = TRUE
#' )
#'
#' @importFrom utils read.csv
#' @export
getRegulon <- function(output_filepath,this_cluster,flag.co.grn=TRUE){

regulon_filepath <- file.path(output_filepath,"cellOracle/data/GRN",paste(this_cluster,"csv",sep="."))
# check if regulon exists
if(!file.exists(regulon_filepath)){
  warning("no cell oracle grn found for this cluster")
  next
}

regulon_this_cluster <- read.csv(regulon_filepath)
regulon_this_cluster <- regulon_this_cluster[,-1]

if(flag.co.grn){
  #trim GRN using CellOracle results
  regulon_this_cluster <- trimGRN(
    grn_df = regulon_this_cluster,
    pValue = 0.01,
    topEdges = 20000,
    minTargets = 20)
}

return(regulon_this_cluster)
}



#' Add Random Gene Regulatory Networks to a Cluster's GRN
#'
#' This function augments a gene regulatory network (GRN) for a specific cluster with random GRNs.
#' The random GRNs are generated based on a reference GRN to simulate additional regulatory
#' relationships. These random GRNs are combined with the original GRN to expand the dataset,
#' which might be useful for simulation or robustness testing.
#'
#' @param decipher_seurat_this_cluster A Seurat object specific to a cluster, from which gene
#'        names are extracted to generate random GRNs.
#' @param regulon_this_cluster_capped A data frame representing the capped GRN of the cluster.
#'        This GRN is used as the reference to generate random GRNs.
#'
#' @return A data frame that combines the original capped GRN with additional random GRNs.
#'
#' @examples
#' # Assuming 'decipher_seurat_cluster' is a Seurat object for a specific cluster and
#' # 'regulon_cluster_capped' is the existing GRN for that cluster:
#' enhanced_grn <- addRandomGRNs(
#'   decipher_seurat_this_cluster = decipher_seurat_cluster,
#'   regulon_this_cluster_capped = regulon_cluster_capped
#' )
#'
#' @importFrom stats rbind
#' @export
addRandomGRNs <- function(decipher_seurat_this_cluster,regulon_this_cluster_capped){

  random_grns <- generateRandomGRNsFromReference(
    all_genes = rownames(decipher_seurat_this_cluster),
    reference_grns = regulon_this_cluster_capped
  )

  regulon_this_cluster_capped <- rbind(regulon_this_cluster_capped,random_grns)

  return(regulon_this_cluster_capped)
}

#' Get Regulons for All Clusters
#'
#' This function retrieves regulons for each cluster in a Seurat object from a specified file.
#'
#' @param filepath A string representing the path to the file containing regulon data.
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#'
#' @return A list where each element corresponds to a cluster and contains the regulons for that cluster.
#'
#' @details The function iterates through each unique cluster in the Seurat object and applies the `getRegulon` function to retrieve regulons for that cluster from the specified file. The results are stored in a list, with each element corresponding to a cluster.
#'
#' @examples
#' \dontrun{
#' filepath <- "path/to/regulon/file"
#' seurat_object <- CreateSeuratObject(counts = your_counts_matrix)
#' regulons <- getRegulonsAllClusters(filepath, seurat_object)
#' }
#'
#' @export
getRegulonsAllClusters <- function(filepath, seurat_object) {
  regulons_all_clusters <- list()
  for(this_cluster in unique(seurat_object$cluster)){

    regulons_all_clusters[[this_cluster]] <- getRegulon(filepath, this_cluster)

  }
  return(regulons_all_clusters)
}

#' Cap Regulons for All Clusters
#'
#' This function processes and caps the regulons for each cluster, normalizing the Seurat object if specified.
#'
#' @param regulons_all_clusters A list of regulons for each cluster, typically obtained from `getRegulonsAllClusters`.
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param flag.normalize.non.log A logical flag indicating whether to normalize non-log-transformed data.
#'
#' @return A list where each element corresponds to a cluster and contains the capped regulons for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if `flag.normalize.non.log` is TRUE, caps the regulons to the top 40 features, and adds random gene regulatory networks (GRNs) to the capped regulons.
#'
#' @examples
#' \dontrun{
#' regulons_all_clusters <- getRegulonsAllClusters(filepath, seurat_object)
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' capped_regulons <- capRegulonsAllClusters(regulons_all_clusters, decipher_seurat, TRUE)
#' }
#'
#' @export
capRegulonsAllClusters <- function(regulons_all_clusters, decipher_seurat, flag.normalize.non.log) {
  capped_regulons_all_clusters <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){
    regulon_this_cluster <- regulons_all_clusters[[this_cluster]]
    # main object
    decipher_seurat_this_cluster <- decipher_seurat[, which(decipher_seurat$cluster == this_cluster), seed=NULL]

    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    regulon_this_cluster_capped <- capRegulon(regulon_this_cluster, n_top = 40)
    regulon_this_cluster_capped <- addRandomGRNs(decipher_seurat_this_cluster, regulon_this_cluster_capped)

    capped_regulons_all_clusters[[this_cluster]] <- regulon_this_cluster_capped
  }

  return(capped_regulons_all_clusters)
}





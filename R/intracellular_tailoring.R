
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





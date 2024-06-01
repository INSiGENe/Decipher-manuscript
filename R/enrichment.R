
#' Get Differentially Expressed Targets for Regulons
#'
#' This function identifies differentially expressed targets for a given set of regulons
#' in a Seurat object.
#'
#' @param seuratObj A Seurat object containing gene expression data.
#' @param regulonNames A vector of regulon names to analyze.
#' @param logFcThreshold Log fold change threshold for identifying markers.
#' @param grnDf Gene regulatory network data frame.
#' @param targetCt Target count.
#' @return A list containing markers for each regulon.
getDifferentiallyExpressedTargetsForRegulons <- function(seuratObj, regulonNames, logFcThreshold, grnDf, targetCt) {
  Idents(seuratObj) <- seuratObj$condition

  regulonWithDiffTargets <- list()

  for (regulon in regulonNames) {
    markers <- findMarkersForRegulon(seuratObj, regulon, grnDf, logFcThreshold)
    if (nrow(markers) > 0) {
      markers <- annotateMarkers(markers, targetCt, regulon)
      regulonWithDiffTargets[[regulon]] <- markers
    }
  }

  return(regulonWithDiffTargets)
}

findMarkersForRegulon <- function(seuratObj, regulon, grnDf, logFcThreshold) {
  targetGenes <- getTargetGenesForRegulon(regulon, grnDf, seuratObj)
  FindMarkers(
    seuratObj,
    ident.1 = "case",
    ident.2 = "control",
    features = targetGenes,
    logfc.threshold = logFcThreshold,
    only.pos = FALSE
  )
}

getTargetGenesForRegulon <- function(regulon, grnDf, seuratObj) {
  targets <- grnDf$target[grnDf$source == regulon]
  intersect(rownames(seuratObj), targets)
}

annotateMarkers <- function(markers, targetCt, regulon) {
  markers$ct <- targetCt
  markers$regulon <- regulon
  markers$tg_gene <- rownames(markers)
  markers
}


#' Calculate Enrichment Statistics
#'
#' This function calculates enrichment statistics, including a p-value, for a specified gene set
#' against another gene set. It uses the hypergeometric test to compute the p-value, which indicates
#' the probability of observing the given overlap between the two gene sets.
#'
#' @param other_gene_set A vector of gene symbols representing the comparison gene set.
#' @param my_gene_set A vector of gene symbols representing the gene set of interest.
#' @param base_num A numeric value representing the total number of possible genes.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{p_value}{The p-value calculated using the hypergeometric test.}
#'   \item{genes}{A semicolon-separated string of genes found in both gene sets.}
#'   \item{gene_in_my_set_in_other_set}{The number of genes in the gene set of interest that are also in the comparison set.}
#'   \item{gene_in_other_set}{The total number of genes in the comparison set.}
#'   \item{gene_in_my_set}{The total number of genes in the gene set of interest.}
#' }
#'
#' @examples
#' # Example usage:
#' my_genes <- c("gene1", "gene2", "gene3")
#' other_genes <- c("gene2", "gene4", "gene5")
#' result <- calculateEnrichmentStatistics(other_genes, my_genes, 1000)
#'
#' @export
calculateEnrichmentStatistics <- function(other_gene_set, my_gene_set, base_num) {
  # Validate input
  if (!is.numeric(base_num) || base_num <= 0) {
    stop("base_num must be a positive number.")
  }

  # Calculate intersection and lengths
  intersect_genes <- intersect(my_gene_set, other_gene_set)
  num_intersect = length(intersect_genes)
  num_my_set = length(my_gene_set)
  num_other_set = length(other_gene_set)
  num_not_in_my_set = base_num - num_my_set
  num_not_in_my_set_in_other_set = num_other_set - num_intersect

  # Calculate the p-value using hypergeometric test
  p_value <- sum(dhyper(num_intersect:num_other_set,
                        num_other_set,
                        num_not_in_my_set,
                        num_my_set))

  return(list(
    p_value = p_value,
    genes = paste(intersect_genes, collapse = ";"),
    gene_in_my_set_in_other_set = num_intersect,
    gene_in_other_set = num_other_set,
    gene_in_my_set = num_my_set
  ))
}


#' Enrichment Analysis of Differential Expression and Regulon Results
#'
#' This function performs enrichment analysis using predefined gene sets from differential
#' expression markers and regulon analysis. It separates positive and negative markers, identifies
#' top markers based on log fold changes, and performs enrichment using a given Enrichr database.
#' The function also handles multiple enrichment databases and adjusts p-values for multiple testing.
#'
#' @param de_markers_this_cluster A dataframe containing differential expression markers with
#'        columns for average log2 fold change ('avg_log2FC') and gene names ('gene').
#' @param significant_regulon_deltas_this_cluster A dataframe containing significant regulons with
#'        at least a column for regulon names ('name').
#' @param regulon_this_cluster A dataframe listing target genes for each regulon source.
#' @param enrichr_database A list of databases, each an enrichment dataset to be used for analysis.
#'
#' @return A list where each element corresponds to a regulon and contains a dataframe of combined
#'         enrichment results across databases, including adjusted p-values and associated terms.
#'
#' @examples
#' # Assuming appropriate data structures are already created:
#' results <- enrichResults(
#'   de_markers_this_cluster = differential_markers,
#'   significant_regulon_deltas_this_cluster = significant_regulons,
#'   regulon_this_cluster = regulon_targets,
#'   enrichr_database = loaded_enrichr_dbs
#' )
#'
#' @importFrom dplyr filter slice_max slice_min select
#' @importFrom stats p.adjust
#' @export
enrichResults <- function(de_markers_this_cluster,significant_regulon_deltas_this_cluster,regulon_this_cluster,enrichr_database){
  all_pos <- de_markers_this_cluster %>%
    filter(avg_log2FC > 0) %>%
    slice_max(avg_log2FC,n=300)%>%
    select(gene)

  all_neg <- de_markers_this_cluster %>%
    filter(avg_log2FC < 0) %>%
    slice_min(avg_log2FC,n=300)%>%
    select(gene)

  all_pos <- all_pos$gene
  all_neg <- all_neg$gene
  all_pos_neg <- c(all_pos,all_neg)

  my_gene_sets <- list()
  for(this_regulon in significant_regulon_deltas_this_cluster$name){
    all_genes <- regulon_this_cluster$target[regulon_this_cluster$source == this_regulon]
    my_gene_sets[[this_regulon]] <- all_genes
  }

  my_gene_sets[["all_pos"]] <- all_pos
  my_gene_sets[["all_neg"]] <- all_neg
  my_gene_sets[["all_pos_neg"]] <- all_pos_neg

  gene_set_results <- list()
  for(this_gene_set in names(my_gene_sets)){
    dbs_results <- list()
    for(this_dbs_name in names(enrichr_database)){
      term_results <- list()
      this_dbs <- enrichr_database[[this_dbs_name]]
      for(this_term in names(this_dbs)){
        term_results[[this_term]] <- calculateEnrichmentStatistics(this_dbs[[this_term]],my_gene_sets[[this_gene_set]],20000)
      }
      dbs_results[[this_dbs_name]] <- term_results
    }
    gene_set_results[[this_gene_set]] <- dbs_results
  }

  regulon_results_df <- list()
  dbs_results_df <- list()
  for(this_regulon in names(gene_set_results)){
    for(this_dbs_name in names(enrichr_database)){
      dbs_results_df[[this_dbs_name]]  <- do.call(rbind.data.frame, gene_set_results[[this_regulon]][[this_dbs_name]])
      dbs_results_df[[this_dbs_name]]$database  <- rep(this_dbs_name,nrow(dbs_results_df[[this_dbs_name]]))
      dbs_results_df[[this_dbs_name]]$p_value_adjusted  <- p.adjust(dbs_results_df[[this_dbs_name]]$p_value, method = "BH", n = length(dbs_results_df[[this_dbs_name]]$p_value))
      dbs_results_df[[this_dbs_name]]$Term <- rownames(dbs_results_df[[this_dbs_name]])
    }
    regulon_results_df[[this_regulon]] <- do.call(rbind.data.frame, dbs_results_df)
  }
  return(regulon_results_df)
}




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


#' Generate Random GRNs From Reference
#'
#' This function generates random Gene Regulatory Networks (GRNs)
#' based on a provided reference set of GRNs. For each unique source
#' in the reference, it samples a specified number of genes and
#' constructs a data frame with default statistical values.
#'
#' @param all_genes A vector of all possible genes to be sampled.
#' @param reference_grns A data frame representing the reference GRNs,
#'        containing at least the 'source' column.
#'
#' @return A data frame representing the randomly generated GRNs,
#'         with each row corresponding to a gene from the sampled genes,
#'         and columns for source, target, coef_mean, coef_abs, p, and X.logp.
#'
#' @examples
#' # Example usage:
#' all_genes_example <- c("gene1", "gene2", "gene3", "gene4")
#' reference_grns_example <- data.frame(
#'   source = c("A", "B", "A", "C"),
#'   other_columns = 1:4
#' )
#' random_grns <- generateRandomGRNsFromReference(all_genes_example, reference_grns_example)
#'
#' @export
generateRandomGRNsFromReference <- function(all_genes, reference_grns) {
  # Extract unique sources from the reference GRNs
  unique_sources <- unique(reference_grns$source)

  # Initialize a list to hold data frames for each unique source
  list_sample_grn_df <- vector("list", length(unique_sources))

  # Loop through each unique source
  for (i in seq_along(unique_sources)) {
    # Determine the number of samples for the current source
    num_samples <- table(reference_grns$source)[unique_sources[i]]

    # Randomly sample genes based on the number of samples
    sample_grn <- sample(all_genes, num_samples)

    # Create a data frame for the current source with default values
    list_sample_grn_df[[i]] <- data.frame(
      source = paste("sample", i, sep = ""), # Source name, e.g., "sample1", "sample2", etc.
      target = sample_grn,                  # Randomly sampled genes
      coef_mean = 0,                        # Default coefficient mean
      coef_abs = 0,                         # Default absolute coefficient
      p = 0,                                # Default p-value
      X.logp = 0                            # Default log p-value
    )
  }

  # Combine all data frames in the list into one data frame
  all_sample_grn_df <- do.call(rbind, list_sample_grn_df)

  # Return the combined data frame
  return(all_sample_grn_df)
}


#' Calculate Differential Regulation for Regulons
#'
#' This function computes the difference in median values for each regulon
#' (row) in a matrix between two classes. The classes are defined in the
#' provided class vector, distinguishing between 'case' and other categories.
#'
#' @param regulonMatrix A numeric matrix where rows represent regulons (or pathways)
#' and columns represent different samples.
#' @param classVector A character vector indicating the class for each sample.
#' It should have the same length as the number of columns in `regulonMatrix`.
#' The function calculates the difference between the medians of the 'case' class
#' and other classes.
#' @return A data frame where each row corresponds to a regulon and contains
#' the differential regulation value (`deltaPagoda`) and the name of the regulon.
#' @export
#'
#' @examples
#' # regulonMatrix is a matrix of regulon activities
#' # classVector is a vector of class labels for each sample
#' deltas <- getRegulonDeltas(regulonMatrix, classVector)
getRegulonDeltas <- function(regulonMatrix, classVector) {
  # Validate arguments
  if (!is.matrix(regulonMatrix)) {
    stop("regulonMatrix must be a matrix")
  }
  if (!is.vector(classVector) || length(classVector) != ncol(regulonMatrix)) {
    stop("classVector must be a vector of the same length as the number of columns in regulonMatrix")
  }

  # Identify indices for 'case' and 'control'
  caseIndices <- which(classVector == "case")
  controlIndices <- setdiff(seq_len(ncol(regulonMatrix)), caseIndices)

  # Calculate delta values using vectorized operations
  deltaValues <- matrixStats::rowMedians(regulonMatrix[, caseIndices]) - matrixStats::rowMedians(regulonMatrix[, controlIndices])

  # Create output data frame
  deltaDataFrame <- data.frame(
    deltaPagoda = deltaValues,
    name = rownames(regulonMatrix)
  )
  rownames(deltaDataFrame) <- deltaDataFrame$name

  deltaDataFrame <- deltaDataFrame %>%
    mutate(class = ifelse(stringr::str_detect(name,"sample"),"random","real"))

  return(deltaDataFrame)
}

#' Calculate Regulon Scores Using Pagoda2
#'
#' This function calculates regulon scores for a given Seurat object and gene regulatory network (GRN) data frame.
#' It involves creating a Pagoda2 object, preparing the GRN environment, running the Pagoda2 analysis, and processing the results.
#'
#' @param seuratObject A Seurat object containing single-cell RNA sequencing data.
#' @param grn_df A data frame representing the gene regulatory network with columns specifying gene interactions.
#' @return A matrix (RMatrixSC) containing the regulon scores calculated from the Pagoda2 analysis.
#'
#' @import pagoda2
#'
#' @examples
#' # seuratObject is a pre-existing Seurat object
#' # grn_df is a data frame representing the gene regulatory network
#' regulonScores <- getRegulonScores(seuratObject, grn_df)
#'
#' @note The function assumes that the 'RNA' assay and the 'data' slot in the Seurat object are appropriate for the analysis.
#' It also utilizes several helper functions: 'createAndPreparePagoda2', 'prepGRN', 'prepareGRNEnvironment', and 'processPagodaResults'.
#' These functions need to be defined and available in the environment.
getRegulonScores <- function(seuratObject, grn_df) {
  # Validate arguments
  if (!inherits(seuratObject, "Seurat")) {
    stop("seuratObject must be a Seurat object")
  }
  if (!is.data.frame(grn_df)) {
    stop("grn_df must be a data frame")
  }

  # Access counts data
  # TODO: check if data slot is appropriate slot
  counts <- Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "data")

  # Pagoda2 object creation and initial steps
  pagodaObject <- createAndPreparePagoda2(counts)

  # Prepare GRN and create environment
  grnSplit <- prepGRN(grn_df = grn_df)
  goEnv <- prepareGRNEnvironment(grnSplit)

  # Run Pagoda2 analysis
  tryCatch({
    pagodaObject$testPathwayOverdispersion(
      goEnv, verbose = FALSE, correlation.distance.threshold = 0.2, recalculate.pca = FALSE, top.aspects = 15
    )
  }, error = function(e) {
    message("Error in Pagoda2 analysis: ", e$message)
  })

  # Process results
  RMatrixSC <- processPagodaResults(pagodaObject,grnSplit)

  return(RMatrixSC)
}



#' Create and prepare pagoda2
#'
#' @param counts
#'
#' @return
#' @export
#'
#' @import pagoda2
#' @importFrom igraph infomap.community multilevel.community
#'
#' @examples
createAndPreparePagoda2 <- function(counts) {
  pagodaObject <- Pagoda2$new(counts, log.scale = TRUE, n.cores = 1)
  pagodaObject$adjustVariance(plot = TRUE, gam.k = 10)
  pagodaObject$calculatePcaReduction(nPcs = 50, n.odgenes = 3e3)
  pagodaObject$makeKnnGraph(k = 40, type = 'PCA', center = TRUE, distance = 'cosine')
  pagodaObject$getKnnClusters(method = infomap.community, type = 'PCA')
  return(pagodaObject)
}
prepGRN <- function(grn_df){
  grn_df_to_split <- grn_df[,c("source","target")]
  grn_split <- split(grn_df_to_split,grn_df_to_split$source)
  return(grn_split)
}
prepareGRNEnvironment <- function(grnSplit) {
  goEnv <- new.env(parent = globalenv())
  for (i in seq_along(grnSplit)) {
    genes <- as.character(grnSplit[[i]]$target)
    name <- as.character(names(grnSplit)[i])
    assign(name, genes, envir = goEnv)
  }
  return(goEnv)
}
processPagodaResults <- function(pagodaObject,grnSplit) {
  # Ensure the object is of the correct type
  if (!inherits(pagodaObject, "Pagoda2")) {
    stop("pagodaObject must be a Pagoda2 object")
  }

  # Get all pathways
  allPathways <- pagodaObject$misc$pathwayODInfo$name

  # Process each pathway
  for (pathway in allPathways) {
    pathwayData <- pagodaObject$misc$pwpca[[pathway]]
    if (is.null(pathwayData$xp$scores)) {
      updatedScores <- getUpdatedPathwayScores(pagodaObject, pathwayData, pathway,grnSplit)
      pagodaObject$misc$pwpca[[pathway]]$xp$scores <- updatedScores
    }
  }

  # Compile results
  RMatrixSC <- compileResults(pagodaObject, allPathways)

  return(RMatrixSC)
}
getUpdatedPathwayScores <- function(pagodaObject, pathwayData, pathway,grnSplit) {
  # Validate input types
  if (!inherits(pagodaObject, "Pagoda2")) {
    stop("pagodaObject must be a Pagoda2 object")
  }

  # Extract necessary components from pathwayData
  pcs <- pathwayData$xp
  z <- pathwayData$z

  # Compute adjusted variance
  avar <- pmax(0, (pcs$d^2 - mean(z[, 1]^2)) / sd(z[, 1]^2))

  # Extract counts and adjust by size factors
  x <- pagodaObject$counts
  sizeFactors <- pagodaObject$misc[['varinfo']][colnames(x), 'gsf']
  x@x <- x@x * rep(sizeFactors, diff(x@p))

  # Compute column means and label match
  colMeansX <- Matrix::colMeans(x)
  isTargetLabel <- colnames(x) %in% grnSplit[[pathway]]$target

  # Calculate scores
  pcs$scores <- calculateScores(x, pcs, colMeansX, isTargetLabel)

  # Flip orientations to correspond with the means
  pcs$scores <- adjustScoresOrientation(pcs$scores, x, pcs$rotation, isTargetLabel)

  return(pcs$scores)
}

#' Title
#'
#' @param x
#' @param pcs
#' @param colMeansX
#' @param isTargetLabel
#'
#'@importFrom Matrix t
#' @return
#' @export
#'
#' @examples
calculateScores <- function(x, pcs, colMeansX, isTargetLabel) {
  # Implementation of score calculation
  newScores <- as.matrix(t(x[,isTargetLabel] %*% pcs$rotation) - as.numeric((colMeansX[isTargetLabel] %*% pcs$rotation)))
  return(newScores)
}

#' Title
#'
#' @param scores
#' @param x
#' @param rotation
#' @param isTargetLabel
#'
#'@importFrom Matrix colMeans
#'
#' @return
#'
#' @examples
adjustScoresOrientation <- function(scores, x, rotation, isTargetLabel) {
  # Implementation of score orientation adjustment
  cs <- unlist(lapply(seq_len(nrow(scores)), function(i) sign(cor(scores[i,], Matrix::colMeans(t(x[, isTargetLabel, drop = FALSE])*abs(rotation[, i]))))))
  rotatedScores <- scores*cs
  return(rotatedScores)
}
compileResults <- function(pagodaObject, allPathways) {
  # Initialize an empty list for results
  resultsList <- list()

  for (pathway in allPathways) {
    pathwayScores <- pagodaObject$misc$pwpca[[pathway]]$xp$scores

    if (length(pathwayScores) == 0) {
      warning("No scores for pathway: ", pathway)
      next
    }

    resultsList[[pathway]] <- pathwayScores
  }

  # Combine results and assign row names
  RMatrixSC <- do.call(rbind, resultsList)
  rownames(RMatrixSC) <- allPathways

  return(RMatrixSC)
}


#simulator functions
getRegulonScoresSimulator <- function(seuratObject, grn_df) {
  # Validate arguments
  if (!inherits(seuratObject, "Seurat")) {
    stop("seuratObject must be a Seurat object")
  }
  if (!is.data.frame(grn_df)) {
    stop("grn_df must be a data frame")
  }

  # Access counts data
  # TODO: check if data slot is appropriate slot
  counts <- Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "data")

  # Pagoda2 object creation and initial steps
  pagodaObject <- createAndPreparePagoda2(counts)

  # Prepare GRN and create environment
  grnSplit <- prepGRN(grn_df = grn_df)
  goEnv <- prepareGRNEnvironment(grnSplit)

  # Run Pagoda2 analysis
  tryCatch({
    pagodaObject$testPathwayOverdispersion(
      goEnv, verbose = FALSE, correlation.distance.threshold = 0.2, recalculate.pca = FALSE, top.aspects = 15,
      z.score = 0
    )
  }, error = function(e) {
    message("Error in Pagoda2 analysis: ", e$message)
  })

  # Process results
  RMatrixSC <- processPagodaResults(pagodaObject,grnSplit)

  return(RMatrixSC)
}


#' Identify Significant Regulons Based on Delta Values
#'
#' This function filters regulons by identifying significant changes in their delta values compared to a random distribution.
#' It computes the density of delta values for random regulons, then uses the 2.5th and 97.5th percentiles as thresholds
#' to identify real regulons with delta values that are significantly higher or lower than these thresholds.
#'
#' @param regulon_deltas_this_cluster A data frame containing delta values for regulons,
#'        which should include columns 'class' (indicating whether the regulon is 'real' or 'random')
#'        and 'deltaPagoda' (the delta values to be analyzed).
#'
#' @return A data frame of regulons classified as 'real' that have deltaPagoda values falling outside
#'         the upper and lower thresholds defined by the random regulons' density.
#'
#' @examples
#' # Assuming 'regulon_deltas_cluster' is a data frame with the necessary structure:
#' significant_regulons <- getSignificantRegulons(regulon_deltas_cluster)
#'
#' @importFrom dplyr filter arrange
#' @importFrom stats density quantile
#' @export
getSignificantRegulons <- function(regulon_deltas_this_cluster){
  random.density <- density(subset(regulon_deltas_this_cluster, class == "random")$deltaPagoda, n = 2^10)
  upper_threshold_random <- quantile(random.density,probs = 0.975,normalize = FALSE)
  lower_threshold_random <- quantile(random.density,probs = 0.025,normalize = FALSE)

  significant_regulon_deltas_this_cluster <-  regulon_deltas_this_cluster %>%
    filter(class == "real") %>%
    filter(deltaPagoda > upper_threshold_random | deltaPagoda < lower_threshold_random)%>%
    arrange(deltaPagoda)

  return(significant_regulon_deltas_this_cluster)
}

#' Get Regulon Scores for All Clusters
#'
#' This function calculates regulon scores for each cluster in a Seurat object based on capped regulons.
#'
#' @param capped_regulons_all_clusters A list of capped regulons for each cluster, typically obtained from `capRegulonsAllClusters`.
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#'
#' @return A list where each element corresponds to a cluster and contains the regulon scores for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates the regulon scores for the capped regulons of that cluster.
#'
#' @examples
#' \dontrun{
#' capped_regulons_all_clusters <- capRegulonsAllClusters(regulons_all_clusters, decipher_seurat, TRUE)
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' regulon_scores <- getRegulonScoresAllClusters(capped_regulons_all_clusters, decipher_seurat)
#' }
#'
#' @export
getRegulonScoresAllClusters <- function(capped_regulons_all_clusters, decipher_seurat) {

  regulon_scores_all_cluster <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){
    regulon_this_cluster_capped <- capped_regulons_all_clusters[[this_cluster]]
    # main object
    decipher_seurat_this_cluster <- decipher_seurat[, which(decipher_seurat$cluster == this_cluster), seed=NULL]
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    regulon_scores_this_cluster <- getRegulonScores(
      seuratObject = decipher_seurat_this_cluster,
      grn_df = regulon_this_cluster_capped)

    regulon_scores_all_cluster[[this_cluster]] <- regulon_scores_this_cluster
  }

  return(regulon_scores_all_cluster)
}


#' Get Regulon Deltas for All Clusters
#'
#' This function calculates the differences (deltas) in regulon scores between conditions for each cluster in a Seurat object.
#'
#' @param regulon_scores_all_clusters A list of regulon scores for each cluster, typically obtained from `getRegulonScoresAllClusters`.
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#'
#' @return A list where each element corresponds to a cluster and contains the regulon deltas for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates the deltas in regulon scores for the conditions in that cluster.
#'
#' @examples
#' \dontrun{
#' regulon_scores_all_clusters <- getRegulonScoresAllClusters(capped_regulons_all_clusters, decipher_seurat)
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' regulon_deltas <- getRegulonDeltasAllClusters(regulon_scores_all_clusters, decipher_seurat)
#' }
#'
#' @export
getRegulonDeltasAllClusters <- function(regulon_scores_all_clusters, decipher_seurat) {

  regulon_deltas_all_cluster <- list()

  for(this_cluster in unique(decipher_seurat$cluster)){
    regulon_scores_this_cluster <- regulon_scores_all_clusters[[this_cluster]]
    # main object
    decipher_seurat_this_cluster <- decipher_seurat[, which(decipher_seurat$cluster == this_cluster), seed=NULL]
    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    regulon_deltas_this_cluster <- getRegulonDeltas(
      regulon_scores_this_cluster,
      decipher_seurat_this_cluster$condition)

    regulon_deltas_all_cluster[[this_cluster]] <- regulon_deltas_this_cluster
  }

  return(regulon_deltas_all_cluster)
}


#' Get Significant Regulons for All Clusters
#'
#' This function identifies significant regulon deltas for each cluster.
#'
#' @param regulon_deltas_all_clusters A list of regulon deltas for each cluster, typically obtained from `getRegulonDeltasAllClusters`.
#'
#' @return A list where each element corresponds to a cluster and contains the significant regulon deltas for that cluster.
#'
#' @details The function iterates through each cluster in the `regulon_deltas_all_clusters` list, applies the `getSignificantRegulons` function to identify significant regulon deltas, and stores the results in a list with each element corresponding to a cluster.
#'
#' @examples
#' \dontrun{
#' regulon_deltas_all_clusters <- getRegulonDeltasAllClusters(regulon_scores_all_clusters, decipher_seurat)
#' significant_regulons <- getSignificantRegulonsAllClusters(regulon_deltas_all_clusters)
#' }
#'
#' @export
getSignificantRegulonsAllClusters <- function(regulon_deltas_all_clusters) {

  significant_regulon_deltas_all_clusters <- list()

  for(this_cluster in names(regulon_deltas_all_clusters)){
    regulon_deltas_this_cluster <- regulon_deltas_all_clusters[[this_cluster]]

    significant_regulon_deltas_this_cluster <- getSignificantRegulons(regulon_deltas_this_cluster)

    significant_regulon_deltas_all_clusters[[this_cluster]] <- significant_regulon_deltas_this_cluster
  }

  return(significant_regulon_deltas_all_clusters)
}


#' Get Regulon Scores for All Clusters
#'
#' This function calculates regulon scores for each cluster in a Seurat object based on capped regulons.
#'
#' @param capped_regulons_all_clusters A list of capped regulons for each cluster, typically obtained from `capRegulonsAllClusters`.
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param paramPairings
#'
#' @return A list where each element corresponds to a cluster and contains the regulon scores for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates the regulon scores for the capped regulons of that cluster.
#'
#' @examples
#' \dontrun{
#' capped_regulons_all_clusters <- capRegulonsAllClusters(regulons_all_clusters, decipher_seurat, TRUE)
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' regulon_scores <- getRegulonScoresAllClusters(capped_regulons_all_clusters, decipher_seurat)
#' }
#'
#' @export
getRegulonScoresAllClustersWParamPairings <- function(capped_regulons_all_clusters, decipher_seurat, paramPairings) {

  regulon_scores_all_cluster <- list()

  all_clusters <- unique(decipher_seurat$cluster)

  for(this_row in c(1:nrow(paramPairings))){
    case_cluster <- paramPairings$case[this_row]
    control_cluster <- paramPairings$control[this_row]
    if(case_cluster %in% all_clusters){} else {next}
    if(control_cluster %in% all_clusters){} else {next}

    regulon_this_cluster_capped <- capped_regulons_all_clusters[[case_cluster]]

    cells <- decipher_seurat@meta.data %>%
      filter((cluster == case_cluster & condition == "case") | (cluster == control_cluster & condition == "control")) %>%
      pull(cell)

    # main object
    #decipher_seurat_this_cluster <- subset(decipher_seurat,cells = cells)
    decipher_seurat_this_cluster <- decipher_seurat[, cells, seed=NULL]


    # set identity
    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    regulon_scores_this_cluster <- getRegulonScores(
      seuratObject = decipher_seurat_this_cluster,
      grn_df = regulon_this_cluster_capped)

    regulon_scores_all_cluster[[case_cluster]] <- regulon_scores_this_cluster
  }

  return(regulon_scores_all_cluster)
}


#' Get Regulon Deltas for All Clusters
#'
#' This function calculates the differences (deltas) in regulon scores between conditions for each cluster in a Seurat object.
#'
#' @param regulon_scores_all_clusters A list of regulon scores for each cluster, typically obtained from `getRegulonScoresAllClusters`.
#' @param decipher_seurat A Seurat object containing single-cell RNA-seq data with cluster and condition metadata.
#' @param paramPairings
#'
#' @return A list where each element corresponds to a cluster and contains the regulon deltas for that cluster.
#'
#' @details The function iterates through each unique cluster in the `decipher_seurat` object, subsets the Seurat object for the cluster, normalizes the data if necessary, and calculates the deltas in regulon scores for the conditions in that cluster.
#'
#' @examples
#' \dontrun{
#' regulon_scores_all_clusters <- getRegulonScoresAllClusters(capped_regulons_all_clusters, decipher_seurat)
#' decipher_seurat <- CreateSeuratObject(counts = your_counts_matrix)
#' regulon_deltas <- getRegulonDeltasAllClusters(regulon_scores_all_clusters, decipher_seurat)
#' }
#'
#' @export
getRegulonDeltasAllClustersWParamPairings <- function(regulon_scores_all_clusters, decipher_seurat,paramPairings) {

  regulon_deltas_all_cluster <- list()

  all_clusters <- unique(decipher_seurat$cluster)

  for(this_row in c(1:nrow(paramPairings))){
    case_cluster <- paramPairings$case[this_row]
    control_cluster <- paramPairings$control[this_row]
    if(case_cluster %in% all_clusters){} else {next}
    if(control_cluster %in% all_clusters){} else {next}

    cells <- decipher_seurat@meta.data %>%
      filter((cluster == case_cluster & condition == "case") | (cluster == control_cluster & condition == "control")) %>%
      pull(cell)

    # main object
    #decipher_seurat_this_cluster <- subset(decipher_seurat,cells = cells)
    decipher_seurat_this_cluster <- decipher_seurat[, cells, seed=NULL]

    SeuratObject::Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition


    regulon_scores_this_cluster <- regulon_scores_all_clusters[[case_cluster]]
    # set identity

    if(flag.normalize.non.log){
      decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster, normalization.method = "RC", scale.factor = 100000)
    }

    regulon_deltas_this_cluster <- getRegulonDeltas(
      regulon_scores_this_cluster,
      decipher_seurat_this_cluster$condition)

    regulon_deltas_all_cluster[[case_cluster]] <- regulon_deltas_this_cluster
  }

  return(regulon_deltas_all_cluster)
}



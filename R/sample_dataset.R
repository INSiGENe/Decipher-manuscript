#' Generate a Sample Seurat Object from ExperimentHub
#'
#' This function loads a dataset from ExperimentHub, specifically searching for a dataset related to "Kang".
#' It processes the data, converts it to a Seurat object, and applies various preprocessing steps, including
#' normalization and filtering based on specified metadata criteria. It also ensures that clusters meet a
#' specified minimum number of cells per condition.
#'
#' @param min_cells_per_cluster_condition The minimum number of cells required per cluster-condition combination.
#'
#' @return A preprocessed Seurat object with clusters filtered to meet the cell count criteria.
#'
#' @examples
#' # Generate a Seurat object with a minimum of 100 cells per cluster-condition combination
#' seurat_object <- generateSampleSeuratFromExperimentHub(min_cells_per_cluster_condition = 100)
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom Seurat NormalizeData
#' @importFrom SeuratObject as.Seurat DefaultAssay
#' @importFrom dplyr filter
#' @importFrom AnnotationHub query
#' @import muscData
#' @export
generateSampleSeuratFromExperimentHub <- function(min_cells_per_cluster_condition,case_condition,control_condition){
  # Load the ExperimentHub library and create an ExperimentHub object
  eh <- ExperimentHub()
  # Search for datasets related to "Kang" in the ExperimentHub
  AnnotationHub::query(eh, "Kang")
  # Retrieve the specific dataset with ID "EH2259"
  sce <- eh[["EH2259"]]

  # Filter out cells where the total counts across all genes are zero
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]

  # Convert the SingleCellExperiment object to a Seurat object
  # 'counts' specifies which assay data to use for the initial data of the Seurat object
  kang.seurat <- Seurat::as.Seurat(sce, counts="counts", data = NULL)

  # copy 'originalexp'  to the "RNA" slot
  kang.seurat[["RNA"]] <- kang.seurat[["originalexp"]]

  # Set the default assay to "RNA" for downstream analysis
  SeuratObject::DefaultAssay(kang.seurat) <- "RNA"

  # remove 'originalexp' assay
  kang.seurat[["originalexp"]] <- NULL

  #add some fields that are missing from the conversion
  kang.seurat@meta.data$nCount_RNA <- kang.seurat@meta.data$nCount_originalexp
  kang.seurat@meta.data$nFeature_RNA <- kang.seurat@meta.data$nFeature_originalexp
  kang.seurat[["percent.mt"]] <- 0

  #this seurat is already clustered, so we don't need to do much pre-processing except some filtering
  kang.seurat <- kang.seurat[, which(kang.seurat$nFeature_RNA > 200 & kang.seurat$nFeature_RNA < 2500 & kang.seurat$percent.mt < 5), seed=NULL]

  dim(kang.seurat) #expect 18890 rows (genes) by 28869 columns (cells)

  #normalize the data
  kang.seurat <- NormalizeData(kang.seurat, normalization.method = "LogNormalize", scale.factor = 10000)

  #dataset specific meta data terms being formatted for Decipher analysis
  kang.seurat$condition <- kang.seurat$stim
  kang.seurat$orig.condition <- kang.seurat$condition
  kang.seurat$cluster <- kang.seurat$cell
  kang.seurat$original_cluster <- kang.seurat$cluster
  kang.seurat$cell <- NULL
  kang.seurat$cluster <- cleanSymbols(kang.seurat$cluster)
  kang.seurat <- KeepClustersWithMtNCellsPerCondition(kang.seurat,N = min_cells_per_cluster_condition)
  kang.seurat <- mapConditionsInSeurat(kang.seurat,"condition",case_condition,control_condition)

  return(kang.seurat)
}


#' Create standard project directory layout for a dataset
#'
#' Creates a set of commonly used directories for dataset analyses and returns a
#' named character vector of the created paths (invisibly). The function will
#' ensure the top-level `dataset_path` exists and create these subdirectories:
#' `pre_processing`, `pre_processing/h5ad_by_cluster` (returned as `co_input`),
#' `data`, `figures`, and `importances`. Existing directories are left intact.
#'
#' @param dataset_path Character(1). Path to the top-level dataset directory to create.
#'
#' @return Named character vector (invisible) with elements:
#' \describe{
#'   \item{dataset}{path to the top-level dataset directory}
#'   \item{pre_processing}{path to preprocessing directory}
#'   \item{co_input}{path to preprocessing/h5ad_by_cluster directory}
#'   \item{data}{path to data directory}
#'   \item{figures}{path to figures directory}
#'   \item{importances}{path to importances directory}
#' }
#'
#' @details
#' Directories are created with \code{dir.create(..., recursive = TRUE)} and
#' \code{showWarnings = FALSE} so the call is idempotent (no error if a directory
#' already exists). The vector of paths is returned invisibly to allow use in
#' pipeline expressions without printing.
#'
#' @examples
#' \dontrun{
#' # create the standard layout for a dataset
#' paths <- create_project_dirs("data/my_dataset")
#' # verify the figures directory exists
#' dir.exists(paths["figures"])
#' }
#'
#' @export
create_project_dirs <- function(dataset_path) {
  dirs <- c(
    dataset = dataset_path,
    pre_processing = file.path(dataset_path, "pre_processing"),
    co_input = file.path(dataset_path,'pre_processing',"h5ad_by_cluster"),
    data = file.path(dataset_path, "data"),
    figures = file.path(dataset_path, "figures"),
    importances = file.path(dataset_path, "importances")
  )

  dir.create(dataset_path, recursive = TRUE, showWarnings = FALSE)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  invisible(dirs)   # returns the paths (named) invisibly
}
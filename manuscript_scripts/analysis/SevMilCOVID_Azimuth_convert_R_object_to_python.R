#restore environment
library(renv)
renv::restore()

############
#libraries ----
############

# Load necessary library
library(Seurat)
library(dplyr)
library(babelgene)
library(devtools)
library(SummarizedExperiment)
library(zellkonverter)
library(SingleCellExperiment)
library(basilisk)
library(HGNChelper)
library(openxlsx)
library(tibble)
library(stringr)
library(org.Hs.eg.db)
library(reticulate)

############
#python env setup ----
############

    # Choose the AnnData version (e.g., latest "0.10.9")
    zellkonverter_env <- zellkonverter:::zellkonverterAnnDataEnv("0.10.9")

    #Force Basilisk to set up the environment that zellkonverter needs
    # Manually trigger environment setup
    basilisk::basiliskStart(zellkonverter_env)
    #basilisk::basiliskStop(zellkonverter_env)


#############
#functions ----
############

# Function to read and create Seurat object from file prefixes
create_seurat_object <- function(matrix_file, features_file, barcodes_file,data_dir) {
  matrix <- Read10X(data.dir = data_dir,
                    gene.column = 2,
                    barcode.column = 1)
  seurat_obj <- CreateSeuratObject(counts = matrix, project = "GSE244126_RAW")
  return(seurat_obj)
}

# Function to create directories if they don't exist
create_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

writeH5ADObjects <- function(seurat_object, pre_processing_path) {
  # Create a new directory for h5ad files if it doesn't exist
  h5ad_dir_path <- file.path(pre_processing_path, "h5ad_by_cluster")
  if (!dir.exists(h5ad_dir_path)) {
    dir.create(h5ad_dir_path)
  }

  # Process each cluster found in the Seurat object
  for (this_cluster in unique(seurat_object$cluster)) {
    #check if file already exists
    h5ad_file_path <- file.path(h5ad_dir_path, paste0(this_cluster, ".h5ad"))

    # Skip if the file already exists
    if (file.exists(h5ad_file_path)) {
      message("Skipping ", this_cluster, " — h5ad file already exists.")
      next
    }
    # Subset the Seurat object for the current cluster
    seurat_object_this_cluster <- seurat_object[,which(seurat_object$cluster == this_cluster),seed=NULL]
    # Convert to SingleCellExperiment
    sce.object <- Seurat::as.SingleCellExperiment(seurat_object_this_cluster)

    # Remove the logcounts assay if it exists
    if ("logcounts" %in% names(SummarizedExperiment::assays(sce.object))) {
      SummarizedExperiment::assays(sce.object)[["logcounts"]] <- NULL
    }

    # Write the SCE object to an h5ad file
    zellkonverter::writeH5AD(sce.object,
                             h5ad_file_path,
                             X_name = "counts")
  }
}

cleanSymbols <- function(string) {
  # Remove or replace various symbols with safer alternatives
  string <- stringr::str_remove_all(string, "/")
  string <- stringr::str_replace_all(string, " ", "_")
  string <- stringr::str_replace_all(string, "\\+", "_plus_")
  string <- stringr::str_replace_all(string, "\\-", "_minus_")
  string <- stringr::str_replace_all(string, "\\(", "")
  string <- stringr::str_replace_all(string, "\\)", "")
  string <- stringr::str_replace_all(string, "%", "_percent_")
  string <- stringr::str_replace_all(string, "\\.", "_dot_")
  string <- stringr::str_replace_all(string, ",", "_comma_")
  string <- stringr::str_replace_all(string, ":", "_colon_")
  string <- stringr::str_replace_all(string, ";", "_semicolon_")
  string <- stringr::str_replace_all(string, "&", "_and_")
  string <- stringr::str_replace_all(string, "\\?", "_question_")
  string <- stringr::str_replace_all(string, "!", "_exclamation_")
  string <- stringr::str_replace_all(string, "\"", "_quote_")
  string <- stringr::str_replace_all(string, "'", "_apostrophe_")
  string <- stringr::str_replace_all(string, "=", "_equals_")
  string <- stringr::str_replace_all(string, "\\*", "_asterisk_")
  string <- stringr::str_replace_all(string, "#", "_hash_")
  string <- stringr::str_replace_all(string, "@", "_at_")
  string <- stringr::str_replace_all(string, "\\$", "_dollar_")
  string <- stringr::str_replace_all(string, "\\^", "_caret_")
  string <- stringr::str_replace_all(string, "<", "_less_than_")
  string <- stringr::str_replace_all(string, ">", "_greater_than_")
  string <- stringr::str_replace_all(string, "\\[", "_lbracket_")
  string <- stringr::str_replace_all(string, "\\]", "_rbracket_")
  string <- stringr::str_replace_all(string, "\\{", "_lbrace_")
  string <- stringr::str_replace_all(string, "\\}", "_rbrace_")
  string <- stringr::str_replace_all(string, "\\|", "_pipe_")
  string <- stringr::str_replace_all(string, "\\\\", "_backslash_")
  string <- stringr::str_replace_all(string, "/", "_slash_")
  string <- stringr::str_replace_all(string, "__", "_")  # Remove double underscores
  return(string)
}

get_column <- function(df, possible_names) {
  found <- intersect(possible_names, colnames(df))
  if (length(found) == 0) {
    stop("None of the expected column names found: ", paste(possible_names, collapse = ", "))
  }
  df[[found[1]]]
}


set.seed(123)
h5ad_path <- "data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.h5ad"
if (!file.exists(h5ad_path)) {
  seurat_object_oi <- readRDS("data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")
  sce.object = as.SingleCellExperiment(seurat_object_oi)
  sce.object@assays@data[["logcounts"]] <- NULL

 # writeH5AD(
 #   sce.object,
 #   h5ad_path,
 #   X_name = "counts"    )
}else {
  message("h5ad file already exists at ", h5ad_path, ", skipping.")
}

message("conversion completed")


#Write counts (genes × cells)
counts_mat <- as.matrix(assay(sce.object, "counts"))
write.csv(counts_mat, file = "counts.csv", row.names = TRUE)

# Write observation metadata (cells × covariates)
obs_df <- as.data.frame(colData(sce.object)[, c("sample_id", "severity_group", "predicted.celltype.l2", "orig.ident")])
write.csv(obs_df, file = "obs.csv", row.names = TRUE)

# Write variable metadata (genes × attributes, optional)
var_df <- data.frame(gene_id = rownames(sce.object))
write.csv(var_df, file = "var.csv", row.names = FALSE)


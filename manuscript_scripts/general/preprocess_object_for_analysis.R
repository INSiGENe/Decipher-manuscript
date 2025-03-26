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

############
#parameters ----
############

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a dataset key, e.g. 'cz_rcc'")
}
dataset_key <- args[1]

config <- jsonlite::fromJSON(txt = "scripts/config.json")

if (!dataset_key %in% names(config)) {
  stop(paste("Dataset key not found in config:", dataset_key))
}

cfg <- config[[dataset_key]][["pre_processing_for_analysis"]]

if (is.null(cfg)) {
  stop(paste("No 'pre_processing_for_analysis' section found for dataset:", dataset_key))
}


############
#data and analysis ----
############

#set a seed
set.seed(123)


#directories
input_path <- cfg$input_path
results_path <- cfg$output_path
pre_processing_path <- file.path(results_path, "pre_processing")
# Path to cached mapping file
mapping_path <- file.path(pre_processing_path, "ensembl_to_hgnc_mapping.tsv")


case_condition <- cfg$case_condition
control_condition <- cfg$control_condition

cytosig_path <- file.path(results_path,"cytosig")

liana_filepath <- file.path(results_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(results_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(cytosig_path,recursive=TRUE)
dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)
dir.create(pre_processing_path,recursive=TRUE)

############
#Decipher pre-processing ----
############
seurat_object_rds_path <- file.path(pre_processing_path, "seurat_object_oi.rds")

if (!file.exists(seurat_object_rds_path)) {
  #if the input comes from a pre-processed anndata object I want to have a flag in the config.json file that takes that into account
  matrix <- Matrix::readMM(file.path(input_path, "output_matrix.mtx"))
  matrix <- Matrix::t(matrix)

  gene_names <- read.csv(file.path(input_path, "gene_names.csv"))
  meta_data <- read.csv(file.path(input_path, "filtered_metadata.csv"), header = TRUE)
  cell_names <- read.csv(file.path(input_path, "cell_names.csv"))

  # Use dynamic column selection
  rownames(matrix) <- get_column(gene_names, c("X0", "name", "gene", "gene_name","index","gene_ids"))
  colnames(matrix) <- get_column(cell_names, c("X0", "name", "barcode", "cell_id","index","cellId"))
  rownames(meta_data) <- get_column(meta_data, c("X", "barcode","cellId"))

  # Convert Ensembl IDs to HGNC gene symbols
  # Generate or load mapping
  if (file.exists(mapping_path)) {
    message("Using existing Ensembl → HGNC mapping file: ", mapping_path)
    mapping <- read.delim(mapping_path)
  } else {
    message("Generating Ensembl → HGNC mapping...")
    ensembl_ids <- rownames(matrix)
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ensembl_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    )

    mapping <- mapping[!is.na(mapping$SYMBOL) & !duplicated(mapping$SYMBOL), ]

    write.table(mapping, mapping_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # Subset and update matrix
  matrix <- matrix[mapping$ENSEMBL, ]
  rownames(matrix) <- mapping$SYMBOL


  meta_data <- meta_data[,-1]

  meta_data$cluster <- meta_data[[cfg$cluster_meta_field]]
  meta_data$cluster <- cleanSymbols(meta_data$cluster)

  meta_data$condition <- meta_data[[cfg$condition_meta_field]]
  meta_data$sample_key <- paste(meta_data[[cfg$individual_meta_field]], meta_data$condition, sep="_")

  seurat_object_oi <- Seurat::CreateSeuratObject(counts = matrix,meta.data = meta_data)

  seurat_object_oi@meta.data <- seurat_object_oi@meta.data %>%
    mutate(condition_original = if_else(condition == case_condition,"case","control"),
          cluster = cluster)

  saveRDS(seurat_object_oi,seurat_object_rds_path)
} else {
  message("Seurat object already exists at ", seurat_object_rds_path, ", skipping creation.")
}

seurat_object_oi <- readRDS(seurat_object_rds_path)
writeH5ADObjects(seurat_object_oi,pre_processing_path)

############
# CYTOSIG pre-processing ----
############

seurat_object_oi <- readRDS(seurat_object_rds_path)

Idents(seurat_object_oi) <- seurat_object_oi$cluster

for(this_cluster in unique(seurat_object_oi$cluster)){

  cytosig_cluster_path <- file.path(cytosig_path,this_cluster)
  dir.create(cytosig_cluster_path,recursive = TRUE)

  output_file <- file.path(cytosig_cluster_path, "differential_profile.tsv.gz")

  if (file.exists(output_file)) {
    message("Skipping ", this_cluster, " — output already exists.")
    next
  }

  # Assume seurat_object_oi is a Seurat object
  seurat_object_oi_this_cluster <- subset(seurat_object_oi,idents = this_cluster)

  # Filter genes with at least 10 counts
  #seurat_object_oi_this_cluster <- subset(seurat_object_oi_this_cluster, subset = nCount_RNA >= 10)
  # Calculate the total counts per gene
  gene_counts <- rowSums(seurat_object_oi_this_cluster@assays$RNA@counts)

  # keep genes with counts (as per cytosig vignette)
  genes_to_keep <- names(gene_counts[gene_counts > 0])

  # Subset the Seurat object to keep only these genes
  seurat_object_oi_this_cluster <- subset(seurat_object_oi_this_cluster, features = genes_to_keep)

  # Normalize the total counts to a target sum of 1e5
  count_data <- seurat_object_oi_this_cluster@assays$RNA@counts
  total_counts <- colSums(count_data)
  normalized_counts <- t(t(count_data) / total_counts * 1e5)

  # Apply logarithmic transformation with base 2 to the gene expression data
  log_transformed_counts <- log2(normalized_counts+1)
  rm(normalized_counts)

  # Extract the data for the control condition ("NE")
  control_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == control_condition)
  control_data <- log_transformed_counts[,control_cells]

  # Calculate the mean expression for each gene
  control_mean_expression <- rowMeans(control_data)
  rm(control_data)
  gc()

  # Extract the data for the experimental condition ("E")
  case_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == case_condition)
  case_data <- log_transformed_counts[,case_cells]

  # Calculate the differential expression profile by subtracting the control mean from the experimental data
  differential_profile <- case_data-control_mean_expression
  rm(case_data)
  rm(log_transformed_counts)
  gc()

  # Write the differential profile to a tab-separated .gz file
  gz1 <- gzfile(file.path(cytosig_cluster_path,"differential_profile.tsv.gz"), "w")
  write.table(differential_profile, gz1, sep = "\t", col.names = NA, quote = FALSE)
  close(gz1)

  rm(differential_profile)
  gc()

}

############
# LIANA pre-processing ----
############

liana_h5ad_path <- file.path(liana_data_filepath, "seurat_object_oi.h5ad")
if (!file.exists(liana_h5ad_path)) {
  seurat_object_oi <- readRDS(seurat_object_rds_path)
  seurat_object_oi$condition <- seurat_object_oi@meta.data %>%
    mutate(liana_condition = if_else(condition == case_condition,"stim","ctrl"))%>%
    dplyr::select(liana_condition)
  sce.object = as.SingleCellExperiment(seurat_object_oi)
  sce.object@assays@data[["logcounts"]] <- NULL
  writeH5AD(sce.object, liana_h5ad_path,X_name = "counts")
  rm(sce.object)
}else {
  message("LIANA h5ad file already exists at ", liana_h5ad_path, ", skipping.")
}

############
# NATMI pre-processing ----
############

#load data
seurat_object_oi = readRDS(seurat_object_rds_path)
DefaultAssay(seurat_object_oi) <- "RNA"
seurat_object_oi <- Seurat::NormalizeData(seurat_object_oi)

#data pre-processing
Idents(seurat_object_oi) <- seurat_object_oi$cluster

natmi_case_em_path <- file.path(natmi_data_filepath, "case", "em.txt")
natmi_case_metadata_path <- file.path(natmi_data_filepath, "case/metadata.txt")
if (!file.exists(natmi_case_em_path) || !file.exists(natmi_case_metadata_path)) {
  message("Preparing NATMI input for case condition...")
  
  #data pre-processing
  dir.create(file.path(natmi_data_filepath,"case"),recursive = TRUE)
  seurat_object_oi_subset <- seurat_object_oi[, which(seurat_object_oi$condition == case_condition),seed=NULL]
  if (!file.exists(natmi_case_em_path)){
    data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
    data_matrix <- expm1(data_matrix)
    data_matrix <- 100*data_matrix

    write.table(data_matrix, file.path(natmi_data_filepath,"case/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))

    rm(data_matrix)
    gc()
  }
  

  meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var = "barcode")

  # Explicitly drop 'annotation' if it exists
  if ("annotation" %in% colnames(meta_data)) {
    meta_data <- meta_data[, setdiff(colnames(meta_data), "annotation")]
  }

  # Now rename and select
  meta_data <- meta_data %>%
    dplyr::rename(annotation = cluster) %>%
    dplyr::select(barcode, annotation)

  write.table(meta_data,natmi_case_metadata_path, quote = F,sep="\t",row.names=FALSE,col.names=TRUE)
}

natmi_control_em_path <- file.path(natmi_data_filepath, "control", "em.txt")
natmi_control_metadata_path <- file.path(natmi_data_filepath, "control/metadata.txt")
if (!file.exists(natmi_control_em_path) || !file.exists(natmi_control_metadata_path)) {
  message("Preparing NATMI input for control condition...")
  
  #data pre-processing
  dir.create(file.path(natmi_data_filepath, "control"),recursive = TRUE)
  seurat_object_oi_subset <- seurat_object_oi[, which(seurat_object_oi$condition == control_condition),seed=NULL]
  if (!file.exists(natmi_control_em_path)){
    data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
    data_matrix <- expm1(data_matrix)
    data_matrix <- 100*data_matrix

    write.table(data_matrix, natmi_control_em_path, quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))

    rm(data_matrix)
    gc()
  }
  
  meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var = "barcode")

  # Explicitly drop 'annotation' if it exists
  if ("annotation" %in% colnames(meta_data)) {
    meta_data <- meta_data[, setdiff(colnames(meta_data), "annotation")]
  }

  # Now rename and select
  meta_data <- meta_data %>%
    dplyr::rename(annotation = cluster) %>%
    dplyr::select(barcode, annotation)

  write.table(meta_data,natmi_control_metadata_path, quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

}

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")



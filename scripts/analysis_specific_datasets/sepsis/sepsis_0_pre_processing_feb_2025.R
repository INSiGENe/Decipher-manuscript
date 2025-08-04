############
#package installation----
############

renv::restore()
#say yes and 1 and yes

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
#library(tidyr)
library(data.table)
#library(tidyverse)

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

writeH5ADObjects <- function(seurat_object, pre_processing_path) {
  # Create a new directory for h5ad files if it doesn't exist
  h5ad_dir_path <- file.path(pre_processing_path, "h5ad_by_cluster")
  if (!dir.exists(h5ad_dir_path)) {
    dir.create(h5ad_dir_path)
  }

  # Process each cluster found in the Seurat object
  for (this_cluster in unique(seurat_object$cluster)) {
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
                             file.path(h5ad_dir_path, paste(this_cluster, ".h5ad", sep = "")),
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


############
#data and analysis ----
############

#set a seed
set.seed(123)

#directories
dataset_path <- "data/sepsis"
results_path <- "data/TNBC"
pre_processing_path <- file.path(results_path,"pre_processing")
case_condition <- "ICU-SEP"
control_condition <- "ICU-NoSEP"
cytosig_path <- file.path(results_path,"cytosig")
liana_filepath <- file.path(results_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(results_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)
dir.create(pre_processing_path,recursive=TRUE)


# RAW DATA PRE-PROCESSING to Seurat ----
# please download and unzip scp_gex_matrix_raw.csv.gz
# please download scp_meta_updated.txt
# from Broad Institute Single Cell Portal SCP548
# place files in pre_processing_path
# print (pre_processing_path)
# Initialize lists to store sparse matrices and gene names
sparse_matrix_list <- list()
gene_names_list <- list()

# Read the first 2000 rows of the raw matrix CSV file
matrix = fread(file.path(dataset_path,"scp_gex_matrix_raw.csv"), header = TRUE, nrows = 2000, skip = 0)

# Extract gene names from the first column
gene_names <- matrix[, 1]$GENE

# Remove the first column (gene names) from the matrix
matrix <- matrix[, -1]

# Convert the matrix to a sparse matrix
matrix <- as.matrix(matrix)
matrix <- as(matrix, "sparseMatrix")

# Store the sparse matrix and gene names in their respective lists
sparse_matrix_list[[1]] <- matrix
gene_names_list[[1]] <- gene_names

# Loop through the remaining chunks of the raw matrix CSV file
for (i in 2:12) {
  matrix = fread(file.path(dataset_path,"scp_gex_matrix_raw.csv"), header = FALSE, nrows = 2000, skip = (i - 1) * 2000)
  gene_names <- matrix[, 1]$V1
  matrix <- matrix[, -1]
  matrix <- as.matrix(matrix)
  matrix <- as(matrix, "sparseMatrix")
  sparse_matrix_list[[i]] <- matrix
  gene_names_list[[i]] <- gene_names
}

# Initialize the main sparse matrix and gene names vector with the first chunk
main_sparse_matrix <- sparse_matrix_list[[1]]
gene_names_vector <- gene_names_list[[1]]

# Append the remaining chunks to the main sparse matrix and gene names vector
for (i in 2:12) {
  main_sparse_matrix <- rbind(main_sparse_matrix, sparse_matrix_list[[i]])
  gene_names_vector <- c(gene_names_vector, gene_names_list[[i]])
}

# Set the row names of the main sparse matrix to the gene names
rownames(main_sparse_matrix) <- gene_names_vector

# Save the main sparse matrix to an RDS file
saveRDS(main_sparse_matrix, file.path(pre_processing_path,"main_sparse_matrix.rds"))

# Read metadata from a text file
meta_data <- read.table(file.path(dataset_path,"scp_meta_updated.txt"), quote = "\"", sep = "\t", header = TRUE, skip = 0)

# Set the row names of the metadata to the NAME column
rownames(meta_data) <- meta_data$NAME

# Display the dimensions of the metadata
dim(meta_data)

# Create a Seurat object using the main sparse matrix and metadata
seurat_object <- Seurat::CreateSeuratObject(counts = main_sparse_matrix, meta.data = meta_data)

# Normalize the data in the Seurat object
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Save the Seurat object to an RDS file
saveRDS(seurat_object, file.path(pre_processing_path,"seurat_object.rds"))

# DECIPHER PRE-PROCESSING ----
## load data
seurat_object <- readRDS(file.path(pre_processing_path,"seurat_object.rds"))
dim(seurat_object)
table(seurat_object$Cohort)


##subset to cohorts of interest
seurat_object_oi <- subset(seurat_object,subset = Cohort %in% c(case_condition,control_condition))
#seurat_object <- subset(seurat_object,subset = biosample_id %in% c("CD45"))

seurat_object_oi@meta.data <- seurat_object_oi@meta.data %>%
  mutate(condition_original = if_else(Cohort == case_condition,"case","control"),
         condition = Cohort,
         cluster = paste(Cell_Type,Cell_State))

seurat_object_oi$cluster <- cleanSymbols(seurat_object_oi$cluster)

saveRDS(seurat_object_oi,file.path(pre_processing_path,"seurat_object_oi.rds"))
writeH5ADObjects(seurat_object_oi,pre_processing_path)

# CYTOSIG pre-processing ----
seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))
Idents(seurat_object_oi) <- seurat_object_oi$cluster

for(this_cluster in unique(seurat_object_oi$cluster)){

  cytosig_cluster_path <- file.path(cytosig_path,this_cluster)
  dir.create(cytosig_cluster_path,recursive = TRUE)

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

  # Extract the data for the control condition ("NE")
  control_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == control_condition)
  control_data <- log_transformed_counts[,control_cells]

  # Calculate the mean expression for each gene
  control_mean_expression <- rowMeans(control_data)

  # Extract the data for the experimental condition ("E")
  case_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == case_condition)
  case_data <- log_transformed_counts[,case_cells]

  # Calculate the differential expression profile by subtracting the control mean from the experimental data
  differential_profile <- case_data-control_mean_expression

  # Write the differential profile to a tab-separated .gz file
  gz1 <- gzfile(file.path(cytosig_cluster_path,"differential_profile.tsv.gz"), "w")
  write.table(as.matrix(differential_profile), gz1, sep = "\t", col.names = NA, quote = FALSE)
  close(gz1)
}

# LIANA pre-processing ----
seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

seurat_object_oi$condition <- seurat_object_oi@meta.data %>%
  mutate(liana_condition = if_else(condition == case_condition,"stim","ctrl"))%>%
  select(liana_condition)
sce.object = as.SingleCellExperiment(seurat_object_oi)
sce.object@assays@data[["logcounts"]] <- NULL
writeH5AD(sce.object, file.path(liana_data_filepath,"seurat_object_oi.h5ad"),X_name = "counts")

# NATMI pre-processing ----
#load data
seurat_object_oi = readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))
DefaultAssay(seurat_object_oi) <- "RNA"
seurat_object_oi <- Seurat::NormalizeData(seurat_object_oi)

#data pre-processing
Idents(seurat_object_oi) <- seurat_object_oi$cluster

#data pre-processing
dir.create(file.path(natmi_data_filepath,"case"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == case_condition)


data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- exp(as.matrix(data_matrix))- 1

write.table(100 * data_matrix, file.path(natmi_data_filepath,"case/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=TRUE)
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  dplyr::rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"case/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"control"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == control_condition)

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- exp(as.matrix(data_matrix)) - 1
write.table(100 * data_matrix, file.path(natmi_data_filepath,"control/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=TRUE)
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  dplyr::rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"control/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")

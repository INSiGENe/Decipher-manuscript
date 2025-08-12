############
#package installation----
############

#renv::restore()
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

dataset_path <- "data/TNBC"
results_path <- "results/TNBC"
pre_processing_path <- file.path(results_path,"pre_processing")
case_condition <- "E" #expanding
control_condition <- "NE" #non-expanding
cytosig_path <- file.path(results_path,"cytosig")
liana_filepath <- file.path(results_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(results_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)
dir.create(pre_processing_path,recursive=TRUE)

cohort_1_count_matrix <- readRDS(file.path(dataset_path,"1863-counts_cells_cohort1.rds"))
cohort_1_metadata <- read.csv(file.path(dataset_path,"1872-BIOKEY_metaData_cohort1_web.csv.xls"))
cells_of_interest <- cohort_1_metadata %>%
  filter(
    timepoint == "On",
    BC_type == "TNBC",
    expansion %in% c("E","NE")
  ) %>%
  select(Cell) %>%
  unlist()

rownames(cohort_1_metadata) <- cohort_1_metadata$Cell
cohort_1_metadata <- cohort_1_metadata[,-1]
seurat_object <- Seurat::CreateSeuratObject(counts = cohort_1_count_matrix[,cells_of_interest],meta.data = cohort_1_metadata[cells_of_interest,])
rm(cohort_1_count_matrix)
# Normalize the data in the Seurat object
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(seurat_object, file.path(pre_processing_path,"seurat_object.rds"))
rm(seurat_object)
# DECIPHER PRE-PROCESSING ----

seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object.rds"))
seurat_object_oi@meta.data <- seurat_object_oi@meta.data %>%
  mutate(condition_original = if_else(expansion == case_condition,"case","control"),
         condition = expansion,
         cluster = paste(cellType))

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
data_matrix <- expm1(data_matrix)
write.table(100 * data_matrix, file.path(natmi_data_filepath,"case/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  dplyr::rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"case/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"control"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == control_condition)

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- expm1(data_matrix)
write.table(100 * data_matrix, file.path(natmi_data_filepath,"control/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  dplyr::rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"control/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")

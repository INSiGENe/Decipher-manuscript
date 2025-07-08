#L3 cohort refers to the cohort of cases and controls that are age matched and equal in number of cases of Asian and European ancestry in processing batch 4 and their replicates in other batches. Processing batch 4 refers to the L3 cohort samples within processing batch 4 only.

#libraries ----
# library(dplyr)
# library(tidyr)
# library(basilisk)
# library(zellkonverter)
library(devtools)
load_all()

#global settings ----
set.seed(123)

#create a case matrix
dataset_path <- ("manuscript_analysis/lupus")
pre_processing_path <- file.path(dataset_path,"pre_processing")
case_condition <- "Managed" #expanding
control_condition <- "Healthy" #non-expanding
cytosig_path <- file.path(dataset_path,"cytosig")
liana_filepath <- file.path(dataset_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(dataset_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)

#Raw data pre-processing ----
#please run lupus_0a_pre_processing.py first on the raw GEO h5ad objects
matrix <- Matrix::readMM(file.path(pre_processing_path,"output_matrix.mtx"))

matrix <- Matrix::t(matrix)
gene_names <- read.csv(file.path(pre_processing_path,"gene_names.csv"))
meta_data <- read.csv(file.path(pre_processing_path,"filtered_metadata.csv"),header = TRUE)
rownames(matrix) <- gene_names$X0
colnames(matrix) <- meta_data$X
rownames(meta_data) <- meta_data$X
meta_data <- meta_data[,-1]

cg_mapping <- c(
  "B" = "B cell",
  "cDC" = "cDC",
  "cM"="Classical Monocyte",
  "ncM"="Non-classical Monocyte",
  "NK"="NK",
  "PB" = "Plasmablast",
  "pDC"="pDC",
  "Progen"="Progenitor cell",
  "Prolif"="Proliferating cell",
  "T4"="CD4 T",
  "T8"="CD8 T"
)

meta_data$cluster <- cg_mapping[meta_data$cg_cov]
meta_data$cluster <- cleanSymbols(meta_data$cluster)

meta_data$condition <- meta_data$Status
meta_data$sample_key <- paste(meta_data$ind_cov,meta_data$condition,sep="_")

seurat_object <- Seurat::CreateSeuratObject(counts = matrix,meta.data = meta_data)
saveRDS(seurat_object,file.path(pre_processing_path,"seurat_object.rds"))

#Decipher data pre-processing ----
seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object.rds"))
seurat_object_oi@meta.data <- seurat_object_oi@meta.data %>%
  mutate(condition_original = if_else(Status == case_condition,"case","control"),
         condition = Status,
         cluster = cluster)

seurat_object_oi$cluster <- cleanSymbols(seurat_object_oi$cluster)

saveRDS(seurat_object_oi,file.path(pre_processing_path,"seurat_object_oi.rds"))
writeH5ADObjects(seurat_object,pre_processing_path)

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


  #try to convert to sparse matrix before saving

  # Write the differential profile to a tab-separated .gz file
  gz1 <- gzfile(file.path(cytosig_cluster_path,"differential_profile.tsv.gz"), "w")
  write.table(differential_profile, gz1, sep = "\t", col.names = NA, quote = FALSE)
  close(gz1)

  rm(differential_profile)
  gc()

}

# LIANA pre-processing ----
seurat_object_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

seurat_object_oi$condition <- seurat_object_oi@meta.data %>%
  mutate(liana_condition = if_else(condition == case_condition,"stim","ctrl"))%>%
  select(liana_condition)
sce.object = as.SingleCellExperiment(seurat_object_oi)
sce.object@assays@data[["logcounts"]] <- NULL
writeH5AD(sce.object, file.path(liana_data_filepath,"seurat_object_oi.h5ad"),X_name = "counts")
rm(sce.object)
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

set.seed(123)
if(dim(seurat_object_oi_subset)[2] > 50000){
  all_cells <- colnames(seurat_object_oi_subset)
  rand_cells <- sample(all_cells,size=50000)
  seurat_object_oi_subset<- subset(seurat_object_oi_subset,cells = rand_cells)
}

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- expm1(data_matrix)
data_matrix <- 100*data_matrix
write.table(data_matrix, file.path(natmi_data_filepath,"case/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"case/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"control"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == control_condition)

#optional if not enough ram
set.seed(123)
if(dim(seurat_object_oi_subset)[2] > 50000){
  all_cells <- colnames(seurat_object_oi_subset)
  rand_cells <- sample(all_cells,size=50000)
  seurat_object_oi_subset<- subset(seurat_object_oi_subset,cells = rand_cells)
}

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- expm1(data_matrix)
data_matrix <- 100*data_matrix
gc()
write.table( data_matrix, file.path(natmi_data_filepath,"control/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"control/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")



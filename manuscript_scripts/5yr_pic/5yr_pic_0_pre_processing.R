library(devtools)
load_all()

#set a seed
set.seed(123)

#directories
dataset_path <- "manuscript_analysis/5yr_pic"
pre_processing_path <- file.path(dataset_path,"pre_processing")
case_condition <- "PIC"
control_condition <- "CTRL"
cytosig_path <- file.path(dataset_path,"cytosig")
liana_filepath <- file.path(dataset_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(dataset_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)

# DECIPHER PRE-PROCESSING ----
## load data
seurat_object <- readRDS(file.path(pre_processing_path,"seurat_object_integrated_annotated.rds"))
DefaultAssay(seurat_object) <- "RNA"
seurat_object[["integrated"]] <- NULL

cells_to_keep <-  seurat_object@meta.data %>%
  filter(Stimuli %in% c("CTRL", "PIC") & Age == "5yr_PBMC") %>%
  rownames()

matrix <- attr(seurat_object[['RNA']],"layers")$counts
cells <- attr(seurat_object[['RNA']],"cells")
genes <- attr(seurat_object[['RNA']],"features")
colnames(matrix) <- rownames(cells@.Data)
rownames(matrix) <- rownames(genes@.Data)

matrix_oi <- matrix[,cells_to_keep]
meta_data_oi <- seurat_object@meta.data[cells_to_keep,]

seurat_object_oi <- CreateSeuratObject(counts = matrix_oi,meta.data = meta_data_oi)
rm(seurat_object)
rm(matrix_oi)
rm(matrix)
rm(meta_data_oi)

seurat_object_oi@meta.data <- seurat_object_oi@meta.data %>%
  mutate(condition_original = if_else(Stimuli == case_condition,"case","control"),
         condition = Stimuli,
         cluster = predicted.celltype.l1)

seurat_object_oi$cluster <- cleanSymbols(seurat_object_oi$cluster)

CpC_data <- generateQCDataByClusterAndCondition(seurat_object_oi,max(stringr::str_length(unique(seurat_object_oi$cluster))))
CpC_data

seurat_object_oi <- subset(seurat_object_oi,subset=cluster != c("CD8_T"))
seurat_object_oi <- subset(seurat_object_oi,subset=cluster != c("DC"))
seurat_object_oi <- subset(seurat_object_oi,subset=cluster != c("other"))
seurat_object_oi <- subset(seurat_object_oi,subset=cluster != c("other_T"))


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
  if(length(which(seurat_object_oi_this_cluster$condition == case_condition)) == 0 || length(which(seurat_object_oi_this_cluster$condition == control_condition)) == 0 ){
      next
  }

  # Assume seurat_object_oi is a Seurat object
  seurat_object_oi_this_cluster <- subset(seurat_object_oi,idents = this_cluster)

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
rm(case_data)
rm(control_data)
rm(count_data)
rm(normalized_counts)
rm(differential_profile)
rm(log_transformed_counts)

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
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"case/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"control"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == control_condition)

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- exp(as.matrix(data_matrix)) - 1
write.table(100 * data_matrix, file.path(natmi_data_filepath,"control/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=TRUE)
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"control/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")

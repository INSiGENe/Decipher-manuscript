# Load necessary library
library(Seurat)
library(dplyr)
library(babelgene)
library(devtools)
load_all()

#set a seed
set.seed(123)

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


#directories
dataset_path <- "data/BCG"
results_path <- "results/BCG"
pre_processing_path <- file.path(results_path,"pre_processing")
tenx_path <- file.path(dataset_path,"GSE244126_RAW")
case_condition <- "D21" #vaccinated
control_condition <- "D0" #unvaccinated
cytosig_path <- file.path(results_path,"cytosig")
liana_filepath <- file.path(results_path,"liana")
liana_data_filepath <- file.path(liana_filepath,"data")
natmi_filepath <- file.path(results_path,"natmi")
natmi_data_filepath <- file.path(natmi_filepath,"data")

dir.create(natmi_data_filepath,recursive=TRUE)
dir.create(liana_data_filepath,recursive=TRUE)


# List of prefixes for the different conditions
prefixes <- c("GSM7807682_Lung_D0",
              "GSM7807683_Lung_D21")

prefixes <- file.path(tenx_path,prefixes)

# Loop through each prefix and move corresponding files
for (prefix in prefixes) {

  # Create a directory for the prefix
  create_dir(prefix)

  # Define the file patterns to match
  file_patterns <- c("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz")

  # Move the files to the new directory
  for (pattern in file_patterns) {
    file_name <- paste0(prefix, pattern)
    base_name <- basename(file_name)
    file.rename(file_name, file.path(prefix, gsub("_", "", pattern)))
  }
}
# Initialize list to store Seurat objects
seurat_objects <- list()

# Loop through each prefix and create Seurat objects
for (tenx_data in prefixes) {
  base_name <- basename(tenx_data)
  matrix <- Read10X(data.dir = tenx_data,
                    gene.column = 2,
                    cell.column = 1)

  # Use the function to convert the gene symbols in the Seurat object from mouse to human
  src_genes <- rownames(matrix)
  orthologs <- babelgene::orthologs(src_genes, species = "mouse",human = FALSE, top = TRUE)

  orthologs_most_support <- orthologs %>%
    # Remove rows where symbol or human_symbol is NA
    filter(!is.na(symbol) & !is.na(human_symbol)) %>%
    group_by(symbol) %>%
    filter(support_n == max(support_n)) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(human_symbol) %>%
    filter(support_n == max(support_n)) %>%
    slice(1) %>%
    ungroup()

  matrix_w_orthologs <- matrix[orthologs_most_support$symbol,]
  rownames(matrix_w_orthologs) <- orthologs_most_support$human_symbol

  seurat_obj <- CreateSeuratObject(counts = matrix_w_orthologs, project = "GSE232186_RAW",min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  seurat_objects[[base_name]] <- seurat_obj
}

seurat_objects$GSM7807682_Lung_D0 <- RenameCells(seurat_objects$GSM7807682_Lung_D0, new.names = paste(colnames(seurat_objects$GSM7807682_Lung_D0),"1",sep="_"))
seurat_objects$GSM7807682_Lung_D0@meta.data$timepoint <- "D0"

seurat_objects$GSM7807683_Lung_D21 <- RenameCells(seurat_objects$GSM7807683_Lung_D21, new.names = paste(colnames(seurat_objects$GSM7807683_Lung_D21),"2",sep="_"))
seurat_objects$GSM7807683_Lung_D21@meta.data$timepoint <- "D21"


# Merge all Seurat objects into one, if necessary
merged_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1], project = "GSE244126_RAW")

#sSeurat
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
# merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # make some filtering based on QC metrics visualizations, see Seurat tutorial: https://satijalab.org/seurat/articles/merged_seurat3k_tutorial.html
merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
merged_seurat <- ScaleData(merged_seurat, features = rownames(merged_seurat))
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

# Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(merged_seurat)

# cluster and visualize
#elbow plot suggests around 15
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:15)
#suggested between 0.4 to 1.2 for sc datasets around 3k (increases for larger datasets)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.4)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:15)
DimPlot(merged_seurat, reduction = "umap")


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus

# prepare gene sets
library(HGNChelper)
gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(merged_seurat[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(merged_seurat[["RNA"]]$scale.data) else as.matrix(merged_seurat[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. merged_seurat[["RNA"]]$counts or merged_seurat[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(merged_seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(merged_seurat@meta.data[merged_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
saveRDS(sctype_scores,file.path(pre_processing_path,"immune_system_sctype_cluster_scores.rds"))


tissue <- "Lung" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus

# prepare gene sets
library(HGNChelper)
gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(merged_seurat[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(merged_seurat[["RNA"]]$scale.data) else as.matrix(merged_seurat[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. merged_seurat[["RNA"]]$counts or merged_seurat[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(merged_seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(merged_seurat@meta.data[merged_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
saveRDS(sctype_scores,file.path(pre_processing_path,"lung_sctype_cluster_scores.rds"))

lung_sctype_cluster_scores <- readRDS(file.path(pre_processing_path,"lung_sctype_cluster_scores.rds"))
immune_sctype_cluster_scores <- readRDS(file.path(pre_processing_path,"immune_system_sctype_cluster_scores.rds"))

sctype_cluster_scores <- left_join(immune_sctype_cluster_scores,lung_sctype_cluster_scores,by = "cluster")
print(sctype_cluster_scores,n=nrow(sctype_cluster_scores))


# Assuming sctype_cluster_scores is your data frame
# Create a named vector with cell type labels having the highest scores
ct_labels <- sctype_cluster_scores %>%
  rowwise() %>%
  mutate(max_score = if_else(scores.x >= scores.y, scores.x, scores.y),
         cell_type = if_else(scores.x >= scores.y, type.x, type.y)) %>%
  ungroup() %>%
  select(cluster, cell_type)

# Print the resulting named vector
print(ct_labels)

ind_match <- match(merged_seurat$seurat_clusters,ct_labels$cluster)
merged_seurat$cluster <- ct_labels$cell_type[ind_match]
table(merged_seurat$cluster)

# DECIPHER PRE-PROCESSING ----
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  mutate(original_condition = if_else(timepoint == case_condition,"case","control"),
         condition = timepoint)

merged_seurat$cluster <- cleanSymbols(merged_seurat$cluster)
merged_seurat <- DietSeurat(merged_seurat,counts=TRUE,data=TRUE,scale.data=FALSE)
saveRDS(merged_seurat,file.path(pre_processing_path,"seurat_object_oi.rds"))
writeH5ADObjects(merged_seurat,pre_processing_path)
rm(merged_seurat)

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

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- expm1(data_matrix)
write.table(100 * data_matrix, file.path(natmi_data_filepath,"case/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"case/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"control"),recursive = TRUE)
seurat_object_oi_subset <- subset(seurat_object_oi,subset = condition == control_condition)

data_matrix <- GetAssayData(object = seurat_object_oi_subset, assay = "RNA", slot = "data")
data_matrix <- expm1(data_matrix)
write.table(100 * data_matrix, file.path(natmi_data_filepath,"control/em.txt"), quote = F, sep = "\t",row.names=TRUE,col.names=colnames(data_matrix))
meta_data <- seurat_object_oi_subset@meta.data %>%
  rownames_to_column(var="barcode") %>%
  rename(annotation=cluster)%>%
  select(barcode,annotation)
write.table(meta_data,file.path(natmi_data_filepath,"control/metadata.txt"), quote = F,sep="\t",row.names=FALSE,col.names=TRUE)

dir.create(file.path(natmi_data_filepath,"diff"))
print("pre-processing finalized")




renv::restore()

# GEO accession number GSE171964
# Load necessary libraries
library(dplyr)
library(Seurat)
library(Matrix)

set.seed(123)

dataset_path <- "data/covid"
results_path <- "results/covid"
pre_processing_path <- file.path(results_path,"pre_processing")
dir.create(results_path, recursive = TRUE)
dir.create(pre_processing_path, recursive = TRUE)


sparse_matrix <- Matrix::readMM(file.path(dataset_path,"GSE171964_unzipped/matrix.mtx"))
barcodes <- read.table(file.path(dataset_path,"GSE171964_unzipped/barcodes.tsv"))
features <- read.table(file.path(dataset_path,"GSE171964_unzipped/features.tsv"))

meta_data <- read.csv(file.path(dataset_path,"GSE171964_geo_pheno_v2.csv"))

rownames(sparse_matrix) <- features$V2
colnames(sparse_matrix) <- barcodes$V2

barcodes_to_retain <- meta_data %>% filter(day %in% c(0,22)) %>% select(barcode) %>% unlist()
ind_barcodes_to_retain <- match(barcodes_to_retain,colnames(sparse_matrix))
sparse_matrix_subset <- sparse_matrix[,ind_barcodes_to_retain]
rm(sparse_matrix)

seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix_subset,meta.data = meta_data)


saveRDS(seurat_object, file.path(pre_processing_path, "seurat_object.rds"))

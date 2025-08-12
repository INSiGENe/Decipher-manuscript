#main Decipher analysis pipeline
#V2
#load decipher package -----
library(devtools)
load_all()

#global options ----
set.seed(123)

#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
species <-  "human"
#for sample dataset condition_name is "condition", case_condition is "stim" and control_condition is "ctrl"
condition_name <- "condition"
case_condition = "stim"
control_condition = "ctrl"

#Parameters: directories ----
dataset_path <- "sample_analysis"
dir.create(dataset_path)
pre_processing_path <- file.path(dataset_path,"pre_processing")
reference_filepath <- "reference_data"
output_filepath <- dataset_path
output_data_filepath <- file.path(output_filepath,"data")
output_figures_filepath <- file.path(output_filepath,"figures")
output_importances_filepath <- file.path(output_filepath,"importances")
#directory set up----
dir.create(pre_processing_path)
dir.create(file.path(pre_processing_path,"h5ad_by_cluster"))
dir.create(output_data_filepath,recursive=TRUE)
dir.create(output_figures_filepath,recursive=TRUE)
dir.create(output_importances_filepath,recursive=TRUE)

#Parameters: analysis ----
flag.normalize.non.log <- FALSE

#create sample dataset ----
#including seurat object and h5ad objects
seurat_oi <- generateSampleSeuratFromExperimentHub(min_cells_per_cluster_condition,case_condition,control_condition)
##save outputs for Decipher analysis
saveRDS(seurat_oi,file.path("sample_analysis/pre_processing","seurat_object_oi.rds"))

seurat_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

#in addition, we need to create python-compatible h5ad objects for the CO pipeline, here, I've opted against it
#as they are not necessary for this script
for(this_cluster in unique(kang.seurat$cluster)){
   seurat_object_oi_this_cluster <- subset(kang.seurat,subset = cluster == this_cluster)
   sce.object = as.SingleCellExperiment(seurat_object_oi_this_cluster)
   sce.object@assays@data[["logcounts"]] <- NULL
   writeH5AD(sce.object, file.path("pre_processing/h5ad_by_cluster",paste(this_cluster,".h5ad",sep="")),X_name = "counts")
 }

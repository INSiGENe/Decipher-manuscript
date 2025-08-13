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

#set a seed
set.seed(123)
dataset_path <- "sample_analysis"
pre_processing_path <- file.path(dataset_path,"pre_processing")

kang.seurat <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))

#in addition, we need to create python-compatible h5ad objects for the CO pipeline, here, I've opted against it
#as they are not necessary for this script
for(this_cluster in unique(kang.seurat$cluster)){
   seurat_object_oi_this_cluster <- subset(kang.seurat,subset = cluster == this_cluster)
   sce.object = as.SingleCellExperiment(seurat_object_oi_this_cluster)
   sce.object@assays@data[["logcounts"]] <- NULL
   writeH5AD(sce.object, file.path("sample_analysis/pre_processing/h5ad_by_cluster",paste(this_cluster,".h5ad",sep="")),X_name = "counts")
 }

renv::restore()
library(devtools)
load_all()

#set up zellkonverter
zellkonverter_env <- zellkonverter:::zellkonverterAnnDataEnv("0.10.9")
basilisk::basiliskStart(zellkonverter_env)

#global options ----
i = 123
set.seed(i)
selected_random_seed <- i

#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
#for sample dataset condition_name is "condition", case_condition is "stim" and control_condition is "ctrl"
condition_name <- "condition"
case_condition = "stim"
control_condition = "ctrl"


#Parameters: directories ----
reference_filepath <- "reference_data"
dataset_path <- "results/sample_analysis"
paths <- create_project_dirs("results/sample_analysis")
#not sure what this is about: dir.create(file.path(paths['pre_processing'],"validity/h5ad_by_cluster"))

# #create sample dataset ----
# #including seurat object and h5ad objects
#will ask you to set up a cache, say yes
seurat_oi <- generateSampleSeuratFromExperimentHub(min_cells_per_cluster_condition,case_condition,control_condition)

# sampling so we can quickly go through the pipeline 
cells_to_keep <- seurat_oi[[]] %>%
  tibble::rownames_to_column("cell_barcode") %>%
  dplyr::group_by(condition, cluster) %>%
  dplyr::filter(dplyr::n() >= 300) %>%
  dplyr::slice_sample(n = 500, replace = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::pull(cell_barcode)

seurat_oi <- subset(seurat_oi, cells = cells_to_keep)

# #save outputs for subsequent analysis
saveRDS(seurat_oi,file.path(paths["data"],"seurat_object_oi.rds"))

#in addition, we need to create python-compatible h5ad objects for the CO pipeline
for(this_cluster in unique(seurat_oi$cluster)){
  seurat_object_oi_this_cluster <- subset(seurat_oi,subset = cluster == this_cluster)
  sce.object = Seurat::as.SingleCellExperiment(seurat_object_oi_this_cluster)
  sce.object@assays@data[["logcounts"]] <- NULL
  zellkonverter::writeH5AD(sce.object, file.path(paths["co_input"],paste(this_cluster,".h5ad",sep="")),X_name = "counts")
}

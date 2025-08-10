##########################
#### libraries
##########################
library(Matrix)
library(Seurat)
library(Azimuth)

##########################
#### analysis
##########################
set.seed(123)

query <- readRDS("data/SevMilCOVID/combined_seurat_for_processing.rds")

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
query <- SCTransform(query, verbose = FALSE)
DefaultAssay(query) <- "SCT"

# Load a reference (e.g., PBMC 10k reference; adjust if using another)
reference <- readRDS("reference_data/pbmc_multimodal.rds")
common.features <- intersect(rownames(reference), rownames(query))

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",       # depends on how your data was normalized
  reference.reduction = "refDR",
  features = common.features
)

seurat_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2"),
  reference.reduction = "refDR",
  reduction.model = "refUMAP"  # this one is used for plotting
)

seurat_mapped$severity_group <- ifelse(
  is.na(seurat_mapped$disease_severity) | seurat_mapped$disease_severity == "",
  "Healthy",
  seurat_mapped$disease_severity
)


# Save the new object
saveRDS(seurat_mapped, "data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")

#### end of analysis










query <- readRDS("data/SevMilCOVID/combined_seurat_for_processing.rds")


remotes::install_github("satijalab/azimuth@release/0.4.6", upgrade = "never")
BiocManager::install("glmGamPoi",update = FALSE)
remotes::install_github("satijalab/seurat-data@v0.2.1", upgrade = "never")
BiocManager::install("glmGamPoi",update = FALSE)
remotes::install_github("satijalab/azimuth@release/0.4.6", upgrade = "never")

query <- readRDS("data/SevMilCOVID/combined_seurat_for_processing.rds")
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
query <- SCTransform(query, verbose = FALSE)
DefaultAssay(query) <- "SCT"

# Load a reference (e.g., PBMC 10k reference; adjust if using another)
reference <- readRDS("reference_data/pbmc_multimodal.rds")
common.features <- intersect(rownames(reference), rownames(query))

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",       # depends on how your data was normalized
  reference.reduction = "refDR",
  features = common.features
)

seurat_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2"),
  reference.reduction = "refDR",
  reduction.model = "refUMAP"  # this one is used for plotting
)


docker run -it -v "$(pwd):/app" -w /app satijalab/azimuth R
docker run -it -v "$(pwd):/app" -w /app azimuth-custom R

install.packages("devtools")
install.packages("remotes")

for (p in c("spatstat.core","spatstat.geom","spatstat.random","spatstat.explore",
            "spatstat.sparse","spatstat.univar","spatstat.data","spatstat")){
  if (p %in% rownames(installed.packages())) remove.packages(p)
}

# 2) Install the meta-package so you get a consistent, compatible set
install.packages("spatstat", type = "source")


remotes::install_github("satijalab/azimuth@release/0.4.6", upgrade = "never")

wget https://zenodo.org/records/4546839/files/ref.Rds?download=1 -O reference_data/pbmc_multimodal.rds
wget https://zenodo.org/records/4546839/files/idx.annoy?download=1 -O reference_data/idx.annoy

install.packages("RcppAnnoy")
library(RcppAnnoy)

# Replace 50 with the number of dimensions used to build the index
# Suppose your refDR has 64 dimensions and used Euclidean
annoy_loaded <- new(AnnoyAngular, 50)
annoy_loaded$load("reference_data/idx.annoy")

annoy_neighbors <- annoy_loaded$getNNsByItemList(0, 31, search_k = -1, include_distances = TRUE)

# This returns a named list with:
# $item      -> integer vector of neighbor indices
# $distance  -> numeric vector of distances

# Access values like this:
annoy_neighbors$item      # neighbor indices (0-based)
annoy_neighbors$distance  # distances to those neighbors





########################
######## functions #####
########################

fix_refdr_colnames <- function(seurat_obj, reduction_name = "refDR", key_prefix = "refdr_") {
  if (!reduction_name %in% Reductions(seurat_obj)) {
    stop(sprintf("Reduction '%s' not found in Seurat object", reduction_name))
  }
  
  # Access the slot directly without triggering validation
  dr <- slot(seurat_obj[[reduction_name]], "feature.loadings")
  
  # Rename columns
  colnames(dr) <- paste0(key_prefix, seq_len(ncol(dr)))
  
  # Assign back into the DimReduc
  slot(seurat_obj[[reduction_name]], "feature.loadings") <- dr
  
  # Set the key
  Key(seurat_obj[[reduction_name]]) <- key_prefix
  
  return(seurat_obj)
}

#there's an issue with Matrix1.7.0 so please re-install 1.6.1.1 first
#you might have to pick a mirror
remotes::install_version("Matrix","1.6.1", type="source")
#now run what needs to be run
packageVersion("Matrix")

#install these things
remotes::install_github("satijalab/seurat-data@v0.2.1", upgrade = "never")
BiocManager::install("glmGamPoi",update = FALSE)
remotes::install_github("satijalab/azimuth@release/0.4.6", upgrade = "never")

library(Matrix)
library(Seurat)
library(Azimuth)

query <- readRDS("data/SevMilCOVID/combined_seurat_for_processing.rds")
query <- CreateSeuratObject(counts = object@assays$RNA@counts)



DefaultAssay(query) <- "RNA"
#v1
query_mapped <- RunAzimuth(query, reference = "pbmcref",verbose = TRUE)

DefaultAssay(query) <- "RNA"
query_mapped <- RunAzimuth(query, reference = "pbmcref",assay = "RNA",umap.name = "ref.umap",verbose = TRUE)

#v2
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
query <- SCTransform(query, verbose = FALSE)
DefaultAssay(query) <- "SCT"
query_mapped <- RunAzimuth(query, reference = "pbmcref",assay = "SCT",umap.name = "ref.umap",verbose = TRUE)


options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
query <- SCTransform(query, verbose = FALSE)
DefaultAssay(query) <- "SCT"

# Load a reference (e.g., PBMC 10k reference; adjust if using another)
reference <- readRDS("reference_data/pbmc_multimodal.rds")
reference <- fix_refdr_colnames(reference)

if ("map" %in% names(reference@tools)) {
  reference@tools$map <- NULL
}

reference@tools <- list()   

common.features <- intersect(rownames(reference), rownames(query))

# Step 2: Find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",       # depends on how your data was normalized
  reference.reduction = "refDR",
  features = common.features
)

query <- SCTransform(
  object = query,
  assay = if ("SCT" %in% Assays(query)) "SCT" else "RNA",
  new.assay.name = "refAssay",
  reference.SCT.model = reference[["refAssay"]]@SCTModel.list$refmodel,
  method = "glmGamPoi",
  ncells = 2000, n_genes = 2000,
  do.correct.umi = FALSE, do.scale = FALSE, do.center = TRUE,
  verbose = TRUE
)

# Step 3: Map query
seurat_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2"),
  reference.reduction = "refDR",
  reduction.model = "refUMAP"  # this one is used for plotting
)

seurat_mapped$severity_group <- ifelse(
  is.na(seurat_mapped$disease_severity) | seurat_mapped$disease_severity == "",
  "Healthy",
  seurat_mapped$disease_severity
)


# Save the new object
saveRDS(seurat_mapped, "data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")


test <- RunAzimuth(query, reference = "pbmcref",assay = "RNA")


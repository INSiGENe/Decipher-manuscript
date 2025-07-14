docker run -it -v "$(pwd):/app" -w /app satijalab/azimuth R
docker run -it -v "$(pwd):/app" -w /app azimuth-custom R

install.packages("devtools")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data@v0.2.1", upgrade = "never")
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



library(Seurat)
library(Azimuth)
query <- readRDS("data/SevMilCOVID/combined_seurat_for_processing.rds")

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
query <- SCTransform(query, verbose = FALSE)
DefaultAssay(query) <- "SCT"

# Load a reference (e.g., PBMC 10k reference; adjust if using another)
reference <- readRDS("reference_data/pbmc_multimodal.rds")

common.features <- intersect(rownames(reference), rownames(query))

# Step 2: Find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",       # depends on how your data was normalized
  reference.reduction = "refDR",
  features = common.features
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

# Check predicted labels
table(seurat_mapped$predicted.celltype.l1)
table(seurat_mapped$predicted.celltype.l2)
table(seurat_mapped$severity_group,seurat_mapped$predicted.celltype.l1)


# Save the new object
saveRDS(seurat_mapped, "data/SevMilCOVID/combined_seurat_for_processing_azimuth_mapped.rds")
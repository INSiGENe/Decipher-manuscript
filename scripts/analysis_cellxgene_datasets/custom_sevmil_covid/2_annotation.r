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

library(Seurat)
library(Matrix)

# Set base path
base_path <- "data/SevMilCOVID"

# Get list of cov folders (e.g., cov01, cov02, ...)
sample_dirs <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[grepl("^cov\\d+", sample_dirs)]  # Just covXX folders

# Load features file once
features_file <- file.path(base_path, "GSE155673_features.tsv")
features <- read.delim(features_file, header = FALSE)
gene_names <- make.unique(as.character(features$V2))

# List to store Seurat objects
seurat_list <- list()

# Loop through each folder
for (sample_id in sample_dirs) {
  cat("Processing", sample_id, "\n")
  
  # File paths
  mtx_file <- file.path(base_path, sample_id, paste0("GSE155673_", sample_id, "_matrix.mtx"))
  barcodes_file <- file.path(base_path, sample_id, paste0("GSE155673_", sample_id, "_barcodes.tsv"))
  
  # Read data
  expression_matrix <- readMM(mtx_file)
  barcodes <- readLines(barcodes_file)
  rownames(expression_matrix) <- gene_names
  colnames(expression_matrix) <- barcodes
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = expression_matrix, project = sample_id)
  seu$sample_id <- toupper(sample_id)  # Add sample_id like "COV01"
  
  seurat_list[[sample_id]] <- seu
}

# Merge all Seurat objects
combined_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = toupper(names(seurat_list))
)

# View result
combined_seurat


#meta data
sample_metadata <- data.frame(
  sample_id = c("COV01", "COV02", "COV03", "COV04", "COV07", "COV08", "COV09", "COV10", "COV11", "COV12", "COV17", "COV18"),
  gsm = c("GSM4712885", "GSM4712887", "GSM4712889", "GSM4712891", "GSM4712893", "GSM4712895", "GSM4712897",
          "GSM4712899", "GSM4712901", "GSM4712903", "GSM4712905", "GSM4712907"),
  sample_name = c("nCOV3EUHM", "nCOV7EUHM", "nCOV1EUHM", "nCOV6EUHM", "280", "259", "279", 
                  "nCOV021EUHM", "nCOV024EUHM", "nCOV0029EUHM", "265", "258"),
  age = c(75, 53, 75, 59, 84, 68, 38, 60, 48, 47, 90, 70),
  sex = c("F", "F", "F", "M", "F", "F", "M", "F", "M", "F", "M", "F"),
  disease_status = c("COVID-19", "COVID-19", "COVID-19", "COVID-19", "Healthy", "Healthy", "Healthy",
                     "COVID-19", "COVID-19", "COVID-19", "Healthy", "Healthy"),
  disease_severity = c("Severe", "Moderate", "Moderate", "Severe", NA, NA, NA, "Severe", "Severe", "Moderate", NA, NA),
  days_since_symptom_onset = c(15, 9, 2, 16, NA, NA, NA, 15, 8, 9, NA, NA),
  stringsAsFactors = FALSE
)

# View the table
print(sample_metadata)

#join  sample_metadata table to the combined_seurat@meta.data using the sample_id column.

# Store original cell names
meta <- combined_seurat@meta.data
meta$cell <- rownames(meta)

# Merge
merged_meta <- merge(
  meta,
  sample_metadata,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore cell names
rownames(merged_meta) <- merged_meta$cell
merged_meta$cell <- NULL

# Reassign to Seurat object
combined_seurat@meta.data <- merged_meta


# Confirm it's there
head(combined_seurat@meta.data)

###############################
#replicate their pre-processing
############################### 
set.seed(123)
library(Seurat)
library(dplyr)

DefaultAssay(combined_seurat) <- "RNA"

# Filter cells with >25% mitochondrial reads
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
combined_seurat <- subset(combined_seurat, subset = percent.mt <= 25)

# Normalize RNA counts
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify top 2000 variable genes
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)

# Scale data
combined_seurat <- ScaleData(combined_seurat)

# Run PCA
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))

# Use first 25 PCs
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:25)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.4)  # Louvain clustering
combined_seurat <- RunUMAP(combined_seurat, dims = 1:25)

# Check how many clusters you get
length(unique(Idents(combined_seurat)))

# Example: remove clusters with high mito and low features
cluster_qc <- combined_seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(
    mean_mito = mean(percent.mt),
    mean_features = mean(nFeature_RNA),
    n = n()
  )
print(cluster_qc,n=25)

# Suppose cluster "X" meets dead-cell criteria:
combined_seurat <- subset(combined_seurat, idents = "9",invert=TRUE)
combined_seurat <- subset(combined_seurat, idents = "4",invert=TRUE)
combined_seurat <- subset(combined_seurat, idents = "11",invert=TRUE)
combined_seurat <- subset(combined_seurat, idents = "15",invert=TRUE)



ref_deg <- read.csv("reference_data/SevMilCOVID19_cluster_labels.csv", stringsAsFactors = FALSE)

combined_seurat$severity_group <- ifelse(
  is.na(combined_seurat$disease_severity) | combined_seurat$disease_severity == "",
  "Healthy",
  combined_seurat$disease_severity
)

Idents(combined_seurat) <- combined_seurat$disease_status

covid_vs_healthy <- FindMarkers(
  combined_seurat,
  ident.1 = "COVID-19",
  ident.2 = "Healthy",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

library(dplyr)

# Save cluster IDs for looping
clusters <- levels(combined_seurat$seurat_clusters)

# Initialize list
cluster_markers <- list()

Idents(combined_seurat) <- combined_seurat$seurat_clusters

for (clust in clusters) {
  # Subset to cells from this cluster
  cells_in_cluster <- WhichCells(combined_seurat, idents = clust)
  
  # Subset Seurat object and set identities to severity
  sub_obj <- subset(combined_seurat, cells = cells_in_cluster)
  Idents(sub_obj) <- sub_obj$disease_status
  
  # Run comparison: COVID vs Healthy
  cluster_markers[[clust]] <- FindMarkers(
    sub_obj,
    ident.1 = "COVID-19",
    ident.2 = "Healthy",
    min.pct = 0.25,
    logfc.threshold = 0.25
  ) %>% mutate(cluster = clust, gene = rownames(.))
}

# Combine all results into one data frame
cluster_markers_df <- bind_rows(cluster_markers)


library(dplyr)

# Reduce reference to gene + cluster + logFC
ref_summary <- ref_deg %>%
  filter(!is.na(gene)) %>%
  select(ref_cluster = clust, gene, ref_logFC = avg_logFC)

# Reduce your markers to gene + cluster + logFC
my_summary <- cluster_markers_df %>%
  select(my_cluster = cluster, gene, my_logFC = avg_log2FC)

# Inner join on gene
merged <- inner_join(my_summary, ref_summary, by = "gene")

# Correlate logFCs between your cluster and each ref cluster
mapping_results <- merged %>%
  group_by(my_cluster, ref_cluster) %>%
  summarize(correlation = cor(my_logFC, ref_logFC, use = "complete.obs")) %>%
  ungroup()

# Get the best match for each cluster
best_matches <- mapping_results %>%
  group_by(my_cluster) %>%
  slice_max(order_by = correlation, n = 1) %>%
  ungroup()

print(best_matches)

ref_cluster_labels <- ref_deg %>%
  select(ref_cluster = clust, ref_cluster_name = final) %>%
  distinct()

best_matches <- best_matches %>%
  left_join(ref_cluster_labels, by = "ref_cluster")

cluster_map <- setNames(best_matches$ref_cluster_name, best_matches$my_cluster)


library(dplyr)

meta <- combined_seurat@meta.data %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  left_join(
    tibble(
      seurat_clusters = names(cluster_map),
      mapped_cluster = unname(cluster_map)
    ),
    by = "seurat_clusters"
  )

# Keep rownames
rownames(meta) <- rownames(combined_seurat@meta.data)
combined_seurat@meta.data <- meta

# Explicitly get cells to keep
cells_to_keep <- rownames(combined_seurat@meta.data)[!is.na(combined_seurat$mapped_cluster)]

# Subset using cells
combined_seurat <- subset(combined_seurat, cells = cells_to_keep)

Idents(combined_seurat) <- combined_seurat$mapped_cluster

saveRDS(combined_seurat, file = file.path(base_path, "combined_seurat_for_processing.rds"))

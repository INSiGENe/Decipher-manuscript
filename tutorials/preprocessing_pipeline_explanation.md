# Explanation of the Pre-processing Script (`2_preprocess_object_for_analysis.R`)

This document breaks down the R script responsible for preparing the initial dataset for all subsequent analyses, including Decipher, CytoSig, LIANA+, and NATMI. This script acts as a bridge between the raw, downloaded data and the structured inputs required by each analysis tool.

## 1. Setup and Initialization

The script begins by setting up its environment.

- **Libraries:** It loads all necessary R packages.
- **Python Environment:** It initializes a specific Python environment using `basilisk` and `zellkonverter` to handle `.h5ad` files.

```R
#libraries ----
library(Seurat)
library(dplyr)
# ... and other libraries

#python env setup ----
zellkonverter_env <- zellkonverter:::zellkonverterAnnDataEnv("0.10.9")
basilisk::basiliskStart(zellkonverter_env)
```

## 2. Helper Functions

The script defines several utility functions to keep the code organized:

- `writeH5ADObjects`: Splits the main Seurat object by cluster and saves each as a separate `.h5ad` file for CellOracle.
- `cleanSymbols`: Removes special characters from metadata strings.
- `get_column`: Flexibly finds a column in a dataframe from a list of possible names.

```R
#functions ----
writeH5ADObjects <- function(seurat_object, pre_processing_path) {
  # ... function implementation
}

cleanSymbols <- function(string) {
  # ... function implementation
}

get_column <- function(df, possible_names) {
  # ... function implementation
}
```

## 3. Parameter Loading and Directory Setup

The script is driven by the `config.json` file.

- It takes a `dataset_key` as a command-line argument.
- It reads the `pre_processing_for_analysis` section from the config to get paths and metadata fields.
- It creates the output directory structure for all analysis tools.

```R
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
dataset_key <- args[1]
config <- jsonlite::fromJSON(txt = "scripts/config.json")
cfg <- config[[dataset_key]][["pre_processing_for_analysis"]]

#directories
input_path <- cfg$input_path
results_path <- cfg$output_path
# ... and other directory creations
```

## 4. Core Pre-processing for Decipher

This is the first major processing step. It checks if a processed Seurat object (`seurat_object_oi.rds`) already exists. If not, it performs the following actions:

1.  **Data Loading:** Reads the raw matrix, gene names, and metadata.
2.  **Gene ID Conversion:** Converts Ensembl IDs to HGNC symbols, caching the mapping for future runs.
3.  **Metadata Cleaning:** Cleans cluster names and creates standardized `condition` and `sample_key` columns.
4.  **Seurat Object Creation & Saving:** Assembles and saves the final Seurat object.

```R
#Decipher pre-processing ----
seurat_object_rds_path <- file.path(pre_processing_path, "seurat_object_oi.rds")

if (!file.exists(seurat_object_rds_path)) {
  matrix <- Matrix::readMM(file.path(input_path, "output_matrix.mtx"))
  # ... (transpose, read metadata)

  # Use dynamic column selection
  rownames(matrix) <- get_column(gene_names, c("X0", "name", "gene", "gene_name","index","gene_ids","ensembl_id"))
  # ... (set colnames, rownames)

  # Convert Ensembl IDs to HGNC gene symbols
  # ... (logic to generate or load mapping)
  matrix <- matrix[mapping$ENSEMBL, ]
  rownames(matrix) <- mapping$SYMBOL

  meta_data$cluster <- meta_data[[cfg$cluster_meta_field]]
  meta_data$cluster <- cleanSymbols(meta_data$cluster)
  meta_data$condition <- meta_data[[cfg$condition_meta_field]]
  # ... (create sample_key)

  seurat_object_oi <- Seurat::CreateSeuratObject(counts = matrix,meta.data = meta_data)

  saveRDS(seurat_object_oi,seurat_object_rds_path)
} else {
  message("Seurat object already exists...")
}

seurat_object_oi <- readRDS(seurat_object_rds_path)
writeH5ADObjects(seurat_object_oi,pre_processing_path)
```

## 5. Tool-Specific Pre-processing

The script then creates specialized input files for other benchmarking tools.

### 5.1. CytoSig Pre-processing

CytoSig requires a specific differential expression profile. The script iterates through each cluster, calculates the profile, and saves it.

```R
# CYTOSIG pre-processing ----
for(this_cluster in unique(seurat_object_oi$cluster)){
  # ... (check if output exists)

  seurat_object_oi_this_cluster <- subset(seurat_object_oi,idents = this_cluster)
  # ... (filter genes, normalize, log-transform)

  # Extract data for control and case conditions
  control_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == control_condition)
  case_cells <- which(seurat_object_oi_this_cluster@meta.data$condition == case_condition)

  # ... (skip if not enough cells)

  # Calculate the differential expression profile
  differential_profile <- case_data-control_mean_expression

  # Write the profile to a compressed file
  gz1 <- gzfile(file.path(cytosig_cluster_path,"differential_profile.tsv.gz"), "w")
  write.table(differential_profile, gz1, sep = "	", col.names = NA, quote = FALSE)
  close(gz1)
}
```

### 5.2. LIANA+ Pre-processing

Converts the main Seurat object into a single `.h5ad` file for LIANA+.

```R
# LIANA pre-processing ----
liana_h5ad_path <- file.path(liana_data_filepath, "seurat_object_oi.h5ad")
if (!file.exists(liana_h5ad_path)) {
  seurat_object_oi <- readRDS(seurat_object_rds_path)
  # ... (add liana_condition column)
  sce.object = as.SingleCellExperiment(seurat_object_oi)
  writeH5AD(sce.object, liana_h5ad_path,X_name = "counts")
}
```

### 5.3. NATMI Pre-processing

NATMI requires separate expression and metadata files for the "case" and "control" conditions.

```R
# NATMI pre-processing ----
seurat_object_oi = readRDS(seurat_object_rds_path)
# ... (normalize data)

# For case condition:
if (!file.exists(natmi_case_em_path) || !file.exists(natmi_case_metadata_path)) {
  seurat_object_oi_subset <- seurat_object_oi[, which(seurat_object_oi$condition == case_condition),seed=NULL]
  # ... (get and write expression matrix)
  # ... (get and write metadata)
}

# For control condition:
if (!file.exists(natmi_control_em_path) || !file.exists(natmi_control_metadata_path)) {
  seurat_object_oi_subset <- seurat_object_oi[, which(seurat_object_oi$condition == control_condition),seed=NULL]
  # ... (get and write expression matrix)
  # ... (get and write metadata)
}
```

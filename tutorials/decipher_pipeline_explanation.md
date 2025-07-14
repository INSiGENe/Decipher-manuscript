# Explanation of the Decipher Analysis Pipeline (`6_decipher_pipeline_v1_modularized.R`)

This document provides a detailed breakdown of the main Decipher analysis script. The script is designed to take a pre-processed Seurat object, run the core Decipher workflow, and save the results.

## 1. Initialization and Setup

The script begins by loading the necessary `Decipher` package and setting a global random seed for reproducibility.

```R
#load decipher package -----
library(devtools)
load_all()

#global options ----
#Set this seed to NULL if you don't need reproducible results
selected_random_seed = 123
set.seed(selected_random_seed)
```

### 1.1. Parameter Loading

The script is configured via command-line arguments and a central `config.json` file.

- **Command-Line Argument:** The script expects a single argument, which is the `dataset_key` (e.g., `cz_rcc`, `cz_influenza`). This key is used to retrieve the correct configuration settings.
- **Configuration File:** It reads `scripts/config.json` to fetch the specific parameters for the given `dataset_key` under the `Decipher_analysis` section.

```R
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a dataset key, e.g. 'cz_rcc'")
}
dataset_key <- args[1]

config <- jsonlite::fromJSON(txt = "scripts/config.json")

if (!dataset_key %in% names(config)) {
  stop(paste("Dataset key not found in config:", dataset_key))
}

dcp <- config[[dataset_key]][["Decipher_analysis"]]

min_cells_per_cluster_condition <- dcp$min_cells_per_cluster_condition
species <- dcp$species
dataset_path <- dcp$dataset_path
condition_name <- dcp$condition_name
case_condition <- dcp$case_condition
control_condition <- dcp$control_condition
k_parameter <- dcp$k_parameter
min_meta_cells_parameter = dcp$min_meta_cells_parameter
```

### 1.2. Directory Setup

The script creates a standardized directory structure within the specified `dataset_path` to store all outputs.

```R
#Parameters: directories ----
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
```

## 2. Data Loading

### 2.1. Seurat Object

The main input is a Seurat object named `seurat_object_oi.rds`.

```R
#load dataset ----
seurat_oi <- readRDS(file.path(pre_processing_path,"seurat_object_oi.rds"))
```

### 2.2. Reference Data

The script loads essential reference databases.

```R
#load reference data ----
L.set <- loadLSet(reference_filepath,species)
enrichr_database <- loadEnrichrDatabase(reference_filepath,species)
cytosig_ligands <- loadCytosigLigands(reference_filepath,species)
```

## 3. Pre-processing and Quality Control (QC)

- **Condition Mapping:** Standardizes the condition information in the Seurat object.
- **Cell Count QC:** Filters out clusters that do not have enough cells in each condition.

```R
#data pre-processing ----
seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi,condition_name,case_condition,control_condition)

##QC ----
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
plotQC_CpC(CpC_data,outputPath=output_figures_filepath)

#PARAM: select the minimum number of cells per cluster + condition
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,minCpc = 100)
seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed=NULL]
plotQC_UpC(seuratObject = seurat_oi,outputPath = output_figures_filepath,id = "_scd")
```

## 4. Meta-Cell Creation

Aggregates single cells into "meta-cells" (or pseudobulks) to manage computational complexity and reduce noise.

```R
##Meta cells ----
decipher_seurat <- metaCellModule(
  seurat_object = seurat_oi,
  min_meta_cells = min_meta_cells_parameter,
  k = k_parameter
)

plotQC_UpC(seuratObject = decipher_seurat,outputPath = output_figures_filepath,id = "_meta")

CpC_data_meta <- generateQCDataByClusterAndCondition(decipher_seurat,max(stringr::str_length(unique(decipher_seurat$cluster))))
plotQC_CpC(CpC_data_meta,outputPath=output_figures_filepath,id = "_meta")

if(is.null(k_parameter)){
  k_parameter = calculate_suggested_number_of_metacell_neighbours(seurat_oi,min_meta_cells_parameter)
}

# ... (parameter saving)

saveRDS(decipher_seurat,file.path(output_data_filepath,"pseudobulk_seurat.rds"))
```

## 5. Core Decipher Analysis

This section contains the main calculations of the Decipher algorithm, performed on the meta-cell object.

### 5.1. Feature Filtering
Identifies which ligands and receptors from the `L.set` are sufficiently expressed in the data.

```R
#data pre-processing: main analysis ----
decipher_seurat_lr <- decipher_seurat[unique(c(L.set$ligand,L.set$receptor)),, seed=NULL]

feature_statistics <- getFeatureStatistics(
  features=unique(c(L.set$ligand,L.set$receptor)),
  seuratObj=decipher_seurat)

expressed_ligands <- getFilteredLigands(
  decipher_seurat,
  L.set,
  param_min_ligand_expr_in_cluster = 0.1)

expressed_receptors_all_clusters <- getExpressedReceptorsForEachCluster(
  decipher_seurat,
  L.set)

L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(
  L.set,
  expressed_ligands,
  expressed_receptors_all_clusters)
```

### 5.2. Regulon Analysis
Infers and analyzes the activity of gene regulatory networks (GRNs).

```R
regulon_grns_by_cluster <- getRegulonsAllClusters(
  output_filepath,
  decipher_seurat)

capped_regulons_all_clusters <- capRegulonsAllClusters(
  regulon_grns_by_cluster,
  decipher_seurat,
  flag.normalize.non.log)

regulon_scores_by_cluster <- getRegulonScoresAllClusters(
  capped_regulons_all_clusters,
  decipher_seurat)

regulon_deltas_by_cluster <- getRegulonDeltasAllClusters(
  regulon_scores_by_cluster,
  decipher_seurat)

significant_regulons_by_cluster <- getSignificantRegulonsAllClusters(
  regulon_deltas_by_cluster)

significant_regulon_markers_by_cluster <- getDifferentiallyExpressedTargetsForRegulonsAllClusters(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_grns_by_cluster,
  flag.normalize.non.log,
  random.seed=selected_random_seed)
```

### 5.3. Interaction Potential Calculation
Calculates the potential for interaction between cell clusters based on ligand-receptor expression.

```R
interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClusters(
  decipher_seurat,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log)

interaction_deltas_by_cluster <- calculateInteractionDeltasAllClusters(
  interaction_potential_by_clusters,
  decipher_seurat_lr)

filtered_interaction_potentials_matrix_all_clusters <- filterIntPotByDeltas(
  interaction_potential_by_clusters,
  interaction_deltas_by_cluster)

interaction_potentials_matrix_clusters_all_clusters <-
  getInteractionPotentialMatrixForRepresentativeInteractionsAllClusters(
    decipher_seurat,
    L_set_relevant_features_all_clusters,
    filtered_interaction_potentials_matrix_all_clusters,
    cytosig_ligands,
    flag.normalize.non.log)
```

### 5.4. Random Forest Modeling (The "Decipher Score")
This is the central step where a Random Forest model links intercellular interactions to intracellular responses. The feature importance of each interaction in predicting a regulon's activity is the **Decipher Score**.

```R
decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClusters(
  decipher_seurat,
  significant_regulons_by_cluster,
  regulon_scores_by_cluster,
  interaction_potentials_matrix_clusters_all_clusters,
  L_set_relevant_features_all_clusters,
  flag.normalize.non.log)
```

### 5.5. Differential Expression
Finds differentially expressed genes for additional context.

```R
lr_markers_by_cluster <- FindLRMarkersAllClusters(
  decipher_seurat,
  decipher_scores_by_regulon_and_cluster,
  flag.normalize.non.log,
  random.seed = selected_random_seed
)

de_markers_by_cluster <- FindMarkersAllClusters(
  decipher_seurat,
  flag.normalize.non.log,
  random.seed= selected_random_seed
)
```

## 6. Finalization and Output

- **Score Aggregation:** The raw Random Forest weights are processed into final scores.
- **Saving Results:** All intermediate and final results are saved to `.rds` files.

```R
#DECIPHER analysis-----
decipher_scores_by_regulon_and_cluster <- lapply(
  decipher_scores_by_regulon_and_cluster,
  FUN = "listOfDFsRenameColumn",
  original_name = "weighted.spearman.cont",
  new_name = "decipher_score")

decipher_scores_by_cluster <- lapply(
  decipher_scores_by_regulon_and_cluster,
  FUN = "calculateScoresByCluster")

decipher_scores_by_cluster <- addListNameToDFElements(
  decipher_scores_by_cluster,
  "receiver_cluster")

#save DECIPHER ----
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(regulon_scores_by_cluster,file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
# ... (many more saveRDS calls)
```

## 7. Visualization

Finally, the script generates a summary plot of the top prioritized interactions.

```R
#plot results ----
plotDecipherPrioritizedMap(dataset_path,top_n=6,dataset_name="dataset_key")
```

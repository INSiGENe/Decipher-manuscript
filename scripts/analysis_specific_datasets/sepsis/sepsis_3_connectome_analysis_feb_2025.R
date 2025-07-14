

############
#libraries ----
############
library(devtools)
library(Connectome)
library(Seurat)
library(dplyr)

#############
#functions ----
############
getForrestLRDatabase <- function(file) {
  if (!file.exists(file)) stop("The specified file does not exist.")
  L.set <- read.csv(file)[, c("Ligand.gene.symbol", "Receptor.gene.symbol")]
  names(L.set) <- c("ligand", "receptor")
  L.set$lr <- paste(L.set$ligand, L.set$receptor, sep = "-")
  return(L.set)
}

replace_clean_with_original <- function(df, seurat_obj) {
  # Create a mapping table of clean cluster names to original names
  cluster_mapping <- unique(data.frame(
    clean_name = seurat_obj$cluster_clean,
    original_name = seurat_obj$cluster
  ))

  # Ensure unique mapping (in case of duplicates)
  cluster_mapping <- cluster_mapping[!duplicated(cluster_mapping$clean_name), ]

  # Function to replace names based on the mapping
  replace_names <- function(column) {
    df[[column]] <- cluster_mapping$original_name[match(df[[column]], cluster_mapping$clean_name)]
    return(df[[column]])
  }

  # Replace in the 'source' and 'target' columns
  df$source <- replace_names("source")
  df$target <- replace_names("target")

  # Update the 'edge' column by reconstructing it with original names
  df$edge <- paste(df$source, df$ligand, df$receptor, df$target, sep = " - ")

  return(df)
}

############
#data and analysis ----
############

## Global Settings ----
set.seed(123)

## Dataset Parameters ----
this_species = "human"
case_condition = "ICU-SEP"
control_condition = "ICU-NoSEP"

## Paths ----
dataset_path <- "results/sepsis"
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- "reference_data"
output_filepath <- file.path(dataset_path,"connectome")
output_data_filepath <- file.path(output_filepath,"data")

## Directory Setup ----
dir.create(output_data_filepath,recursive=TRUE)

## Load Data ----
seuratObj = readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))
DefaultAssay(seuratObj) <- "RNA"

## Data Pre-processing ----
Idents(seuratObj) <- seuratObj$cluster
seuratObj$cluster_clean <- stringr::str_replace_all(seuratObj$cluster,"_"," ")

## Load Ligand-Receptor Database ----
L.set <- getForrestLRDatabase(file.path(reference_filepath, "connectomedb_forrest_lrc2p.csv"))
L.set.for.Connectome <- L.set %>% 
  select(ligand,receptor) %>%
  mutate(modal_cat = 'UNCAT')

## Split Seurat Object by Condition ----
seuratObj.list <- SplitObject(seuratObj,split.by = 'condition')

## Initialize List for Connectome Results ----
seuratObj.con.list <- list()

## Process Case Condition ----
seuratObj.list[[case_condition]] <- seuratObj.list[[case_condition]] %>%
 NormalizeData() %>%
 ScaleData(features = rownames(.))

Idents(seuratObj.list[[case_condition]]) <- seuratObj.list[[case_condition]]$cluster_clean

seuratObj.con.list[[case_condition]] <- Connectome::CreateConnectome(
  seuratObj.list[[case_condition]],
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)

## Process Control Condition ----
seuratObj.list[[control_condition]] <- seuratObj.list[[control_condition]] %>%
  NormalizeData() %>%
  ScaleData(features = rownames(.))

Idents(seuratObj.list[[control_condition]]) <- seuratObj.list[[control_condition]]$cluster_clean

seuratObj.con.list[[control_condition]] <- Connectome::CreateConnectome(
  seuratObj.list[[control_condition]],
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)

## Compute Differential Connectome ----
seuratObj.con.diff <- DifferentialConnectome(
  connect.ref = seuratObj.con.list[[control_condition]], 
  connect.test = seuratObj.con.list[[case_condition]])

## Restore Original Cluster Names ----
seuratObj.con.diff <- replace_clean_with_original(seuratObj.con.diff, seuratObj)

## Save Results ----
saveRDS(seuratObj.con.diff,file.path(output_data_filepath,"connectome_results.rds"))

## Completion Message ----
print("connectome analysis ended")

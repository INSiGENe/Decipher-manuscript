# Script Name: connectome analysis
# Author: Edgar
# Date: 6/12/2023
# Description: Script to run and save differential connectome analysis
library(devtools)
#load_all()

#libraries ----
library(Connectome)
library(Seurat)
#global settings ----
set.seed(123)

#dataset  parameters ----
this_species = "human"
case_condition = "PIC"
control_condition = "CTRL"

#parameters ----
dataset_path <- "results/5yr_pic"
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
reference_filepath <- "reference_data"
output_filepath <- file.path(dataset_path,"connectome")
output_data_filepath <- file.path(output_filepath,"data")

#directory set up----
dir.create(output_data_filepath,recursive=TRUE)

#load data ----
seuratObj = readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))
DefaultAssay(seuratObj) <- "RNA"
#data pre-processing ----
Idents(seuratObj) <- seuratObj$cluster

L.set <- getForrestLRDatabase(file.path(reference_filepath,"connectomedb_forrest_lrc2p.csv"))
L.set.for.Connectome <- L.set %>% select(ligand,receptor) %>%
  mutate(modal_cat = 'UNCAT')

# Split the object by condition:
seuratObj.list <- SplitObject(seuratObj,split.by = 'condition')
rm(seuratObj)
# Normalize, Scale, and create Connectome:
seuratObj.con.list <- list()

seuratObj.list[[case_condition]] <- NormalizeData(seuratObj.list[[case_condition]])

seuratObj.list[[case_condition]] <- ScaleData(seuratObj.list[[case_condition]],features = rownames(seuratObj.list[[case_condition]]))

seuratObj.con.list[[case_condition]] <- CreateConnectome(
  seuratObj.list[[case_condition]],
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)

seuratObj.list[[control_condition]] <- NormalizeData(seuratObj.list[[control_condition]])

seuratObj.list[[control_condition]] <- ScaleData(seuratObj.list[[control_condition]],features = rownames(seuratObj.list[[control_condition]]))

seuratObj.con.list[[control_condition]] <- CreateConnectome(
  seuratObj.list[[control_condition]],
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)


seuratObj.con.diff <- DifferentialConnectome(connect.ref = seuratObj.con.list[[control_condition]], connect.test = seuratObj.con.list[[case_condition]])
saveRDS(seuratObj.con.diff,file.path(output_data_filepath,"connectome_results.rds"))
print("connectome analysis ended")

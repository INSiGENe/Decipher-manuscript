# Script Name: connectome analysis
# Author: Edgar
# Date: 6/12/2023
# Description: Script to run and save differential connectome analysis
library(devtools)
load_all()

#libraries ----
library(tidyverse)
library(Connectome)

#global settings ----
set.seed(123)

#dataset  parameters ----
this_species = "human"
case_condition = "Healthy"
control_condition = "Managed"

#parameters ----
dataset_path <- "manuscript_analysis/lupus"
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

#case condition
seurat_obj_case <- subset(seuratObj,subset = condition == case_condition)
seurat_obj_case <- Seurat::NormalizeData(seurat_obj_case)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj_case)]
seurat_obj_case <- Seurat::ScaleData(seurat_obj_case,features = genes)

seurat_obj_case_connectome <- CreateConnectome(
  seurat_obj_case,
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)

rm(seurat_obj_case)

seurat_obj_control <- subset(seuratObj,subset = condition == control_condition)
seurat_obj_control <- Seurat::NormalizeData(seurat_obj_control)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj_control)]
seurat_obj_control <- Seurat::ScaleData(seurat_obj_control,features = genes)

seurat_obj_control_connectome <- CreateConnectome(
  seurat_obj_control,
  species = this_species,
  p.values = F,
  LR.database = "custom",
  custom.list = L.set.for.Connectome)

rm(seurat_obj_control)

seuratObj.con.diff <- DifferentialConnectome(connect.ref = seurat_obj_control_connectome, connect.test = seurat_obj_case_connectome)
saveRDS(seuratObj.con.diff,file.path(output_data_filepath,"connectome_results.rds"))
print("connectome analysis ended")

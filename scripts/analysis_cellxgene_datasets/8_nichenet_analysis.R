
library(nichenetr) # Please update to v2.0.4
library(tidyverse)
library(Seurat)

#global settings ----
set.seed(123)


############
#parameters ----
############

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a dataset key, e.g. 'cz_rcc'")
}
dataset_key <- args[1]

config <- jsonlite::fromJSON(txt = "scripts/config.json")

if (!dataset_key %in% names(config)) {
  stop(paste("Dataset key not found in config:", dataset_key))
}

nhnt <- config[[dataset_key]][["NicheNet_analysis"]]

if (is.null(nhnt)) {
  stop(paste("No 'pre_processing_for_analysis' section found for dataset:", dataset_key))
}

case_condition <- nhnt$case_condition
control_condition <- nhnt$control_condition
dataset_path <- nhnt$dataset_path
this_species <- nhnt$this_species

#parameters ----
pre_processing_filepath <- file.path(dataset_path,"pre_processing")
output_filepath <- file.path(dataset_path,"nichenet")
output_data_filepath <- file.path(output_filepath,"data")

#directory set up----
dir.create(output_data_filepath,recursive=TRUE)

#load data ----
seuratObj = readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))
DefaultAssay(seuratObj) <- "RNA"
Idents(seuratObj) <- paste(seuratObj$cluster,seuratObj$condition)
#focus on a subset of cells
cell.list <- WhichCells(seuratObj, downsample = 2000)
seuratObj <- subset(seuratObj,cells = cell.list)

#data pre-processing ----
Idents(seuratObj) <- seuratObj$cluster

# For newer Seurat versions, you may need to run the following
seuratObj <- UpdateSeuratObject(seuratObj)
if(this_species == "human"){
  lr_network = readRDS("reference_data/nichenet/lr_network_human_21122021.rds")
  ligand_target_matrix = readRDS("reference_data/nichenet/ligand_target_matrix_nsga2r_final.rds")
  weighted_networks = readRDS("reference_data/nichenet/weighted_networks_nsga2r_final.rds")
}

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

ligand_activities_all_clusters <- list()
prior_table_all_clusters <- list()
ct_cond_tab <- table(seuratObj$cluster, seuratObj$condition)
#for each receiver cell type
for(this_receiver_ct in unique(seuratObj$cluster)){
  receiver = this_receiver_ct
  if (min(ct_cond_tab[this_receiver_ct, ]) < 3) {
      message(sprintf("Skipping %s – fewer than 3 cells in at least one condition", 
                      this_receiver_ct))
      next
  }
 # get genes with non-zero expression in at least 5% of cells on the selected receiver cell type
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.05,assay_oi="RNA")
# treat as background those genes that occur in the ligand-target matrix
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  ## sender, all cell type
  sender_celltypes = unique(seuratObj$cluster)

  # get genes with non-zero expression in at least 5% of cells in each sender cell type
  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05,assay_oi="RNA") # lapply to get the expressed genes of every sender cell type separately here
  # unlist these sender-related genes
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

  # 2. Define a gene set of interest: these are the genes in the “receiver/target” cell population that are potentially affected by ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)

  seurat_obj_receiver= subset(seuratObj, idents = receiver)

  if(min(table(seurat_obj_receiver$condition))<3){
    next
  }
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["condition", drop=TRUE]])

# get differentially expressed genes in the receiver cell type (only those in the ligand-target matrix)
  DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = case_condition, ident.2 = control_condition, min.pct = 0.10,min.cells.group=3) %>% rownames_to_column("gene")
# p value lt 0.05 and log2FC gt 0.25
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

  if(length(geneset_oi) == 0){
    print(paste(this_receiver_ct,"does not have differentially expressed genes, skipping it"))
    next
  }

  # 3. Define a set of potential ligands
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)

  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  # 4. Perform NicheNet ligand activity analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

  ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
  ligand_activities_all_clusters[[receiver]] <- ligand_activities

  #prioritization of ligand-receptor pairs
  # By default, ligand_condition_specificty and receptor_condition_specificty are 0
  prioritizing_weights = c("de_ligand" = 1,
                           "de_receptor" = 1,
                           "activity_scaled" = 2,
                           "exprs_ligand" = 1,
                           "exprs_receptor" = 1,
                           "ligand_condition_specificity" = 0.5,
                           "receptor_condition_specificity" = 0.5)

  lr_network_renamed <- lr_network %>% rename(ligand=from, receptor=to)

  # Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
  DE_table <- calculate_de(seuratObj, celltype_colname = "cluster",
                           condition_colname = "condition", condition_oi = case_condition,
                           features = union(expressed_ligands, expressed_receptors))

  # Average expression information - only for LCMV condition
  expression_info <- get_exprs_avg(seuratObj, "cluster", condition_colname = "condition", condition_oi = case_condition)

  # Calculate condition specificity - only for datasets with two conditions!
  condition_markers <- FindMarkers(object = seuratObj, ident.1 = case_condition, ident.2 = control_condition,
                                   group.by = "condition", min.pct = 0, logfc.threshold = 0,
                                   features = union(expressed_ligands, expressed_receptors)) %>% rownames_to_column("gene")

  # Combine DE of senders and receivers -> used for prioritization
  processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                            senders_oi = sender_celltypes, receivers_oi = receiver)

  processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)

  processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_renamed)

  prior_table <- generate_prioritization_tables(processed_expr_table,
                                                processed_DE_table,
                                                ligand_activities,
                                                processed_condition_markers,
                                                prioritizing_weights = prioritizing_weights)
  prior_table %>%
    select(sender,receiver,ligand,receptor,prioritization_score)

  prior_table_all_clusters[[receiver]] <- prior_table
}

prior_table_all_clusters_bind_rows <- bind_rows(prior_table_all_clusters)


result <- prior_table_all_clusters_bind_rows %>%
  # Create interaction column
  mutate(interaction = paste(ligand, receptor, sep = "-")) %>%
  # Group by receiver
  group_by(receiver) %>%
  # Filter for top 5 interactions based on prioritization_score
  top_n(5, wt = prioritization_score) %>%
  # Select the desired columns
  select(sender, receiver, interaction, prioritization_score) %>%
  # Order by receiver and then by prioritization_score in descending order
  arrange(receiver, desc(prioritization_score)) %>%
  ungroup()

write.csv(result,file.path(output_data_filepath,"nichenet_lr_top_5_interactions_per_rct.csv"))
saveRDS(prior_table_all_clusters,file.path(output_data_filepath,"prior_table_all_clusters.rds"))
saveRDS(ligand_activities_all_clusters,file.path(output_data_filepath,"nichenet_results.rds"))
print("NicheNet analysis complete")

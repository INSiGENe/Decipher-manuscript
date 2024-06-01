#main Decipher analysis pipeline
#Warning: please ensure that the seurat object is pre-processed as per our vignette before running
#load decipher package -----
library(devtools)
load_all()
#libraries ----
library(dplyr) #1.1.2
library(ggplot2) #3.4.2
library(randomForest) #4.7-1.1
library(ggbeeswarm) #0.7.2
library(pagoda2) #1.0.10
library(Seurat) #4.3.0.1
library(enrichR) #3.2
library(BMS) #0.3.5
library(babelgene) #22.9
#library(magrittr) #2.0.3
#library(tibble) #3.2.1
#library(stringr) #1.5.0
#library(Matrix) #1.6-0
#library(SeuratObject) #4.1.3
#library(tidyr) #1.3.0

#the libraries below are to create a sample analysis
library(tidyr) #1.3.0
library(stringr) #1.5.0
library(basilisk) #1.11.2
library(zellkonverter) #1.9.0

#global options ----
set.seed(123)

#user defined functions ----
#to add to package
#meta cell module
metaCellModule <- function(seurat_object,min_meta_cells){
  suggested_number_of_metacell_neighbours <- calculate_suggested_number_of_metacell_neighbours(seurat_object,min_meta_cells)

  MetaCellMatrices <- generateMetaCellMatrices(
    seuratObj = seurat_object,
    paramMaxScCells = 1200*(suggested_number_of_metacell_neighbours+1),
    paramK = suggested_number_of_metacell_neighbours,
    paramMinMetaCells = min_meta_cells)

  parameter_record <- data.frame(
    "k" = k_parameter,
    "min_meta_cells" = min_meta_cells
  )
  write.csv(parameter_record,file.path(output_data_filepath,"parameter_record.csv"))

  seurat_pseudo_bulk <- generatePseudoBulkSeurat(
    pseudobulkList = MetaCellMatrices)

  #plotQC_UpC(seurat_pseudo_bulk,output_figures_filepath)
  rm(MetaCellMatrices)


  seurat_pseudo_bulk <- Seurat::NormalizeData(seurat_pseudo_bulk,normalization.method="RC",scale.factor = 100000)
  return(seurat_pseudo_bulk)
}
#ok so let's pick this parameter automatically
calculate_suggested_number_of_metacell_neighbours <- function(seurat_object,param_min_n_cells){
  UpC <- colSums(seurat_oi@assays$RNA@counts)
  UpCdf <- data.frame(
    cell = colnames(seurat_oi),
    cluster = stringr::str_sub(seurat_oi@meta.data$cluster, start = 1, end = 20),
    UpC = UpC,
    condition = seurat_oi@meta.data$condition
  )

  flag <- FALSE
  this_df <- UpCdf
  while(flag == FALSE){
    median_values <- this_df  %>%
      group_by(cluster, condition) %>%
      summarize(
        median_UpC = median(UpC, na.rm = TRUE),
        count = n())

    min_k <- ceiling(10000/min(median_values$median_UpC))

    median_values$adjusted_counts <- floor(median_values$count/min_k)
    if(min(median_values$adjusted_counts) < param_min_n_cells){
      ind_remove <- which(median_values$adjusted_counts < param_min_n_cells)
      clusters_to_remove <- unique(median_values$cluster[ind_remove])
      this_df <- this_df[-which(this_df$cluster %in% clusters_to_remove),]
    } else {
      #why do I remove one from here? because it includes the meta cell itself
      min_k <- min_k-1
      flag <- TRUE

    }

  }

  # UMI of all groups is above 10,000 (as per PISCES and cite)
  # Calculate the median UpC value by cluster and condition

  return(min_k)
}


#Parameters: dataset ----
min_cells_per_cluster_condition <- 100
species <-  "human"
condition_oi = "ctrl"
condition_reference = "stim"

#Parameters: directories ----
dataset_path <- "sample_analysis"
dir.create(dataset_path)
pre_processing_path <- file.path(dataset_path,"pre_processing")
#this variable is duplicated, remove
pre_processing_filepath <- pre_processing_path
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

#Parameters: analysis ----
flag.normalize.non.log <- FALSE
flag.co.grn <- TRUE
max_n_cells <- 600
param_default_enricher_g_size <- 300
param_min_ligand_expr_in_cluster <- 0.1
param_min_n_cells <- 100
param_max_n_cells <- 600
param_min_receptor_expr_in_cluster <- 0.1
param_n_sample_regulons <- 20

##output objects initialize ----
ligand_scores_result <- list()
decipher_scores_by_regulon_and_cluster <- list()
regulon_scores_by_cluster <- list()
interaction_potential_by_clusters <- list()
regulon_deltas_by_cluster <- list()
significant_regulons_by_cluster <- list()
significant_regulon_markers_by_cluster <- list()
regulon_grns_by_cluster <- list()
lr_markers_by_cluster <- list()
de_markers_by_cluster <- list()
enrichr_results_by_cluster <- list()
interaction_deltas_by_cluster <- list()

#create sample dataset ----
#including seurat object and h5ad objects
seurat_oi <- generateSampleSeuratFromExperimentHub()
##save outputs for Decipher analysis
saveRDS(seurat_oi,file.path("sample_analysis/pre_processing","seurat_object_oi.rds"))
#in addition, we need to create python-compatible h5ad objects for the CO pipeline, here, I've opted against it
#as they are not necessary for this script
# for(this_cluster in unique(kang.seurat$cluster)){
#   seurat_object_oi_this_cluster <- subset(kang.seurat,subset = cluster == this_cluster)
#   sce.object = as.SingleCellExperiment(seurat_object_oi_this_cluster)
#   sce.object@assays@data[["logcounts"]] <- NULL
#   writeH5AD(sce.object, file.path("pre_processing/h5ad_by_cluster",paste(this_cluster,".h5ad",sep="")),X_name = "counts")
# }

#load main data ----
seurat_oi <- readRDS(file.path(pre_processing_filepath,"seurat_object_oi.rds"))

#load reference data ----
L.set <- getForrestLRDatabase(file.path(reference_filepath,"connectomedb_forrest_lrc2p.csv"))
L.set <- L.set %>% mutate(interaction = paste(ligand,receptor,sep="-"),
                          lr = interaction) %>% unique()

cytosig_ligands <- readRDS(file.path(reference_filepath,"cytosig_ligands_human.rds"))

if(species == "human"){
  enrichr_database <- readRDS(file.path(reference_filepath,"enrichr_database_human.rds"))
} else if (species == "mouse"){
  #TODO: check if database in enrichr_database_mouse.rds or enrichr_database_mouse_custom.rds
  enrichr_database <- readRDS(file.path(reference_filepath,"enrichr_database_mouse.rds"))
  cytosig_ligands <- convertHumanSymbolsToMouse(cytosig_ligands)
  L.set <- convertLsetToMouse(L.set)
} else {
  stop("Decipher currently only supports human and mouse species.")
}

#data pre-processing ----
#define case and control
case_condition <- condition_oi
control_condition <- condition_reference

#retain original condition
seurat_oi@meta.data$orig.condition <- seurat_oi@meta.data$condition
#map conditions to case and control
seurat_oi <- mapConditionsInSeurat(seurat_oi,"condition",case_condition,control_condition)

##############
##QC ----
##############
CpC_data <- generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))
plotQC_CpC(CpC_data,outputPath=output_figures_filepath)
param_min_CpC <- 100
#PARAM: select the minimum number of cells per cluster + condition
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data,param_min_CpC)

seurat_oi <- subset(seurat_oi,subset = cluster %in% clusters_passing_CpC_filter)
##############
##Meta cells ----
##############
decipher_seurat <- metaCellModule(
  seurat_object = seurat_oi,
  min_meta_cells = param_min_n_cells
)

saveRDS(decipher_seurat,file.path(output_data_filepath,"pseudobulk_seurat.rds"))

##############
#data pre-processing: main analysis ----
##############
DefaultAssay(decipher_seurat) <- "RNA"
Idents(decipher_seurat) <- decipher_seurat@meta.data$cluster

receptors <- unique(L.set$receptor)
ligands <- unique(L.set$ligand)
all_ligand_receptors <- unique(c(ligands,receptors))

#OUTPUT:
decipher_seurat_lr <- subset(decipher_seurat,features = all_ligand_receptors)

#OUTPUT:
feature_statistics <- getFeatureStatistics(
  features=all_ligand_receptors,
  seuratObj=decipher_seurat)

expressed_ligands <- feature_statistics %>%
  filter(feature %in% ligands & frac.cells.w.counts > param_min_ligand_expr_in_cluster) %>%
  pull(feature) %>%
  unique()

ind_case <- which(decipher_seurat$condition == "case")
data_decipher_seurat_case <- decipher_seurat@assays$RNA@data[,ind_case]
data_decipher_seurat_control <- decipher_seurat@assays$RNA@data[,-ind_case]

#run DECIPHER pipeline-----
all_clusters <- unique(decipher_seurat$cluster)

for(this_cluster in all_clusters){
  print(paste("calculating scores for",this_cluster))
  this_cluster.rds <- paste(this_cluster,"rds",sep=".")
  this_cluster.csv <- paste(this_cluster,"csv",sep=".")

  ###reference objects ----
  # interactions
  selected_receptors <- feature_statistics %>%
    dplyr::filter(cluster == this_cluster) %>%
    dplyr::filter(frac.cells.w.counts > param_min_receptor_expr_in_cluster) %>%
    pull(feature) %>%
    unique()

  L_set_relevant_features <- L.set %>%
    filter(receptor %in% selected_receptors & ligand %in% expressed_ligands)

  # regulons
  if(!file.exists(file.path(output_filepath,"../celloracle/data/GRN",this_cluster.csv))){
    warning("no cell oracle grn found for this cluster")
    next
  }

  regulon_this_cluster <- read.csv(file.path(output_filepath,"../celloracle/data/GRN",this_cluster.csv))
  regulon_this_cluster <- regulon_this_cluster[,-1]

  if(flag.co.grn){
    #trim GRN using CellOracle results
    regulon_this_cluster = trimGRN(
      grn_df = regulon_this_cluster,
      pValue = 0.01,
      topEdges = 20000,
      minTargets = 20)
  }

  regulon_this_cluster_capped <- capRegulon(regulon_this_cluster,n_top = 40)
  #regulon_this_cluster_capped_2 <- capRegulon_2(regulon_this_cluster,n_top = 40)

  random_grns <- generateRandomGRNsFromReference(
    all_genes = rownames(decipher_seurat),
    reference_grns = regulon_this_cluster_capped
  )

  regulon_this_cluster_capped <- rbind(regulon_this_cluster_capped,random_grns)


  ### primary object ----
  decipher_seurat_this_cluster <- subset(decipher_seurat,subset = cluster == this_cluster)

  if(flag.normalize.non.log){
    decipher_seurat_this_cluster <- NormalizeData(decipher_seurat_this_cluster,normalization.method = "RC",scale.factor=100000)
  }

  base_n_cells <- min(table(decipher_seurat_this_cluster$condition))
  if(base_n_cells < param_min_n_cells){
    next
  }
  if(base_n_cells > param_max_n_cells){
    base_n_cells <- param_max_n_cells
  }

  Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster@meta.data$condition
  decipher_seurat_this_cluster_downsampled <- subset(decipher_seurat_this_cluster, downsample = base_n_cells)

  data_this_cluster_downsampled <- decipher_seurat_this_cluster_downsampled@assays$RNA@data

  #dplyr doesn't work with sparse matrix
  #data_this_cluster_downsampled_receptors <- data_this_cluster_downsampled[which(rownames(data_this_cluster_downsampled) %in% selected_receptors),]
  data_this_cluster_downsampled_receptors <- data_this_cluster_downsampled[which(rownames(data_this_cluster_downsampled) %in% unique(L_set_relevant_features$receptor)),]
  ##PAGODA -----
  #silence this function
  regulon_scores_this_cluster <- getRegulonScores(
    seuratObject = decipher_seurat_this_cluster,
    grn_df = regulon_this_cluster_capped)

  ##PAGODA DELTA ----
  regulon_deltas_this_cluster <- getRegulonDeltas(
    regulon_scores_this_cluster,
    decipher_seurat_this_cluster$condition)

  regulon_deltas_this_cluster <- regulon_deltas_this_cluster %>%
    mutate(class = ifelse(stringr::str_detect(name,"sample"),"random","real"))

  random.density <- density(subset(regulon_deltas_this_cluster, class == "random")$deltaPagoda, n = 2^10)
  upper_threshold_random <- quantile(random.density,probs = 0.975,normalize = FALSE)
  lower_threshold_random <- quantile(random.density,probs = 0.025,normalize = FALSE)

  significant_regulon_deltas_this_cluster <-  regulon_deltas_this_cluster %>%
    filter(class == "real") %>%
    filter(deltaPagoda > upper_threshold_random | deltaPagoda < lower_threshold_random)%>%
    arrange(deltaPagoda)

  print(paste(this_cluster,"number of significant regulons:",nrow(significant_regulon_deltas_this_cluster)))
  #### find target genes for each top differentially expressed regulons and calculate diff expr. ----
  Idents(decipher_seurat_this_cluster) <- decipher_seurat_this_cluster$condition
  #wait but this needs to align to my GRN right?
  significant_regulon_markers_by_cluster[[this_cluster]] <- getDifferentiallyExpressedTargetsForRegulons(
    seuratObj = decipher_seurat_this_cluster,
    regulonNames = significant_regulon_deltas_this_cluster$name,
    logFcThreshold = 0.58,
    grnDf = regulon_this_cluster,
    targetCt = this_cluster
  )

  ## calculate Interaction Potential Matrix ----
  interaction_potentials <- list()
  interaction_mapping_table <-  getInteractionMappingTable(
    receptorMatrix = data_this_cluster_downsampled_receptors,
    ligandSet = L_set_relevant_features
    )

  #case and control
  condition_vector <- decipher_seurat_this_cluster_downsampled$condition
  unique_ligands <- unique(L_set_relevant_features$ligand)
  #ligand_means is a data frame with a case and a control column
  ligand_means <- data.frame(
    case = Matrix::rowMeans(data_decipher_seurat_case[unique_ligands,]),
    control = Matrix::rowMeans(data_decipher_seurat_control[unique_ligands,]),
    row.names = unique_ligands
  )

  interaction_potentials_matrix_this_cluster <- calculateInteractionMatrix(
    receptorMatrix = data_this_cluster_downsampled_receptors,
    conditionVector = condition_vector,
    ligandMeans = ligand_means,
    LRSet = L_set_relevant_features
  )

  #remove zeros
  ind_no_information <- which(rowSums(interaction_potentials_matrix_this_cluster) == 0)
  if(length(ind_no_information) > 0){
    interaction_potentials_matrix_this_cluster <- interaction_potentials_matrix_this_cluster[-ind_no_information,]
  }

  new_seurat <- Seurat::CreateSeuratObject(counts = interaction_potentials_matrix_this_cluster,meta.data = decipher_seurat_lr@meta.data[colnames(interaction_potentials_matrix_this_cluster),])
  Idents(new_seurat) <- new_seurat$condition
  interaction_deltas <- FindMarkers(new_seurat,ident.1 = "case",logfc.threshold = 0.1)
  interaction_deltas <- interaction_deltas %>%
    filter(p_val_adj < 0.01)
  interaction_deltas$name = rownames(interaction_deltas)
  interaction_deltas_by_cluster[[this_cluster]] <- interaction_deltas


  interaction_potentials_matrix_this_cluster <- interaction_potentials_matrix_this_cluster[rownames(interaction_deltas),]
  ## subset interaction_potential matrix by correlation clusters for cluster-based RF ----
  #split matrix into interactions comprised of receptors with a unique ligand (one-to-one), and interactions of receptors with multiple ligands (many-to-one)
  one_to_one_interactions <- intersect(getOneToOneInteractions(interaction_mapping_table),rownames(interaction_potentials_matrix_this_cluster))
  many_to_one_interactions <- intersect(getManyToOneInteractions(interaction_mapping_table),rownames(interaction_potentials_matrix_this_cluster))

  interaction_potentials_matrix_this_cluster_oto <- interaction_potentials_matrix_this_cluster[one_to_one_interactions,]
  interaction_potentials_matrix_this_cluster_mto <- interaction_potentials_matrix_this_cluster[many_to_one_interactions,]


  ## correlation clusters for many-to-one interactions ----
  mto_interactions_clusters <- getCorrelationClusters(
    interactionPotentialsMatrixMTO = interaction_potentials_matrix_this_cluster_mto,
    interactionMappingTable = interaction_mapping_table,
    pctMTOReceptors = 1.15,
    correlationMethod = "spearman",
    clusteringMethod = "complete")

  ## representative interaction for each cluster ----
  representative_interactions_mto <- getRepresentativeInteractionsForMTOClusters(
    mtoInteractionsClusters = mto_interactions_clusters,
    interactionMappingTable = interaction_mapping_table,
    prioritizedBenchmarkingLigands = cytosig_ligands
  )

  ## cluster-based matrix from random forest -----
  interaction_potentials_matrix_this_cluster_mto_representative <- interaction_potentials_matrix_this_cluster[representative_interactions_mto$interaction,]
  interaction_potentials_matrix_clusters <- rbind(interaction_potentials_matrix_this_cluster_oto,interaction_potentials_matrix_this_cluster_mto_representative)

  #constrain analysis to differentially-expressed regulons
  regulon_differential_expression <- significant_regulon_deltas_this_cluster

  ## run random forest on each regulon -----
  all_rf_results <- list()
  for(this.tf in significant_regulon_deltas_this_cluster$name){
    ind.this.tf <- which(significant_regulon_deltas_this_cluster$name == this.tf)
    val.this.tf <- significant_regulon_deltas_this_cluster$deltaPagoda[ind.this.tf]
    print(paste("calculating forest for",this.tf))
    tf.merged <- regulon_scores_this_cluster[this.tf,colnames(interaction_potentials_matrix_clusters)]
    rf <- randomForest(x = t(interaction_potentials_matrix_clusters),y=tf.merged, ntree = 100,importance=T)
    #print(rf)
    imp.perm.merged <- importance(rf,type=1, scale = F)
    #head(imp.perm.merged)

    spearman.cor <- cor(t(interaction_potentials_matrix_clusters),tf.merged,method = "spearman")
    pearson.cor <- cor(t(interaction_potentials_matrix_clusters),tf.merged,method = "pearson")

    imp <- importance(rf, scale = F)
    index_match_interaction_mapping_table <- match(rownames(imp),interaction_mapping_table$interaction)


    imp.df <- data.frame(
      interaction = interaction_mapping_table$interaction[index_match_interaction_mapping_table],
      ligand =  interaction_mapping_table$ligand[index_match_interaction_mapping_table],
      receptor =  interaction_mapping_table$receptor[index_match_interaction_mapping_table],
      imp.perm = imp[,1],
      perm.rank = length(imp[,1])-rank(imp[,1]),
      imp.gini = imp[,2],
      gini.rank = length(imp[,2])-rank(imp[,2]),
      gene = rownames(imp),
      regulon = this.tf,
      regulon.val = val.this.tf,
      pearson.cor =  pearson.cor,
      spearman.cor = spearman.cor,
      possible.spearman.cont = spearman.cor*val.this.tf,
      weighted.spearman.cont = imp[,1]*sign(spearman.cor)*val.this.tf
    )

    imp.df <- imp.df[order(imp.df$perm.rank,decreasing=FALSE),]

    all_rf_results[[this.tf]] <- imp.df
    #write.csv(imp.df,file.path("data/importances",file=paste(this_cluster,this.tf,"all_importances.csv",sep="_")))
  }

  #convert interaction_potential list into a matrix
  first.flag <- TRUE
  for(this.tf in names(all_rf_results)){
    this.tf.results <- all_rf_results[[this.tf]]
    if(first.flag){
      interaction_weights <- this.tf.results
      first.flag <- FALSE
    }else {
      interaction_weights <- rbind(interaction_weights,this.tf.results)
    }
  }
  #stuff for visualization
  lr_markers_this_cluster <- FindMarkers(decipher_seurat_this_cluster,
                                         ident.1 = "case",
                                         ident.2 = "control",
                                         feature = unique(c(interaction_weights$ligand,interaction_weights$receptor)),
                                         logfc.threshold = 0,
                                         min.pct = 0,
                                         only.pos = FALSE)

  de_markers_this_cluster <- FindMarkers(decipher_seurat_this_cluster,
                                         ident.1 = "case",
                                         ident.2 = "control",
                                         logfc.threshold = 0.58,
                                         only.pos = FALSE)


  ##enrichr ----
  #enrichr on transcription factors
  de_markers_this_cluster$gene <- rownames(de_markers_this_cluster)

  all_pos <- de_markers_this_cluster %>%
    filter(avg_log2FC > 0) %>%
    slice_max(avg_log2FC,n=300)%>%
    select(gene)

  all_neg <- de_markers_this_cluster %>%
    filter(avg_log2FC < 0) %>%
    slice_min(avg_log2FC,n=300)%>%
    select(gene)


  all_pos <- all_pos$gene
  all_neg <- all_neg$gene
  all_pos_neg <- c(all_pos,all_neg)

  my_gene_sets <- list()
  for(this_regulon in significant_regulon_deltas_this_cluster$name){
    all_genes <- regulon_this_cluster$target[regulon_this_cluster$source == this_regulon]
    my_gene_sets[[this_regulon]] <- all_genes
  }
  my_gene_sets[["all_pos"]] <- all_pos
  my_gene_sets[["all_neg"]] <- all_neg
  my_gene_sets[["all_pos_neg"]] <- all_pos_neg

  gene_set_results <- list()
  for(this_gene_set in names(my_gene_sets)){
    dbs_results <- list()
    for(this_dbs_name in names(enrichr_database)){
      term_results <- list()
      this_dbs <- enrichr_database[[this_dbs_name]]
      for(this_term in names(this_dbs)){
        term_results[[this_term]] <- calculateEnrichmentStatistics(this_dbs[[this_term]],my_gene_sets[[this_gene_set]],20000)
      }
      dbs_results[[this_dbs_name]] <- term_results
    }
    gene_set_results[[this_gene_set]] <- dbs_results
  }

  regulon_results_df <- list()
  dbs_results_df <- list()
  for(this_regulon in names(gene_set_results)){
    for(this_dbs_name in names(enrichr_database)){
      dbs_results_df[[this_dbs_name]]  <- do.call(rbind.data.frame, gene_set_results[[this_regulon]][[this_dbs_name]])
      dbs_results_df[[this_dbs_name]]$database  <- rep(this_dbs_name,nrow(dbs_results_df[[this_dbs_name]]))
      dbs_results_df[[this_dbs_name]]$p_value_adjusted  <- p.adjust(dbs_results_df[[this_dbs_name]]$p_value, method = "BH", n = length(dbs_results_df[[this_dbs_name]]$p_value))
      dbs_results_df[[this_dbs_name]]$Term <- rownames(dbs_results_df[[this_dbs_name]])
    }
    regulon_results_df[[this_regulon]] <- do.call(rbind.data.frame, dbs_results_df)
  }
  enrichr_results_by_cluster[[this_cluster]] <- regulon_results_df

  #save all intermediate data
  regulon_grns_by_cluster[[this_cluster]] <- regulon_this_cluster
  regulon_scores_by_cluster[[this_cluster]] <- regulon_scores_this_cluster
  regulon_deltas_by_cluster[[this_cluster]] <- regulon_deltas_this_cluster
  #TODO: uncomment :)
  interaction_potential_by_clusters[[this_cluster]] <- interaction_potentials_matrix_this_cluster
  decipher_scores_by_regulon_and_cluster[[this_cluster]] <- interaction_weights
  lr_markers_by_cluster[[this_cluster]] <- lr_markers_this_cluster
  de_markers_by_cluster[[this_cluster]] <- de_markers_this_cluster
  significant_regulons_by_cluster[[this_cluster]] <- significant_regulon_deltas_this_cluster


  }

decipher_scores_by_regulon_and_cluster <- lapply(decipher_scores_by_regulon_and_cluster,FUN = "listOfDFsRenameColumn",original_name = "weighted.spearman.cont",new_name = "decipher_score")
decipher_scores_by_cluster <- lapply(decipher_scores_by_regulon_and_cluster,FUN = "calculateScoresByCluster")
decipher_scores_by_cluster <- addListNameToDFElements(decipher_scores_by_cluster,"receiver_cluster")

#save DECIPHER ----
saveRDS(regulon_scores_by_cluster,file.path(output_data_filepath,"regulon_scores_by_cluster.rds"))
saveRDS(decipher_scores_by_regulon_and_cluster,file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(interaction_potential_by_clusters,file.path(output_data_filepath,"interaction_potential_by_clusters.rds"))
saveRDS(regulon_deltas_by_cluster,file.path(output_data_filepath,"regulon_deltas_by_cluster.rds"))
saveRDS(significant_regulons_by_cluster,file.path(output_data_filepath,"significant_regulons_by_cluster.rds"))
saveRDS(significant_regulon_markers_by_cluster,file.path(output_data_filepath,"significant_regulon_markers_by_cluster.rds"))
saveRDS(regulon_grns_by_cluster,file.path(output_data_filepath,"regulon_grns_by_cluster.rds"))
saveRDS(lr_markers_by_cluster,file.path(output_data_filepath,"lr_markers_by_cluster.rds"))
saveRDS(de_markers_by_cluster,file.path(output_data_filepath,"de_markers_by_cluster.rds"))
saveRDS(enrichr_results_by_cluster,file.path(output_data_filepath,"enrichr_results_by_cluster.rds"))
saveRDS(feature_statistics,file.path(output_data_filepath,"feature_statistics.rds"))
saveRDS(decipher_seurat_lr,file.path(output_data_filepath,"decipher_seurat_lr.rds"))
saveRDS(L.set,file.path(output_data_filepath,"L_set.rds"))
saveRDS(decipher_scores_by_regulon_and_cluster, file.path(output_data_filepath,"decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(decipher_scores_by_cluster,file.path(output_data_filepath,"decipher_scores_by_cluster.rds"))
saveRDS(interaction_deltas_by_cluster,file.path(output_data_filepath,"interaction_deltas_by_cluster.rds"))


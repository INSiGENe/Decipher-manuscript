library(R6)

DecipherObject <- R6Class(
  "DecipherObject",
  
  public = list(
    
    #---------------------------------------------------------------------
    # 1) USER-DEFINED PARAMETERS AND PATHS
    #---------------------------------------------------------------------
    dataset_name = NULL,
    dataset_path = NULL,
    condition_name = NULL,
    case_condition = NULL,
    control_condition = NULL,
    k_parameter = NULL,
    species = NULL,
    
    # optional: you may store the random seed in the object
    # so you can refer to it if you want reproducible runs
    selected_random_seed = NULL,
    
    # directories
    pre_processing_path = NULL,
    reference_filepath = NULL,
    output_filepath = NULL,
    output_data_filepath = NULL,
    output_figures_filepath = NULL,
    output_importances_filepath = NULL,
    
    #---------------------------------------------------------------------
    # 2) INTERMEDIATE OBJECTS
    #---------------------------------------------------------------------
    seurat_oi = NULL,
    decipher_seurat = NULL,
    decipher_seurat_lr = NULL,
    
    # Example intermediate pipeline outputs:
    CpC_data = NULL,
    clusters_passing_CpC_filter = NULL,
    L_set = NULL,
    enrichr_database = NULL,
    cytosig_ligands = NULL,
    # ... etc. You can store more as you go along ...
    
    #---------------------------------------------------------------------
    # 3) INITIALIZATION
    #---------------------------------------------------------------------
    initialize = function(dataset_path = "results/BCG",
                          condition_name = "condition",
                          case_condition = "D21",
                          control_condition = "D0",
                          k_parameter = 2,
                          species = "human",
                          selected_random_seed = 123) {
      
      # store parameters
      self$dataset_path        <- dataset_path
      self$condition_name      <- condition_name
      self$case_condition      <- case_condition
      self$control_condition   <- control_condition
      self$k_parameter         <- k_parameter
      self$species             <- species
      self$selected_random_seed <- selected_random_seed
      
      # set the seed if reproducibility is needed
      if (!is.null(selected_random_seed)) {
        set.seed(selected_random_seed)
      }
      
      # set up directories
      self$pre_processing_path         <- file.path(dataset_path, "pre_processing")
      self$reference_filepath          <- "reference_data"
      self$output_filepath             <- dataset_path
      self$output_data_filepath        <- file.path(self$output_filepath, "data")
      self$output_figures_filepath     <- file.path(self$output_filepath, "figures")
      self$output_importances_filepath <- file.path(self$output_filepath, "importances")
      
      # create directories if they don't exist
      dir.create(dataset_path, showWarnings = FALSE)
      dir.create(self$pre_processing_path, showWarnings = FALSE)
      dir.create(file.path(self$pre_processing_path, "h5ad_by_cluster"), showWarnings = FALSE)
      dir.create(self$output_data_filepath, recursive = TRUE, showWarnings = FALSE)
      dir.create(self$output_figures_filepath, recursive = TRUE, showWarnings = FALSE)
      dir.create(self$output_importances_filepath, recursive = TRUE, showWarnings = FALSE)
    },
    
    #---------------------------------------------------------------------
    # 4) METHODS FOR EACH STAGE OF THE PIPELINE
    #---------------------------------------------------------------------
    
    #---------------------------
    # 4a) LOAD THE DATASET
    #---------------------------
    loadDataset = function() {
      # e.g. read in your Seurat object from pre-processing
      message("Loading Seurat object from: ", file.path(self$pre_processing_path,"seurat_object_oi.rds"))
      self$seurat_oi <- readRDS(file.path(self$pre_processing_path,"seurat_object_oi.rds"))
      
      # if you have extra steps like:
      self$seurat_oi$orig.condition <- self$seurat_oi[[self$condition_name]]
      self$seurat_oi <- mapConditionsInSeurat(
        self$seurat_oi, 
        self$condition_name,
        self$case_condition,
        self$control_condition
      )
      return(invisible(self))
    },
    
    #---------------------------
    # 4b) LOAD REFERENCE DATA
    #---------------------------
    loadReferences = function() {
      # example references
      self$L_set <- loadLSet(self$reference_filepath, self$species)
      self$enrichr_database <- loadEnrichrDatabase(self$reference_filepath, self$species)
      self$cytosig_ligands  <- loadCytosigLigands(self$reference_filepath, self$species)
      
      return(invisible(self))
    },
    
    #---------------------------
    # 4c) RUN QC
    #---------------------------
    runQC = function(min_cells_per_cluster_condition = 100) {
      # generate data for QC by cluster and condition
      self$CpC_data <- generateQCDataByClusterAndCondition(
        self$seurat_oi,
        max(stringr::str_length(unique(self$seurat_oi$cluster)))
      )
      
      plotQC_CpC(self$CpC_data, outputPath = self$output_figures_filepath)
      
      # filter clusters
      self$clusters_passing_CpC_filter <- getClustersPassingCpCFilter(
        self$CpC_data, 
        minCpc = min_cells_per_cluster_condition
      )
      
      # subset the seurat object for passing clusters
      self$seurat_oi <- self$seurat_oi[
        , which(self$seurat_oi$cluster %in% self$clusters_passing_CpC_filter),
        seed = NULL
      ]
      
      plotQC_UpC(self$seurat_oi, 
                 outputPath = self$output_figures_filepath, 
                 id = "_sc")
      
      return(invisible(self))
    },
    
    #---------------------------
    # 4d) META CELL MODULE
    #---------------------------
    runMetaCells = function(min_meta_cells_parameter = 100) {
      
      self$decipher_seurat <- metaCellModule(
        seurat_object = self$seurat_oi,
        min_meta_cells = min_meta_cells_parameter,
        k = self$k_parameter
      )
      
      plotQC_UpC(
        seuratObject = self$decipher_seurat, 
        outputPath = self$output_figures_filepath, 
        id = "_meta"
      )
      
      # optional: generate QC data for meta
      CpC_data_meta <- generateQCDataByClusterAndCondition(
        self$decipher_seurat,
        max(stringr::str_length(unique(self$decipher_seurat$cluster)))
      )
      plotQC_CpC(CpC_data_meta, outputPath = self$output_figures_filepath, id = "_meta")
      
      # record the parameters
      parameter_record <- data.frame(
        "k" = self$k_parameter,
        "min_meta_cells" = min_meta_cells_parameter
      )
      
      write.csv(
        parameter_record,
        file.path(self$output_data_filepath, "parameter_record.csv")
      )
      
      # save pseudobulk seurat
      saveRDS(
        self$decipher_seurat,
        file.path(self$output_data_filepath, "pseudobulk_seurat.rds")
      )
      
      return(invisible(self))
    },
    
    #---------------------------
    # 4e) DATA PRE-PROCESSING (MAIN ANALYSIS)
    #---------------------------
    runMainAnalysis = function(flag.normalize.non.log = FALSE) {
      
      # subset to ligands and receptors
      self$decipher_seurat_lr <- self$decipher_seurat[
        unique(c(self$L_set$ligand, self$L_set$receptor)),
        , seed = NULL
      ]
      
      feature_statistics <- getFeatureStatistics(
        features = unique(c(self$L_set$ligand, self$L_set$receptor)),
        seuratObj = self$decipher_seurat
      )
      
      expressed_ligands <- getFilteredLigands(
        self$decipher_seurat,
        self$L_set,
        param_min_ligand_expr_in_cluster = 0.1
      )
      
      expressed_receptors_all_clusters <- getExpressedReceptorsForEachCluster(
        self$decipher_seurat,
        self$L_set
      )
      
      L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(
        self$L_set,
        expressed_ligands,
        expressed_receptors_all_clusters
      )
      
      regulon_grns_by_cluster <- getRegulonsAllClusters(
        self$output_filepath,
        self$decipher_seurat
      )
      
      capped_regulons_all_clusters <- capRegulonsAllClusters(
        regulon_grns_by_cluster,
        self$decipher_seurat,
        flag.normalize.non.log
      )
      
      regulon_scores_by_cluster <- getRegulonScoresAllClusters(
        capped_regulons_all_clusters,
        self$decipher_seurat
      )
      
      regulon_deltas_by_cluster <- getRegulonDeltasAllClusters(
        regulon_scores_by_cluster,
        self$decipher_seurat
      )
      
      significant_regulons_by_cluster <- getSignificantRegulonsAllClusters(
        regulon_deltas_by_cluster
      )
      
      significant_regulon_markers_by_cluster <- getDifferentiallyExpressedTargetsForRegulonsAllClusters(
        self$decipher_seurat,
        significant_regulons_by_cluster,
        regulon_grns_by_cluster,
        flag.normalize.non.log,
        random.seed = self$selected_random_seed
      )
      
      interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClusters(
        self$decipher_seurat,
        L_set_relevant_features_all_clusters,
        flag.normalize.non.log
      )
      
      interaction_deltas_by_cluster <- calculateInteractionDeltasAllClusters(
        interaction_potential_by_clusters,
        self$decipher_seurat_lr
      )
      
      filtered_interaction_potentials_matrix_all_clusters <- filterIntPotByDeltas(
        interaction_potential_by_clusters,
        interaction_deltas_by_cluster
      )
      
      interaction_potentials_matrix_clusters_all_clusters <-
        getInteractionPotentialMatrixForRepresentativeInteractionsAllClusters(
          self$decipher_seurat,
          L_set_relevant_features_all_clusters,
          filtered_interaction_potentials_matrix_all_clusters,
          self$cytosig_ligands,
          flag.normalize.non.log
        )
      
      decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClusters(
        self$decipher_seurat,
        significant_regulons_by_cluster,
        regulon_scores_by_cluster,
        interaction_potentials_matrix_clusters_all_clusters,
        L_set_relevant_features_all_clusters,
        flag.normalize.non.log
      )
      
      lr_markers_by_cluster <- FindLRMarkersAllClusters(
        self$decipher_seurat,
        decipher_scores_by_regulon_and_cluster,
        flag.normalize.non.log,
        random.seed = self$selected_random_seed
      )
      
      de_markers_by_cluster <- FindMarkersAllClusters(
        self$decipher_seurat,
        flag.normalize.non.log,
        random.seed = self$selected_random_seed
      )
      
      # optional: if you want the Enrichr step:
      # enrichr_results_by_cluster <- enrichResultsAllClusters(
      #   de_markers_by_cluster,
      #   significant_regulons_by_cluster,
      #   regulon_grns_by_cluster,
      #   self$enrichr_database
      # )
      
      # rename columns for DECIPHER
      decipher_scores_by_regulon_and_cluster <- lapply(
        decipher_scores_by_regulon_and_cluster,
        FUN = "listOfDFsRenameColumn",
        original_name = "weighted.spearman.cont",
        new_name = "decipher_score"
      )
      
      decipher_scores_by_cluster <- lapply(
        decipher_scores_by_regulon_and_cluster,
        FUN = "calculateScoresByCluster"
      )
      
      decipher_scores_by_cluster <- addListNameToDFElements(
        decipher_scores_by_cluster,
        "receiver_cluster"
      )
      
      # final plot
      plotDecipherPrioritizedMap(
        self$dataset_path,
        top_n = 6,
        dataset_name = dataset_name
      )
      
      # Optionally save the enrichr step
      # saveRDS(enrichr_results_by_cluster, 
      #         file.path(self$output_data_filepath, "enrichr_results_by_cluster.rds"))
      
      return(invisible(self))
    }
  )
)

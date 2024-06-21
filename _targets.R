# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  #Decipher here (#0)
  packages = c("tibble","randomForest","pagoda2","Seurat","babelgene","magrittr","tibble","stringr","Matrix","SeuratObject","tidyr","scde","ExperimentHub","dplyr","ggplot2"), # Packages that your targets need for their tasks.
  seed = 123
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
lapply(list.files("R", full.names = TRUE), source)
# tar_source("other_functions.R") # Source other scripts as needed.

#parameters for pipeline
#TODO: Figure out a more user-friendly way/adaptable way to add these

# Define parameters
 min_cells_per_cluster_condition <- 100  # Minimum number of cells per cluster condition
 species <- "human"                      # Species of the sample
 case_condition <- "stim"                # Case condition (e.g., stimulated)
 control_condition <- "ctrl"             # Control condition (e.g., control)
#
# # Define analysis parameters
 flag.normalize.non.log <- FALSE
# Define directory parameters
 reference_filepath <- "reference_data"

# Define the target list
list(
  # Load parameters
  tar_target(min_cells_per_cluster_condition, 100),
  tar_target(species, "human"),
  tar_target(case_condition, "stim"),
  tar_target(control_condition, "ctrl"),
  tar_target(reference_filepath, "reference_data"),
  tar_target(flag.normalize.non.log, FALSE),

  #load data
  tar_target(seurat_oi, generateSampleSeuratFromExperimentHub(min_cells_per_cluster_condition,case_condition,control_condition)),
  tar_target(L.set,loadLSet(reference_filepath,species)),
  tar_target(enrichr_database,loadEnrichrDatabase(reference_filepath,species)),
  tar_target(cytosig_ligands,loadCytosigLigands(reference_filepath,species)),

  #QC visuals
  tar_target(CpC_data,generateQCDataByClusterAndCondition(seurat_oi,max(stringr::str_length(unique(seurat_oi$cluster))))),
  tar_target(CpC_plot,plotQC_CpC(CpC_data,outputPath=reference_filepath)),
  tar_target(clusters_passing_CpC_filter,getClustersPassingCpCFilter(CpC_data,minCpc = 100)),
  tar_target(seurat_oi_CpC,subset(seurat_oi,subset = cluster %in% clusters_passing_CpC_filter)),

  #meta cells
  tar_target(seurat_oi_meta,metaCellModule(
    seurat_object = seurat_oi_CpC,
    min_meta_cells = 100
  )),
  tar_target(decipher_seurat_lr,subset(seurat_oi_meta,features = unique(c(L.set$ligand,L.set$receptor)))),

  tar_target(feature_statistics,getFeatureStatistics(
    features=unique(c(L.set$ligand,L.set$receptor)),
    seuratObj=seurat_oi_meta)),

  tar_target(expressed_ligands,getFilteredLigands(
    seurat_oi_meta,
    L.set,
    param_min_ligand_expr_in_cluster = 0.1)
  )


)

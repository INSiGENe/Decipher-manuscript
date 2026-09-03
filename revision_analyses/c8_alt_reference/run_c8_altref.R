# R2.5 robustness — rerun the C8 (CD14+BDCA1+PD-L1+) Decipher analysis with an
# ALTERNATIVE monocyte reference population as the control arm.
#
# Usage: Rscript run_c8_altref.R <control_cluster>
#   control_cluster ∈ CD14_plus_monocytes (baseline replication),
#                     CD16_plus_monocytes, cDC2
#
# Mirrors scripts/analysis_specific_datasets/covid/covid_2_decipher_pipeline_v1_modularized.R
# (seed 123, k=2, same functions sourced from the repo) with two deliberate changes:
#   1. paramPairings row for C8 gets `control_cluster` instead of CD14+ monocytes
#      (C8 is the LAST pairing row, so all other clusters' meta-cell draws are
#      identical to the manuscript run).
#   2. Downstream per-cluster scoring is restricted to the C8 pairing (the full
#      decipher_seurat, with every cluster, is still built first so ligand/receptor
#      expression pooling matches the original pipeline).
# FindMarkers side-outputs not needed for the comparison (significant regulon
# markers, LR markers, DE markers) are skipped.
#
# Seurat 5 compat: repo files are sourced through a text patch replacing the
# defunct `slot = "` argument with `layer = "` — behaviour-identical for v3
# assays. Outputs land in outputs/<control_cluster>/.

repo_path  <- "/Users/edgarbasto/Documents/Decipher-manuscript"
drive_covid <- "/Volumes/MegaEdgar/aws_pull_20260813/results/covid"
base_dir   <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/c8_alt_reference"
C8 <- "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells"

args <- commandArgs(TRUE)
stopifnot(length(args) >= 1)
control_cluster <- args[1]
out_dir <- file.path(base_dir, "outputs", control_cluster)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# working dir: results/covid layout expected by getRegulonsAllClusters
work <- file.path(base_dir, "work")
dir.create(file.path(work, "results"), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(work, "results", "covid"))) {
  dir.create(file.path(work, "results", "covid", "cellOracle", "data"), recursive = TRUE)
  file.symlink(file.path(drive_covid, "cellOracle", "data", "GRN"),
               file.path(work, "results", "covid", "cellOracle", "data", "GRN"))
}
setwd(work)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(dplyr); library(tidyr)
  library(tibble); library(stringr); library(igraph); library(magrittr)
  library(Matrix); library(pagoda2); library(randomForest); library(progress)
})
# repo code does direct @assays$RNA@data slot access — build v3 Assays like the
# original pipeline did, not Seurat 5's Assay5
options(Seurat.object.assay.version = "v3")

source_patched <- function(f) {
  txt <- readLines(f, warn = FALSE)
  txt <- gsub('slot = "', 'layer = "', txt, fixed = TRUE)
  eval(parse(text = paste(txt, collapse = "\n")), envir = globalenv())
}
for (f in c("utilities.R", "pre_processing_qc.R", "pre_processing_meta_cells.R",
            "intercellular_tailoring.R", "intercellular_scoring.R",
            "intracellular_tailoring.R", "intracellular_scoring.R",
            "decipher_scoring.R", "integrated_scoring.R")) {
  source_patched(file.path(repo_path, "R", f))
}
stopifnot(exists("generateMetaCellMatricesWPairings"),
          exists("getRandomForestWeightsAllClustersWParamPairings"),
          exists("getRegulonScoresAllClustersWParamPairings"))

message(sprintf("[%s] control_cluster = %s", format(Sys.time(), "%H:%M:%S"), control_cluster))

# ---- pipeline (covid_2 verbatim unless commented) ----
selected_random_seed <- 123
set.seed(selected_random_seed)

min_cells_per_cluster_condition <- 100
species <- "human"
condition_name <- "condition"
case_condition <- 22
control_condition <- 0
k_parameter <- 2
reference_filepath <- file.path(repo_path, "reference_data")
flag.normalize.non.log <- FALSE

seurat_oi <- readRDS(file.path(drive_covid, "pre_processing", "seurat_object_oi.rds"))
seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi, condition_name, case_condition, control_condition)

L.set <- loadLSet(reference_filepath, species)
cytosig_ligands <- loadCytosigLigands(reference_filepath, species)

CpC_data <- generateQCDataByClusterAndCondition(seurat_oi, max(stringr::str_length(unique(seurat_oi$cluster))))
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data, minCpc = 100)
clusters_passing_CpC_filter <- c(clusters_passing_CpC_filter, C8)
seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed = NULL]

min_meta_cells_parameter <- 100

paramPairings <- data.frame(
  case = clusters_passing_CpC_filter,
  control = clusters_passing_CpC_filter
)
paramPairings <- paramPairings %>%
  mutate(control = if_else(case == C8, control_cluster, control))   # <-- the swap

min_counts <- seurat_oi@meta.data %>% group_by(cluster) %>% count(condition) %>% ungroup()
groups <- createGroupsFromPairings(paramPairings)
paramPairings_min_n <- calculateMinimumN(groups, min_counts, paramPairings)

seurat_oi@meta.data$barcode <- rownames(seurat_oi@meta.data)
message(sprintf("[%s] meta-cell generation...", format(Sys.time(), "%H:%M:%S")))
MetaCellMatrices <- generateMetaCellMatricesWPairings(
  seuratObj = seurat_oi,
  paramMaxScCells = 1200 * (k_parameter + 1),
  paramK = k_parameter,
  paramPairings = paramPairings_min_n)

seurat_pseudo_bulk <- generatePseudoBulkSeurat(pseudobulkList = MetaCellMatrices)
decipher_seurat <- Seurat::NormalizeData(seurat_pseudo_bulk, normalization.method = "RC", scale.factor = 100000)
rm(MetaCellMatrices, seurat_pseudo_bulk, seurat_oi); gc(verbose = FALSE)
DefaultAssay(decipher_seurat) <- "RNA"

message(sprintf("[%s] pseudobulk: %d meta-cells", format(Sys.time(), "%H:%M:%S"), ncol(decipher_seurat)))
saveRDS(decipher_seurat@meta.data, file.path(out_dir, "pseudobulk_meta.rds"))

decipher_seurat_lr <- decipher_seurat[unique(c(L.set$ligand, L.set$receptor)), , seed = NULL]
feature_statistics <- getFeatureStatistics(
  features = unique(c(L.set$ligand, L.set$receptor)), seuratObj = decipher_seurat)
expressed_ligands <- getFilteredLigands(decipher_seurat, L.set, param_min_ligand_expr_in_cluster = 0.1)
expressed_receptors_all_clusters <- getExpressedReceptorsForEachCluster(decipher_seurat, L.set)
L_set_relevant_features_all_clusters <- getRelevantFeaturesForEachCluster(
  L.set, expressed_ligands, expressed_receptors_all_clusters)

message(sprintf("[%s] regulons (GRN load + cap)...", format(Sys.time(), "%H:%M:%S")))
regulon_grns_by_cluster <- getRegulonsAllClusters(file.path("results", "covid"), decipher_seurat)
capped_regulons_all_clusters <- capRegulonsAllClusters(
  regulon_grns_by_cluster, decipher_seurat, flag.normalize.non.log)

# restrict downstream to the C8 pairing
pp_c8 <- paramPairings %>% filter(case == C8)
stopifnot(nrow(pp_c8) == 1)

message(sprintf("[%s] PAGODA regulon scoring (C8 pairing)...", format(Sys.time(), "%H:%M:%S")))
regulon_scores_by_cluster <- getRegulonScoresAllClustersWParamPairings(
  capped_regulons_all_clusters, decipher_seurat, paramPairings = pp_c8)
regulon_deltas_by_cluster <- getRegulonDeltasAllClustersWParamPairings(
  regulon_scores_by_cluster, decipher_seurat, pp_c8)
significant_regulons_by_cluster <- getSignificantRegulonsAllClusters(regulon_deltas_by_cluster)

message(sprintf("[%s] interaction potentials...", format(Sys.time(), "%H:%M:%S")))
interaction_potential_by_clusters <- getInteractionPotentialsMatrixAllClustersWParamPairings(
  decipher_seurat, L_set_relevant_features_all_clusters, flag.normalize.non.log, pp_c8)
interaction_deltas_by_cluster <- calculateInteractionDeltasAllClusters(
  interaction_potential_by_clusters, decipher_seurat_lr)
filtered_interaction_potentials_matrix_all_clusters <- filterIntPotByDeltas(
  interaction_potential_by_clusters, interaction_deltas_by_cluster)
interaction_potentials_matrix_clusters_all_clusters <-
  getInteractionPotentialMatrixForRepresentativeInteractionsAllClustersWParamPairings(
    decipher_seurat, L_set_relevant_features_all_clusters,
    filtered_interaction_potentials_matrix_all_clusters,
    cytosig_ligands, flag.normalize.non.log, pp_c8)

message(sprintf("[%s] random forest weights...", format(Sys.time(), "%H:%M:%S")))
# the RF wrapper guards on BOTH pairing clusters being present in
# names(significant_regulon_deltas); with scoring restricted to the C8 pairing
# the control cluster is absent — pad an empty entry so the guard passes
sig_for_rf <- significant_regulons_by_cluster
if (!pp_c8$control %in% names(sig_for_rf)) {
  sig_for_rf[[pp_c8$control]] <- data.frame(deltaPagoda = numeric(0),
                                            name = character(0),
                                            class = character(0))
}
decipher_scores_by_regulon_and_cluster <- getRandomForestWeightsAllClustersWParamPairings(
  decipher_seurat, sig_for_rf, regulon_scores_by_cluster,
  interaction_potentials_matrix_clusters_all_clusters,
  L_set_relevant_features_all_clusters, flag.normalize.non.log, pp_c8)

decipher_scores_by_regulon_and_cluster <- lapply(
  decipher_scores_by_regulon_and_cluster, FUN = "listOfDFsRenameColumn",
  original_name = "weighted.spearman.cont", new_name = "decipher_score")
decipher_scores_by_cluster <- lapply(decipher_scores_by_regulon_and_cluster,
                                     FUN = "calculateScoresByCluster")
decipher_scores_by_cluster <- addListNameToDFElements(decipher_scores_by_cluster, "receiver_cluster")

saveRDS(decipher_scores_by_regulon_and_cluster, file.path(out_dir, "decipher_scores_by_regulon_and_cluster.rds"))
saveRDS(decipher_scores_by_cluster, file.path(out_dir, "decipher_scores_by_cluster.rds"))
saveRDS(significant_regulons_by_cluster, file.path(out_dir, "significant_regulons_by_cluster.rds"))
saveRDS(regulon_deltas_by_cluster, file.path(out_dir, "regulon_deltas_by_cluster.rds"))
saveRDS(interaction_deltas_by_cluster, file.path(out_dir, "interaction_deltas_by_cluster.rds"))
saveRDS(interaction_potentials_matrix_clusters_all_clusters, file.path(out_dir, "interaction_potentials_matrix.rds"))

message(sprintf("[%s] DONE %s — top 10 interactions:", format(Sys.time(), "%H:%M:%S"), control_cluster))
print(head(decipher_scores_by_cluster[[C8]], 10))

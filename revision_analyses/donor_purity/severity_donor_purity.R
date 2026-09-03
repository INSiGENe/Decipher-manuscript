# Donor-composition check for Decipher meta-cells — COVID severity dataset (GSE155673).
# Round-2 revision, Reviewer 2 Major comment 3 (severity arm; vaccination done 13 Aug).
#
# Reproduces the meta-cell generation of the SevCOVID_Azimuthl2 run of
#   Decipher-manuscript/scripts/analysis_cellxgene_datasets/6_decipher_pipeline_v1_modularized.R
# (config.json: condition_name="severity_group", case="Severe", control="Healthy",
#  k=1, min_meta_cells=100, minCpc=100, seed=123, NO paramPairings) with the same
# code, while additionally recording which single cells compose each meta-cell.
#
# The instrumented functions are line-for-line copies of generateMetaCellMatrices /
# calculatePseudoBulkCell / calculatePseudoBulkMatrix; RNG-consuming calls
# (sample_n, sample) are identical and in identical order.
#
# Donor ID: sample_id (COV01-style, 12 donors in GSE155673; this object holds the
# Severe [4 donors] + Healthy [5 donors] subset — Moderate lives in MilCOVID).
# NB meta-cell/seed IDs must be matched VERBATIM downstream (never suffix-strip).

repo_path    <- "/Users/edgarbasto/Documents/Decipher-manuscript"
input_rds    <- "/Volumes/MegaEdgar/aws_pull_20260813/results/SevCOVID_Azimuthl2/pre_processing/seurat_object_oi.rds"
archived_rds <- "/Volumes/MegaEdgar/aws_pull_20260813/results/SevCOVID_Azimuthl2/data/pseudobulk_seurat.rds"
output_path  <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/donor_purity/severity"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(tibble)
  library(stringr)
  library(igraph)
  library(magrittr)
  library(dplyr)   # last, so dplyr::filter/select win
})

source(file.path(repo_path, "R/utilities.R"))
source(file.path(repo_path, "R/pre_processing_qc.R"))
source(file.path(repo_path, "R/pre_processing_meta_cells.R"))

# Seurat 5 compatibility: the repo code uses the defunct `slot=` argument of
# GetAssayData; redefine calculateDistVarF with `layer=` — otherwise identical
# (FindVariableFeatures and the layer contents are unchanged, so the distance
# matrix and all RNG draws are unaffected).
calculateDistVarF <- function(seuratObject) {
  if (!inherits(seuratObject, "Seurat")) {
    stop("seuratObject must be a Seurat object")
  }
  seuratObject <- Seurat::FindVariableFeatures(seuratObject)
  variableFeatures <- Seurat::VariableFeatures(seuratObject, assay = "RNA")
  dataMatrix <- Seurat::GetAssayData(seuratObject, assay = "RNA", layer = "data")
  dataMatrix <- as.matrix(dataMatrix)
  CorDist(dataMatrix[variableFeatures, ])
}

############
# Instrumented copies (identical RNG behaviour, plus membership recording) ----
############

calculatePseudoBulkMatrix_instr <- function(rnaCountsMatrix, distanceMatrix, numNearestNeighbors, numMetaCells) {
  subSamples <- sample(colnames(rnaCountsMatrix), numMetaCells)

  if (numNearestNeighbors > 0) {
    knnMatrix <- getNearestNeighbors(distanceMatrix, numNearestNeighbors, subSamples)
  }
  imputedMatrix <- matrix(0, nrow = nrow(rnaCountsMatrix), ncol = numMetaCells)
  rownames(imputedMatrix) <- rownames(rnaCountsMatrix)
  colnames(imputedMatrix) <- subSamples

  membership <- list()
  for (sample_ in subSamples) {
    if (numNearestNeighbors > 0) {
      neighbors <- knnMatrix[sample_, ]
      members <- c(sample_, colnames(rnaCountsMatrix)[neighbors])
      sampleData <- rnaCountsMatrix[, members]
      imputedMatrix[, sample_] <- rowSums(sampleData)
    } else {
      members <- sample_
      imputedMatrix[, sample_] <- rnaCountsMatrix[, sample_]
    }
    membership[[sample_]] <- members
  }
  list(matrix = imputedMatrix, membership = membership)
}

calculatePseudoBulkCell_instr <- function(seuratObject, numNearestNeighbors, numMetaCells) {
  rnaCountsMatrix <- Seurat::GetAssayData(seuratObject, assay = "RNA", layer = "counts")
  rnaCountsMatrix <- as.matrix(rnaCountsMatrix)
  distanceMatrix <- calculateDistVarF(seuratObject = seuratObject)
  calculatePseudoBulkMatrix_instr(
    rnaCountsMatrix = rnaCountsMatrix,
    distanceMatrix = distanceMatrix,
    numNearestNeighbors = numNearestNeighbors,
    numMetaCells = numMetaCells
  )
}

# instrumented copy of generateMetaCellMatrices (the UNPAIRED variant used by
# metaCellModule for this dataset) — loop order and RNG calls identical
generateMetaCellMatrices_instr <- function(seuratObj, paramMinMetaCells = 100, paramMaxMetaCells = 600,
                                           paramMaxScCells, paramK) {
  B_matrices <- list()
  memberships <- list()
  cells_used <- list()

  for (this_cluster in unique(seuratObj$cluster)) {
    cat("Calculating pseudobulk matrices for cluster:", this_cluster, "\n")
    seuratObjectCluster <- seuratObj[, which(seuratObj$cluster == this_cluster), seed = NULL]

    minCellCount <- min(table(seuratObjectCluster$condition))
    minCellCount <- floor(minCellCount / (paramK + 1))

    if (minCellCount < paramMinMetaCells) {
      next
    }

    minCellCount <- min(minCellCount, paramMaxMetaCells)
    cat("Number of pseudobulk cells:", minCellCount, "\n")

    for (this_condition in unique(seuratObjectCluster$condition)) {
      cat("Calculating pseudobulk matrices for condition:", this_condition, "\n")

      conditionData <- seuratObjectCluster@meta.data %>%
        tibble::rownames_to_column(var = "cell") %>%
        dplyr::filter(condition == this_condition)

      if (nrow(conditionData) > paramMaxScCells) {
        cellsToKeep <- conditionData %>% dplyr::select(cell) %>% dplyr::sample_n(size = paramMaxScCells) %>% unlist(use.names = FALSE)
      } else {
        cellsToKeep <- conditionData$cell
      }

      seuratObjectCondition <- seuratObjectCluster[, cellsToKeep, seed = NULL]

      res <- calculatePseudoBulkCell_instr(
        seuratObject = seuratObjectCondition,
        numNearestNeighbors = paramK,
        numMetaCells = minCellCount
      )

      B_matrices[[this_cluster]][[this_condition]] <- res$matrix
      memberships[[this_cluster]][[this_condition]] <- res$membership
      cells_used[[this_cluster]][[this_condition]] <- data.frame(
        out_cluster = this_cluster,
        condition = this_condition,
        n_cells_available = nrow(conditionData),
        n_cells_used = length(cellsToKeep)
      )
    }
  }

  list(matrices = B_matrices, memberships = memberships, cells_used = cells_used)
}

############
# Replicate SevCOVID_Azimuthl2 pipeline up to meta-cell generation ----
############

selected_random_seed <- 123
set.seed(selected_random_seed)

condition_name    <- "severity_group"
case_condition    <- "Severe"
control_condition <- "Healthy"
k_parameter       <- 1
min_meta_cells    <- 100

seurat_oi <- readRDS(input_rds)

seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi, condition_name, case_condition, control_condition)

CpC_data <- generateQCDataByClusterAndCondition(seurat_oi, max(stringr::str_length(unique(seurat_oi$cluster))))
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data, minCpc = 100)

seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed = NULL]

seurat_oi@meta.data$barcode <- rownames(seurat_oi@meta.data)

res <- generateMetaCellMatrices_instr(
  seuratObj = seurat_oi,
  paramMinMetaCells = min_meta_cells,
  paramMaxScCells = 1200 * (k_parameter + 1),
  paramK = k_parameter
)

############
# Verify against archived pseudobulk (the manuscript run's own output) ----
############

archived <- readRDS(archived_rds)
archived_meta <- archived@meta.data
rm(archived)

regen <- bind_rows(lapply(names(res$matrices), function(cl) {
  bind_rows(lapply(names(res$matrices[[cl]]), function(cond) {
    data.frame(cluster = cl, condition = cond,
               seed = colnames(res$matrices[[cl]][[cond]]), stringsAsFactors = FALSE)
  }))
}))
key_regen <- paste(regen$cluster, regen$condition, regen$seed)
key_arch  <- paste(archived_meta$cluster, archived_meta$condition, rownames(archived_meta))
match_summary <- data.frame(
  n_metacells_regenerated = nrow(regen),
  n_metacells_archived = nrow(archived_meta),
  n_key_overlap = sum(key_regen %in% key_arch),
  identical_sets = setequal(key_regen, key_arch)
)
cat("\n=== Archive comparison (verbatim keys) ===\n")
print(match_summary)
write.csv(match_summary, file.path(output_path, "archive_match_summary.csv"), row.names = FALSE)

############
# Donor composition ----
############

donor_lookup <- seurat_oi@meta.data %>%
  select(barcode, sample_id, cluster, condition)

# full barcode -> donor map for downstream use (R2.4 severity LODO folds)
write.csv(donor_lookup %>% select(barcode, sample_id),
          file.path(output_path, "cell_donor_lookup.csv"), row.names = FALSE)

membership_df <- bind_rows(lapply(names(res$memberships), function(cl) {
  bind_rows(lapply(names(res$memberships[[cl]]), function(cond) {
    mem <- res$memberships[[cl]][[cond]]
    bind_rows(lapply(names(mem), function(seed) {
      data.frame(out_cluster = cl, condition = cond, metacell_seed = seed,
                 member_barcode = mem[[seed]], stringsAsFactors = FALSE)
    }))
  }))
}))

membership_df <- membership_df %>%
  left_join(donor_lookup %>% select(barcode, sample_id), by = c("member_barcode" = "barcode")) %>%
  rename(donor = sample_id) %>%
  group_by(out_cluster, condition, metacell_seed) %>%
  mutate(seed_donor = donor[member_barcode == metacell_seed][1]) %>%
  ungroup()

purity_by_metacell <- membership_df %>%
  group_by(out_cluster, condition, metacell_seed, seed_donor) %>%
  summarize(
    n_members = n(),
    n_donors = n_distinct(donor),
    frac_seed_donor = mean(donor == seed_donor),
    frac_modal_donor = max(table(donor)) / n(),
    .groups = "drop"
  )

purity_summary <- purity_by_metacell %>%
  group_by(out_cluster, condition) %>%
  summarize(
    n_metacells = n(),
    mean_frac_seed_donor = mean(frac_seed_donor),
    median_frac_seed_donor = median(frac_seed_donor),
    mean_frac_modal_donor = mean(frac_modal_donor),
    pct_pure = 100 * mean(n_donors == 1),
    pct_majority_seed_donor = 100 * mean(frac_seed_donor > 0.5),
    mean_n_donors_per_metacell = mean(n_donors),
    .groups = "drop"
  )

# donor-blind null: expected frac_seed_donor if the k non-seed members were
# drawn at random from the (subsampled) pool, per cluster x condition
pool_comp <- membership_df %>%
  distinct(out_cluster, condition, member_barcode, donor) %>%
  count(out_cluster, condition, donor, name = "n_pool_cells")

# null: for a meta-cell seeded by donor d, expected fraction of members from d
# = (1 + k * p_d) / (k + 1); average over seeds weighted by donor pool shares
null_expectation <- pool_comp %>%
  group_by(out_cluster, condition) %>%
  mutate(p_donor = n_pool_cells / sum(n_pool_cells)) %>%
  summarize(null_frac_seed_donor = sum(p_donor * (1 + k_parameter * p_donor) / (k_parameter + 1)),
            .groups = "drop")

observed_vs_null <- purity_summary %>%
  select(out_cluster, condition, mean_frac_seed_donor) %>%
  left_join(null_expectation, by = c("out_cluster", "condition")) %>%
  mutate(enrichment = mean_frac_seed_donor / null_frac_seed_donor)

# composition skew: donor cell counts per arm + effective donor n (inverse Simpson)
donor_cells_by_arm <- donor_lookup %>%
  count(cluster, condition, sample_id, name = "n_cells") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0)

composition_skew <- donor_lookup %>%
  count(cluster, condition, sample_id) %>%
  group_by(cluster, condition) %>%
  summarize(
    n_donors_present = n(),
    n_cells = sum(n),
    max_donor_share = max(n) / sum(n),
    effective_n_donors = 1 / sum((n / sum(n))^2),
    .groups = "drop"
  )

############
# Counts for supplementary table (R2.3 / R2.6) ----
############

cells_used_df <- bind_rows(lapply(res$cells_used, bind_rows))

counts_table <- donor_lookup %>%
  group_by(cluster, condition) %>%
  summarize(n_donors = n_distinct(sample_id), n_cells = n(), .groups = "drop") %>%
  left_join(cells_used_df %>% select(out_cluster, condition, n_cells_used),
            by = c("cluster" = "out_cluster", "condition" = "condition")) %>%
  left_join(purity_summary %>% select(out_cluster, condition, n_metacells),
            by = c("cluster" = "out_cluster", "condition" = "condition"))

metacell_donor_diversity <- membership_df %>%
  group_by(out_cluster, condition) %>%
  summarize(n_donors_represented_in_metacells = n_distinct(donor), .groups = "drop")

counts_table <- counts_table %>%
  left_join(metacell_donor_diversity, by = c("cluster" = "out_cluster", "condition" = "condition"))

############
# Outputs ----
############

write.csv(membership_df, file.path(output_path, "metacell_membership.csv"), row.names = FALSE)
write.csv(purity_by_metacell, file.path(output_path, "donor_purity_by_metacell.csv"), row.names = FALSE)
write.csv(purity_summary, file.path(output_path, "donor_purity_summary.csv"), row.names = FALSE)
write.csv(observed_vs_null, file.path(output_path, "observed_vs_null_purity.csv"), row.names = FALSE)
write.csv(donor_cells_by_arm, file.path(output_path, "donor_cells_by_arm.csv"), row.names = FALSE)
write.csv(composition_skew, file.path(output_path, "donor_composition_skew.csv"), row.names = FALSE)
write.csv(counts_table, file.path(output_path, "counts_donors_cells_metacells.csv"), row.names = FALSE)

cat("\n=== Donor purity summary (per cluster x condition) ===\n")
print(as.data.frame(purity_summary), digits = 3)
cat("\n=== Observed vs null purity ===\n")
print(as.data.frame(observed_vs_null), digits = 3)
cat("\n=== Composition skew ===\n")
print(as.data.frame(composition_skew), digits = 3)
cat("\n=== Counts table ===\n")
print(as.data.frame(counts_table))
cat("\n=== Overall ===\n")
cat("Total donors in dataset:", n_distinct(donor_lookup$sample_id), "\n")
cat("Overall mean fraction of members from seed donor:", round(mean(purity_by_metacell$frac_seed_donor), 3), "\n")
cat("Overall % fully donor-pure meta-cells:", round(100 * mean(purity_by_metacell$n_donors == 1), 1), "%\n")
cat("\nDone. Outputs in:", output_path, "\n")

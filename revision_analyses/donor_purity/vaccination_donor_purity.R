# Donor-composition check for Decipher meta-cells — vaccination dataset (GSE171964)
# Round-2 revision, Reviewer 2 Major comment 3.
#
# Reproduces the meta-cell generation of
#   Decipher-manuscript/scripts/analysis_specific_datasets/covid/covid_2_decipher_pipeline_v1_modularized.R
# with the same code (sourced from the Decipher-manuscript repo) and the same seed,
# while additionally recording which single cells compose each meta-cell.
# It then joins donor IDs (pt_id) and reports donor purity per meta-cell, plus the
# donors/cells/meta-cells counts needed for the supplementary summary table.
#
# The instrumented functions below are line-for-line copies of the originals; the
# only additions are membership bookkeeping. RNG-consuming calls (sample_n, sample)
# are identical and in identical order, so meta-cells match the pipeline exactly.

repo_path    <- "/Users/edgarbasto/Documents/Decipher-manuscript"
scratch      <- "/private/tmp/claude-501/-Users-edgarbasto-Documents-Decipher-revision-2/5cdd2b80-1575-43a1-94a9-f02fa349a71e/scratchpad/covid_donor_check"
input_rds    <- file.path(scratch, "covid/pre_processing/seurat_object_oi.rds")
archived_rds <- file.path(scratch, "covid/data/pseudobulk_seurat.rds")
output_path  <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/donor_purity"
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
  rnaCountsMatrix <- Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "counts")
  rnaCountsMatrix <- as.matrix(rnaCountsMatrix)
  distanceMatrix <- calculateDistVarF(seuratObject = seuratObject)
  calculatePseudoBulkMatrix_instr(
    rnaCountsMatrix = rnaCountsMatrix,
    distanceMatrix = distanceMatrix,
    numNearestNeighbors = numNearestNeighbors,
    numMetaCells = numMetaCells
  )
}

generateMetaCellMatricesWPairings_instr <- function(seuratObj, paramMinMetaCells = 100, paramMaxMetaCells = 600,
                                                    paramMaxScCells, paramK, paramPairings) {
  B_matrices <- list()
  memberships <- list()
  cells_used <- list()

  for (this_row in c(1:nrow(paramPairings))) {
    cat("Calculating pseudobulk matrices for cluster:", paramPairings$case[this_row], "\n")
    case_cluster <- paramPairings$case[this_row]
    control_cluster <- paramPairings$control[this_row]
    case_cells <- seuratObj@meta.data %>% filter(cluster %in% case_cluster & condition == "case") %>% pull(barcode)
    control_cells <- seuratObj@meta.data %>% filter(cluster %in% control_cluster & condition == "control") %>% pull(barcode)
    all_cells <- c(case_cells, control_cells)

    seuratObjectCluster <- seuratObj[, all_cells, seed = NULL]
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
        filter(condition == this_condition)

      if (nrow(conditionData) > paramMaxScCells) {
        cellsToKeep <- conditionData %>% select(cell) %>% sample_n(size = paramMaxScCells) %>% unlist(use.names = FALSE)
      } else {
        cellsToKeep <- conditionData$cell
      }

      seuratObjectCondition <- seuratObjectCluster[, cellsToKeep, seed = NULL]

      res <- calculatePseudoBulkCell_instr(
        seuratObject = seuratObjectCondition,
        numNearestNeighbors = paramK,
        numMetaCells = minCellCount
      )

      this_out_cluster <- unique(seuratObjectCondition$cluster)
      B_matrices[[this_out_cluster]][[this_condition]] <- res$matrix
      memberships[[this_out_cluster]][[this_condition]] <- res$membership
      cells_used[[this_out_cluster]][[this_condition]] <- data.frame(
        pairing_case = case_cluster,
        pairing_control = control_cluster,
        out_cluster = this_out_cluster,
        condition = this_condition,
        n_cells_available = nrow(conditionData),
        n_cells_used = length(cellsToKeep)
      )
    }
  }

  list(matrices = B_matrices, memberships = memberships, cells_used = cells_used)
}

############
# Replicate covid_2 pipeline up to meta-cell generation ----
############

selected_random_seed <- 123
set.seed(selected_random_seed)

min_cells_per_cluster_condition <- 100
condition_name <- "condition"
case_condition <- 22
control_condition <- 0
k_parameter <- 2

seurat_oi <- readRDS(input_rds)

seurat_oi$orig.condition <- seurat_oi[[condition_name]]
seurat_oi <- mapConditionsInSeurat(seurat_oi, condition_name, case_condition, control_condition)

CpC_data <- generateQCDataByClusterAndCondition(seurat_oi, max(stringr::str_length(unique(seurat_oi$cluster))))
clusters_passing_CpC_filter <- getClustersPassingCpCFilter(CpC_data, minCpc = 100)
clusters_passing_CpC_filter <- c(clusters_passing_CpC_filter, "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells")

seurat_oi <- seurat_oi[, which(seurat_oi$cluster %in% clusters_passing_CpC_filter), seed = NULL]

paramPairings <- data.frame(
  case = clusters_passing_CpC_filter,
  control = clusters_passing_CpC_filter
)
paramPairings <- paramPairings %>%
  mutate(control = if_else(case == "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells", "CD14_plus_monocytes", control))

min_counts <- seurat_oi@meta.data %>%
  group_by(cluster) %>%
  count(condition) %>%
  ungroup()

groups <- createGroupsFromPairings(paramPairings)
paramPairings_min_n <- calculateMinimumN(groups, min_counts, paramPairings)

seurat_oi@meta.data$barcode <- rownames(seurat_oi@meta.data)

res <- generateMetaCellMatricesWPairings_instr(
  seuratObj = seurat_oi,
  paramMaxScCells = 1200 * (k_parameter + 1),
  paramK = k_parameter,
  paramPairings = paramPairings_min_n
)

############
# Verify against archived pseudobulk (Jul 2024 vintage) ----
############

archived <- readRDS(archived_rds)
archived_meta <- archived@meta.data
rm(archived)

regen_seeds <- unlist(lapply(names(res$matrices), function(cl) {
  lapply(names(res$matrices[[cl]]), function(cond) colnames(res$matrices[[cl]][[cond]]))
}))
match_summary <- data.frame(
  n_metacells_regenerated = length(regen_seeds),
  n_metacells_archived = nrow(archived_meta),
  n_seed_overlap = length(intersect(regen_seeds, archived_meta$cell)),
  identical_sets = setequal(regen_seeds, archived_meta$cell)
)
cat("\n=== Archive comparison ===\n")
print(match_summary)
write.csv(match_summary, file.path(output_path, "archive_match_summary.csv"), row.names = FALSE)

############
# Donor composition ----
############

donor_lookup <- seurat_oi@meta.data %>%
  select(barcode, pt_id, day, cluster, condition)

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
  left_join(donor_lookup %>% select(barcode, pt_id), by = c("member_barcode" = "barcode")) %>%
  group_by(out_cluster, condition, metacell_seed) %>%
  mutate(seed_donor = pt_id[member_barcode == metacell_seed][1]) %>%
  ungroup()

purity_by_metacell <- membership_df %>%
  group_by(out_cluster, condition, metacell_seed, seed_donor) %>%
  summarize(
    n_members = n(),
    n_donors = n_distinct(pt_id),
    frac_seed_donor = mean(pt_id == seed_donor),
    frac_modal_donor = max(table(pt_id)) / n(),
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

############
# Counts for supplementary table (R2.3 / Minor 1) ----
############

cells_used_df <- bind_rows(lapply(res$cells_used, bind_rows))

counts_table <- donor_lookup %>%
  group_by(cluster, condition) %>%
  summarize(n_donors = n_distinct(pt_id), n_cells = n(), .groups = "drop") %>%
  left_join(cells_used_df %>% select(out_cluster, condition, n_cells_used),
            by = c("cluster" = "out_cluster", "condition" = "condition")) %>%
  left_join(purity_summary %>% select(out_cluster, condition, n_metacells),
            by = c("cluster" = "out_cluster", "condition" = "condition"))

metacell_donor_diversity <- membership_df %>%
  group_by(out_cluster, condition) %>%
  summarize(n_donors_represented_in_metacells = n_distinct(pt_id), .groups = "drop")

counts_table <- counts_table %>%
  left_join(metacell_donor_diversity, by = c("cluster" = "out_cluster", "condition" = "condition"))

############
# Outputs ----
############

write.csv(membership_df, file.path(output_path, "metacell_membership.csv"), row.names = FALSE)
write.csv(purity_by_metacell, file.path(output_path, "donor_purity_by_metacell.csv"), row.names = FALSE)
write.csv(purity_summary, file.path(output_path, "donor_purity_summary.csv"), row.names = FALSE)
write.csv(counts_table, file.path(output_path, "counts_donors_cells_metacells.csv"), row.names = FALSE)

cat("\n=== Donor purity summary (per cluster x condition) ===\n")
print(as.data.frame(purity_summary), digits = 3)
cat("\n=== Counts table ===\n")
print(as.data.frame(counts_table))
cat("\n=== Overall ===\n")
cat("Total donors in dataset:", n_distinct(donor_lookup$pt_id), "\n")
cat("Overall mean fraction of members from seed donor:", round(mean(purity_by_metacell$frac_seed_donor), 3), "\n")
cat("Overall % fully donor-pure meta-cells:", round(100 * mean(purity_by_metacell$n_donors == 1), 1), "%\n")
cat("\nDone. Outputs in:", output_path, "\n")

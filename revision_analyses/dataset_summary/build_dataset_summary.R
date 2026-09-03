#!/usr/bin/env Rscript
# build_dataset_summary.R
# Reviewer-requested per-dataset summary table (round 2):
#   per comparison x cell type (analysed clusters only) x condition:
#   donors, cells, meta-cells, candidate LR pairs, retained TFs, predicted interactions, k.
#
# Rerunnable. Stage 1 caches the (small) metadata of the large Seurat objects under
# dataset_summary/cache/ so stage 2 (table building, donor-column choice) can be iterated
# without re-reading multi-GB RDS files. Delete cache/ to force a full reload.
#
# Usage (from analysis_r2/):
#   caffeinate -i Rscript dataset_summary/build_dataset_summary.R            # all comparisons
#   caffeinate -i Rscript dataset_summary/build_dataset_summary.R ERP covid  # subset (stage 1 only for these)
#
# Environment notes (HANDOVER.md): R 4.6 / Seurat 5 reading Seurat v3/v4 objects ->
# set the v3 assay option before loading and access metadata via obj@meta.data only.

options(Seurat.object.assay.version = "v3")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

ROOT      <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(ROOT) || !dir.exists(file.path(ROOT, "aws_data"))) ROOT <- getwd()
RESULTS   <- file.path(ROOT, "aws_data", "results")
OUT_DIR   <- file.path(ROOT, "dataset_summary")
CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

ALL_COMPARISONS <- c("5yr_pic", "BCG", "cord_pic", "covid", "cz_cf_bronchial_biopsy",
                     "cz_hnscc_hpv", "cz_hpap_t1d_islets", "cz_human_kidney_v1.5",
                     "cz_influenza", "ERP", "lupus", "MilCOVID_Azimuthl2", "sepsis",
                     "SevCOVID_Azimuthl2", "TNBC")

# ---------------------------------------------------------------------------
# Per-dataset configuration: human-readable label, donor column in
# pre_processing/seurat_object_oi.rds@meta.data (NA if none identifiable).
# Filled in after inspecting the cached metadata (see NOTES.md).
# ---------------------------------------------------------------------------
DATASET_CONFIG <- list(
  "5yr_pic"                = list(label = "5-year-old PBMC, poly(I:C) stimulated vs unstimulated",              donor_col = "Sample"),
  "BCG"                    = list(label = "BCG (GSE232186), day 21 vs day 0",                                    donor_col = NA),
  "cord_pic"               = list(label = "Cord blood mononuclear cells, poly(I:C) stimulated vs unstimulated", donor_col = "Sample"),
  "covid"                  = list(label = "COVID-19 mRNA vaccination, day 22 vs day 0 (Arunachalam 2021)",      donor_col = "pt_id"),
  "cz_cf_bronchial_biopsy" = list(label = "Cystic fibrosis bronchial biopsy vs normal airway (CZ CELLxGENE)",   donor_col = "donor_id"),
  "cz_hnscc_hpv"           = list(label = "Oropharyngeal HNSCC tumour vs adjacent normal (CZ CELLxGENE)",       donor_col = "donor_id"),
  "cz_hpap_t1d_islets"     = list(label = "Type 1 diabetes pancreatic islets vs non-diabetic (HPAP, CZ CELLxGENE)", donor_col = "donor_id"),
  "cz_human_kidney_v1.5"   = list(label = "Chronic kidney disease vs reference kidney (KPMP atlas v1.5, CZ CELLxGENE)", donor_col = "donor_id"),
  "cz_influenza"           = list(label = "Severe influenza PBMC vs healthy (CZ CELLxGENE)",                    donor_col = "donor_id"),
  "ERP"                    = list(label = "ER+ breast cancer, anti-PD1 T-cell expansion E vs NE (Bassez 2021)", donor_col = "patient_id"),
  "lupus"                  = list(label = "Systemic lupus erythematosus PBMC vs healthy (CZ CELLxGENE)",        donor_col = "donor_id"),
  "MilCOVID_Azimuthl2"     = list(label = "Moderate COVID-19 PBMC vs healthy (Arunachalam 2020; Azimuth L2)",   donor_col = "sample_id"),
  "sepsis"                 = list(label = "Sepsis: ICU-SEP vs ICU-NoSEP (Reyes 2020)",                          donor_col = "donor_id"),
  "SevCOVID_Azimuthl2"     = list(label = "Severe COVID-19 PBMC vs healthy (Arunachalam 2020; Azimuth L2)",     donor_col = "sample_id"),
  "TNBC"                   = list(label = "TNBC, anti-PD1 T-cell expansion E vs NE (Bassez 2021)",              donor_col = "patient_id")
)

# ---------------------------------------------------------------------------
# Stage 1: cache metadata + small count objects for one comparison
# ---------------------------------------------------------------------------
cache_comparison <- function(comp) {
  cache_file <- file.path(CACHE_DIR, paste0(comp, ".rds"))
  if (file.exists(cache_file)) { message("[", comp, "] cache exists, skipping load"); return(invisible(cache_file)) }
  message("[", comp, "] loading ... ", format(Sys.time(), "%H:%M:%S"))
  d  <- file.path(RESULTS, comp, "data")
  pp <- file.path(RESULTS, comp, "pre_processing")

  param <- read.csv(file.path(d, "parameter_record.csv"), row.names = 1)

  pb <- readRDS(file.path(d, "pseudobulk_seurat.rds"))
  pb_meta <- pb@meta.data
  pb_meta$metacell_id <- rownames(pb_meta)
  rm(pb); gc(verbose = FALSE)

  so <- readRDS(file.path(pp, "seurat_object_oi.rds"))
  sc_meta <- so@meta.data
  sc_meta$cell_id <- rownames(sc_meta)
  rm(so); gc(verbose = FALSE)

  scores  <- readRDS(file.path(d, "decipher_scores_by_cluster.rds"))
  regs    <- readRDS(file.path(d, "significant_regulons_by_cluster.rds"))
  lrel    <- readRDS(file.path(d, "L_set_relevant_features_all_clusters.rds"))
  ideltas <- readRDS(file.path(d, "interaction_deltas_by_cluster.rds"))
  lrmk    <- readRDS(file.path(d, "lr_markers_by_cluster.rds"))
  lset    <- readRDS(file.path(d, "L_set.rds"))

  count_rows <- function(x) if (is.null(x)) NA_integer_ else nrow(x)
  counts <- data.frame(
    cluster = union(names(scores), names(lrel)),
    stringsAsFactors = FALSE
  )
  counts$n_candidate_lr_pairs     <- sapply(counts$cluster, function(cl) if (cl %in% names(lrel)) length(unique(lrel[[cl]]$interaction)) else NA_integer_)
  counts$n_lr_pairs_diff_potential <- sapply(counts$cluster, function(cl) count_rows(ideltas[[cl]]))
  counts$n_lr_genes_de            <- sapply(counts$cluster, function(cl) count_rows(lrmk[[cl]]))
  counts$n_retained_tfs           <- sapply(counts$cluster, function(cl) if (cl %in% names(regs)) length(unique(regs[[cl]]$name)) else NA_integer_)
  counts$n_predicted_interactions <- sapply(counts$cluster, function(cl) count_rows(scores[[cl]]))
  counts$n_unique_predicted_lr    <- sapply(counts$cluster, function(cl) if (cl %in% names(scores)) length(unique(scores[[cl]]$interaction)) else NA_integer_)
  counts$in_scores <- counts$cluster %in% names(scores)

  saveRDS(list(comparison = comp, param = param, pb_meta = pb_meta, sc_meta = sc_meta,
               counts = counts, n_lset_pairs = nrow(lset), scores_clusters = names(scores),
               cached_at = Sys.time()),
          cache_file)
  message("[", comp, "] cached  ", format(Sys.time(), "%H:%M:%S"),
          "  cells=", nrow(sc_meta), " metacells=", nrow(pb_meta), " clusters_scored=", length(scores))
  rm(pb_meta, sc_meta, scores, regs, lrel, ideltas, lrmk, lset); gc(verbose = FALSE)
  invisible(cache_file)
}

# ---------------------------------------------------------------------------
# Stage 2: build per-cluster x condition table from the cache
# ---------------------------------------------------------------------------
# Identify the case/control condition column in the single-cell metadata: the pipeline
# relabels to "case"/"control"; some datasets keep the native labels in `condition` and
# the pipeline labels in `condition_original` (or vice versa).
pick_condition_col <- function(meta) {
  cands <- intersect(c("condition", "condition_original", "original_condition", "Condition"), colnames(meta))
  for (cc in cands) {
    v <- unique(as.character(meta[[cc]]))
    if (all(v %in% c("case", "control"))) return(cc)
  }
  NA_character_
}

summarise_comparison <- function(comp) {
  cf <- readRDS(file.path(CACHE_DIR, paste0(comp, ".rds")))
  cfg <- DATASET_CONFIG[[comp]]
  sc <- cf$sc_meta; pb <- cf$pb_meta
  cond_col <- pick_condition_col(sc)
  if (is.na(cond_col)) stop("[", comp, "] no case/control condition column in single-cell metadata")
  sc$.cond <- as.character(sc[[cond_col]])
  pb$.cond <- as.character(pb$condition)
  sc$.cluster <- as.character(sc$cluster)
  pb$.cluster <- as.character(pb$cluster)

  donor_col <- cfg$donor_col
  if (!is.na(donor_col) && !donor_col %in% colnames(sc)) stop("[", comp, "] donor column ", donor_col, " not in metadata")

  analysed <- cf$scores_clusters
  conds <- c("case", "control")

  rows <- list()
  for (cl in analysed) {
    cnt <- cf$counts[cf$counts$cluster == cl, ]
    for (cd in conds) {
      sc_sub <- sc[sc$.cluster == cl & sc$.cond == cd, , drop = FALSE]
      pb_sub <- pb[pb$.cluster == cl & pb$.cond == cd, , drop = FALSE]
      n_don <- if (is.na(donor_col)) NA_integer_ else length(unique(na.omit(as.character(sc_sub[[donor_col]]))))
      rows[[length(rows) + 1]] <- data.frame(
        comparison = comp,
        dataset_label = cfg$label,
        cell_type = cl,
        condition = cd,
        donor_column = ifelse(is.na(donor_col), NA_character_, donor_col),
        n_donors = n_don,
        n_cells = nrow(sc_sub),
        # generateMetaCellMatrices() subsamples each arm to at most 1200*(k+1) cells
        n_cells_input_to_metacells = min(nrow(sc_sub), 1200 * (cf$param$k + 1)),
        n_metacells = nrow(pb_sub),
        n_candidate_lr_pairs = cnt$n_candidate_lr_pairs,
        n_lr_pairs_diff_potential = cnt$n_lr_pairs_diff_potential,
        n_retained_tfs = cnt$n_retained_tfs,
        n_predicted_interactions = cnt$n_predicted_interactions,
        n_unique_predicted_lr = cnt$n_unique_predicted_lr,
        k = cf$param$k,
        min_meta_cells = cf$param$min_meta_cells,
        condition_column_sc = cond_col,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- bind_rows(rows)

  # consistency check of the meta-cell rule: per cluster, expected meta-cells per arm
  # = min(600, floor(min(arm cells)/(k+1))); flag deviations
  out <- out %>% group_by(cell_type) %>%
    mutate(expected_metacells_per_arm = min(600L, floor(min(n_cells) / (k + 1))),
           metacell_rule_ok = all(n_metacells == expected_metacells_per_arm)) %>%
    ungroup() %>% as.data.frame()

  # comparison-level extras
  n_don_total <- if (is.na(donor_col)) NA_integer_ else length(unique(na.omit(as.character(sc[[donor_col]]))))
  n_don_by_cond <- if (is.na(donor_col)) c(case = NA_integer_, control = NA_integer_) else
    sapply(conds, function(cd) length(unique(na.omit(as.character(sc[[donor_col]][sc$.cond == cd])))))
  extra <- data.frame(
    comparison = comp,
    dataset_label = cfg$label,
    n_cell_types_in_sc_object = length(unique(sc$.cluster)),
    n_cell_types_in_metacell_object = length(unique(pb$.cluster)),
    n_cell_types_analysed = length(analysed),
    n_cells_sc_object_total = nrow(sc),
    n_cells_case_total = sum(sc$.cond == "case"),
    n_cells_control_total = sum(sc$.cond == "control"),
    n_cells_analysed_clusters = sum(sc$.cluster %in% analysed),
    n_metacells_total = nrow(pb),
    n_metacells_analysed_clusters = sum(pb$.cluster %in% analysed),
    donor_column = ifelse(is.na(donor_col), NA_character_, donor_col),
    n_donors_total = n_don_total,
    n_donors_case = unname(n_don_by_cond["case"]),
    n_donors_control = unname(n_don_by_cond["control"]),
    n_lset_pairs_database = cf$n_lset_pairs,
    k = cf$param$k,
    stringsAsFactors = FALSE
  )
  list(rows = out, extra = extra, cond_col = cond_col,
       sc_clusters = sort(unique(sc$.cluster)), pb_clusters = sort(unique(pb$.cluster)))
}

# ---------------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  comps <- if (length(args)) args else ALL_COMPARISONS
  failures <- character()
  for (comp in comps) {
    ok <- tryCatch({ cache_comparison(comp); TRUE },
                   error = function(e) { message("[", comp, "] FAILED stage 1: ", conditionMessage(e)); FALSE })
    if (!ok) failures <- c(failures, comp)
  }

  cached <- ALL_COMPARISONS[file.exists(file.path(CACHE_DIR, paste0(ALL_COMPARISONS, ".rds")))]
  per_row <- list(); per_comp <- list()
  for (comp in cached) {
    res <- tryCatch(summarise_comparison(comp),
                    error = function(e) { message("[", comp, "] FAILED stage 2: ", conditionMessage(e)); failures <<- c(failures, comp); NULL })
    if (is.null(res)) next
    per_row[[comp]] <- res$rows
    per_comp[[comp]] <- res$extra
  }
  if (!length(per_row)) { message("nothing summarised"); return(invisible()) }
  tab <- bind_rows(per_row)
  write.csv(tab, file.path(OUT_DIR, "ST_dataset_summary.csv"), row.names = FALSE)

  roll_counts <- tab %>%
    group_by(comparison) %>%
    summarise(
      n_cells_case = sum(n_cells[condition == "case"]),
      n_cells_control = sum(n_cells[condition == "control"]),
      n_metacells_case = sum(n_metacells[condition == "case"]),
      n_metacells_control = sum(n_metacells[condition == "control"]),
      total_candidate_lr_pairs = sum(n_candidate_lr_pairs[condition == "case"]),
      total_lr_pairs_diff_potential = sum(n_lr_pairs_diff_potential[condition == "case"]),
      total_retained_tfs = sum(n_retained_tfs[condition == "case"]),
      total_predicted_interactions = sum(n_predicted_interactions[condition == "case"]),
      total_unique_predicted_lr = sum(n_unique_predicted_lr[condition == "case"]),
      .groups = "drop"
    )
  roll <- bind_rows(per_comp) %>% left_join(roll_counts, by = "comparison") %>%
    select(comparison, dataset_label, k, n_cell_types_analysed, n_cell_types_in_sc_object,
           donor_column, n_donors_total, n_donors_case, n_donors_control,
           n_cells_analysed_clusters, n_cells_case, n_cells_control, n_cells_sc_object_total,
           n_metacells_analysed_clusters, n_metacells_case, n_metacells_control, n_metacells_total,
           n_lset_pairs_database, total_candidate_lr_pairs, total_lr_pairs_diff_potential,
           total_retained_tfs, total_predicted_interactions, total_unique_predicted_lr)
  write.csv(roll, file.path(OUT_DIR, "ST_dataset_summary_by_comparison.csv"), row.names = FALSE)
  message("wrote ", nrow(tab), " rows / ", nrow(roll), " comparisons")
  if (length(failures)) message("FAILED: ", paste(unique(failures), collapse = ", "))
  invisible(list(tab = tab, roll = roll))
}

if (sys.nframe() == 0L) main()

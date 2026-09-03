# R2.2 ablations on the manuscript's CytoSig AUROC benchmark.
#
# For each dataset, builds ablated Decipher score variants and evaluates them
# with the manuscript's own benchmark machinery (getPredictionsResponsesForMethods
# + plotROCAndExtractAUC from Decipher-manuscript R/vis_data_wrangling.R), so the
# comparison to the published Fig 2e Decipher AUROC is exact and like-for-like.
#
# Variants (reviewer's ablation menu, R2.2):
#   Decipher      — the archive's decipher_scores_by_cluster.rds (full framework:
#                   RF permutation importance x sign(spearman) x deltaPagoda,
#                   mean over significant TFs). Reference.
#   RidgeRegr     — replace RF permutation importance with ridge |beta_j|*sd(x_j)
#                   (cv.glmnet alpha=0, lambda.min); same sign convention, same
#                   TF aggregation. "linear/regularized regression instead of RF".
#   SpearmanCorr  — mean over TFs of possible.spearman.cont = rho(LR,TF) x
#                   deltaPagoda (stored by the pipeline; no model at all).
#                   "simple LR-TF correlation instead of RF regression".
#   LROnly        — differential LR signalling potential alone (the pipeline's
#                   interaction_deltas avg_log2FC), no TF layer at all.
#                   "LR potential alone without downstream TF info".
#
# Usage: Rscript run_ablation_auroc.R [n_cores] [dataset_key ...]
#        (no keys = all 15 benchmark datasets)
# Output: results/ablation_auroc.csv (dataset x method x threshold), resumable
#         per dataset; ROC PNGs under figures/<ds>/.

repo_path  <- "/Users/edgarbasto/Documents/Decipher-manuscript"
drive_root <- "/Volumes/MegaEdgar/aws_pull_20260813/results"
local_root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/aws_data"
out_dir    <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/ablations"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "results", "ablation_auroc.csv")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(magrittr); library(data.table)
  library(pROC); library(ggplot2); library(glmnet); library(parallel)
  library(Seurat)
})

source(file.path(repo_path, "R/vis_data_wrangling.R"))
source(file.path(repo_path, "R/intercellular_tailoring.R"))

# dataset key -> results/ folder name (from 1_load_all_results.r)
DATASETS <- c(
  "5yr_pic" = "5yr_pic", "bcg" = "BCG", "cord_pic" = "cord_pic",
  "covid" = "covid", "erp" = "ERP", "lupus" = "lupus", "sepsis" = "sepsis",
  "tnbc" = "TNBC", "cz_influenza" = "cz_influenza",
  "cz_hpap_t1d_islets" = "cz_hpap_t1d_islets", "cz_hnscc_hpv" = "cz_hnscc_hpv",
  "cz_human_kidney_v1.5" = "cz_human_kidney_v1.5",
  "cz_cf_bronchial_biopsy" = "cz_cf_bronchial_biopsy",
  "SevCOVID_Azimuthl2" = "SevCOVID_Azimuthl2",
  "MilCOVID_Azimuthl2" = "MilCOVID_Azimuthl2"
)

args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
keys <- if (length(args) >= 2) args[-1] else names(DATASETS)
stopifnot(all(keys %in% names(DATASETS)))

# prefer the local aws_data copy of the data store, fall back to the drive
data_dir_for <- function(folder) {
  local <- file.path(local_root, paste0(folder, "_data"))
  if (dir.exists(local)) local else file.path(drive_root, folder, "data")
}

mapping_table <- read.csv(file.path(repo_path, "reference_data", "cytosig_mapping_table_ligands_genes.csv"), header = TRUE)
L.set <- getForrestLRDatabase(file.path(repo_path, "reference_data", "connectomedb_forrest_lrc2p.csv")) %>%
  mutate(interaction = paste(ligand, receptor, sep = "-"), lr = interaction) %>%
  unique()

as_score_df <- function(df, cl) {
  df %>% mutate(sender_cluster = "mixed", receiver_cluster = cl) %>%
    select(interaction, score, ligand, receptor, sender_cluster, receiver_cluster)
}

done <- if (file.exists(out_csv)) unique(read.csv(out_csv)$dataset) else character(0)

for (key in keys) {
  if (key %in% done) { message(key, ": already done"); next }
  folder <- DATASETS[[key]]
  ddir <- data_dir_for(folder)
  zdir <- file.path(drive_root, folder, "cytosig", "0_outputs", "z_score")
  if (!dir.exists(zdir)) { message("SKIP ", key, ": no cytosig z_score dir"); next }
  message(sprintf("[%s] %s (data: %s)", format(Sys.time(), "%H:%M:%S"), key, ddir))

  dec_by_cluster <- readRDS(file.path(ddir, "decipher_scores_by_cluster.rds"))
  dec_by_reg    <- readRDS(file.path(ddir, "decipher_scores_by_regulon_and_cluster.rds"))
  deltas        <- readRDS(file.path(ddir, "interaction_deltas_by_cluster.rds"))
  potentials    <- readRDS(file.path(ddir, "interaction_potentials_matrix_clusters_all_clusters.rds"))
  regulons      <- readRDS(file.path(ddir, "regulon_scores_by_cluster.rds"))
  sig_regulons  <- readRDS(file.path(ddir, "significant_regulons_by_cluster.rds"))
  seurat_object_oi <- readRDS(file.path(ddir, "pseudobulk_seurat.rds"))

  clusters <- names(dec_by_cluster)

  # interaction -> ligand/receptor map (from the pipeline's own tables; robust
  # to hyphens inside gene names, which make string-splitting ambiguous)
  lig_map <- bind_rows(lapply(dec_by_cluster, function(d)
    d %>% select(interaction, ligand, receptor))) %>% distinct()

  # --- Decipher (reference) ---
  v_decipher <- setNames(lapply(clusters, function(cl) {
    dec_by_cluster[[cl]] %>% rename(score = decipher_score) %>% as_score_df(cl)
  }), clusters)

  # --- SpearmanCorr ---
  v_spearman <- setNames(lapply(clusters, function(cl) {
    dec_by_reg[[cl]] %>%
      group_by(interaction) %>%
      summarize(score = mean(possible.spearman.cont), .groups = "drop") %>%
      left_join(lig_map, by = "interaction") %>% as_score_df(cl)
  }), clusters)

  # --- LROnly ---
  v_lronly <- setNames(lapply(clusters, function(cl) {
    deltas[[cl]] %>%
      transmute(interaction = name, score = avg_log2FC) %>%
      inner_join(lig_map, by = "interaction") %>% as_score_df(cl)
  }), clusters)

  # --- RidgeRegr (refits; parallel over clusters) ---
  ridge_list <- mclapply(clusters, function(cl) {
    pot <- potentials[[cl]]
    x <- t(pot)
    sds <- apply(x, 2, sd)
    edge_info <- dec_by_reg[[cl]] %>%
      select(interaction, regulon, spearman.cor, regulon.val) %>% distinct()
    tfs <- sig_regulons[[cl]]$name
    contribs <- lapply(tfs, function(tf) {
      y <- as.numeric(regulons[[cl]][tf, colnames(pot)])
      set.seed(123)
      fit <- tryCatch(cv.glmnet(x, y, alpha = 0, nfolds = 5), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      beta <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      data.frame(interaction = rownames(pot), regulon = tf,
                 imp = abs(beta) * sds, stringsAsFactors = FALSE)
    })
    bind_rows(contribs) %>%
      inner_join(edge_info, by = c("interaction", "regulon")) %>%
      mutate(contrib = imp * sign(spearman.cor) * regulon.val) %>%
      group_by(interaction) %>%
      summarize(score = mean(contrib), .groups = "drop") %>%
      left_join(lig_map, by = "interaction") %>% as_score_df(cl)
  }, mc.cores = n_cores)
  v_ridge <- setNames(ridge_list, clusters)

  # --- CytoSig responses + AUROC via the manuscript's own functions ---
  z_score_folder <- paste0(zdir, "/")
  cytosig_sig <- summarizeZScores(list.files(z_score_folder), z_score_folder, mapping_table)

  figpath <- file.path(out_dir, "figures", key)
  dir.create(figpath, recursive = TRUE, showWarnings = FALSE)

  results_to_compare <- list(Decipher = v_decipher, RidgeRegr = v_ridge,
                             SpearmanCorr = v_spearman, LROnly = v_lronly)

  pr <- getPredictionsResponsesForMethods(results_to_compare, cytosig_sig,
                                          L.set = L.set, seurat_object_oi, figpath)
  auc <- plotROCAndExtractAUC(pr$predictions, pr$responses, figpath, dataset_name = key)

  res <- map_dfr(names(auc), function(thr) map_dfr(names(auc[[thr]]), function(m) {
    v <- auc[[thr]][[m]]
    if (is.null(v$auc)) return(NULL)
    data.frame(dataset = key, method = m, threshold = as.numeric(thr),
               auroc = v$auc, n_true = v$n_true, stringsAsFactors = FALSE)
  }))
  write.table(res, out_csv, sep = ",", row.names = FALSE,
              col.names = !file.exists(out_csv), append = file.exists(out_csv))
  message(sprintf("[%s] %s done — threshold-2 AUROCs: %s", format(Sys.time(), "%H:%M:%S"), key,
                  paste(sprintf("%s=%.3f", res$method[res$threshold == 2],
                                res$auroc[res$threshold == 2]), collapse = " ")))
  rm(seurat_object_oi, potentials, regulons, dec_by_reg); gc(verbose = FALSE)
}
message("ABLATION AUROC DONE")

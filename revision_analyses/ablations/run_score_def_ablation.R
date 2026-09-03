# R2.2 extension — alternative Decipher score DEFINITIONS on the CytoSig AUROC
# benchmark. All variants recomputed from the archived per-regulon tables
# (no refitting): published score = mean over TFs of imp.perm x sign(rho) x deltaPagoda.
#   Decipher     — published (reference, from decipher_scores_by_cluster)
#   SumOverTFs   — sum instead of mean (the R2.8 Methods-vs-code discrepancy)
#   NoSpearmanSign — mean(imp.perm x deltaPagoda), sign(rho) dropped
#   UnsignedMagnitude — mean(|score|), pure magnitude ranking
#   ImpOnly      — mean(imp.perm), deltaPagoda dropped too
# Usage: Rscript run_score_def_ablation.R  -> results/ablation_score_defs.csv

repo_path  <- "/Users/edgarbasto/Documents/Decipher-manuscript"
drive_root <- "/Volumes/MegaEdgar/aws_pull_20260813/results"
local_root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/aws_data"
out_dir    <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/ablations"
out_csv <- file.path(out_dir, "results", "ablation_score_defs.csv")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(magrittr); library(data.table)
  library(pROC); library(ggplot2); library(Seurat)
})
source(file.path(repo_path, "R/vis_data_wrangling.R"))
source(file.path(repo_path, "R/intercellular_tailoring.R"))

DATASETS <- c("5yr_pic"="5yr_pic","bcg"="BCG","cord_pic"="cord_pic","covid"="covid",
  "erp"="ERP","lupus"="lupus","sepsis"="sepsis","tnbc"="TNBC","cz_influenza"="cz_influenza",
  "cz_hpap_t1d_islets"="cz_hpap_t1d_islets","cz_hnscc_hpv"="cz_hnscc_hpv",
  "cz_human_kidney_v1.5"="cz_human_kidney_v1.5","cz_cf_bronchial_biopsy"="cz_cf_bronchial_biopsy",
  "SevCOVID_Azimuthl2"="SevCOVID_Azimuthl2","MilCOVID_Azimuthl2"="MilCOVID_Azimuthl2")

data_dir_for <- function(folder) {
  local <- file.path(local_root, paste0(folder, "_data"))
  if (dir.exists(local)) local else file.path(drive_root, folder, "data")
}

mapping_table <- read.csv(file.path(repo_path, "reference_data", "cytosig_mapping_table_ligands_genes.csv"))
L.set <- getForrestLRDatabase(file.path(repo_path, "reference_data", "connectomedb_forrest_lrc2p.csv")) %>%
  mutate(interaction = paste(ligand, receptor, sep = "-"), lr = interaction) %>% unique()

as_score_df <- function(df, cl) df %>% mutate(sender_cluster="mixed", receiver_cluster=cl) %>%
  select(interaction, score, ligand, receptor, sender_cluster, receiver_cluster)

done <- if (file.exists(out_csv)) unique(read.csv(out_csv)$dataset) else character(0)

for (key in names(DATASETS)) {
  if (key %in% done) { message(key, ": done"); next }
  folder <- DATASETS[[key]]; ddir <- data_dir_for(folder)
  zdir <- file.path(drive_root, folder, "cytosig", "0_outputs", "z_score")
  if (!dir.exists(zdir)) { message("SKIP ", key); next }
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), key))

  dec_by_cluster <- readRDS(file.path(ddir, "decipher_scores_by_cluster.rds"))
  dec_by_reg <- readRDS(file.path(ddir, "decipher_scores_by_regulon_and_cluster.rds"))
  seurat_object_oi <- readRDS(file.path(ddir, "pseudobulk_seurat.rds"))
  clusters <- names(dec_by_cluster)
  lig_map <- bind_rows(lapply(dec_by_cluster, function(d) d %>% select(interaction, ligand, receptor))) %>% distinct()

  variant <- function(fn) setNames(lapply(clusters, function(cl) {
    dec_by_reg[[cl]] %>% group_by(interaction) %>% summarize(score = fn(cur_data()), .groups="drop") %>%
      left_join(lig_map, by="interaction") %>% as_score_df(cl)
  }), clusters)

  results_to_compare <- list(
    Decipher = setNames(lapply(clusters, function(cl)
      dec_by_cluster[[cl]] %>% rename(score = decipher_score) %>% as_score_df(cl)), clusters),
    SumOverTFs        = variant(function(d) sum(d$decipher_score)),
    NoSpearmanSign    = variant(function(d) mean(d$imp.perm * d$regulon.val)),
    UnsignedMagnitude = variant(function(d) mean(abs(d$decipher_score))),
    ImpOnly           = variant(function(d) mean(d$imp.perm))
  )

  zf <- paste0(zdir, "/")
  cytosig_sig <- summarizeZScores(list.files(zf), zf, mapping_table)
  figpath <- file.path(out_dir, "figures_score_defs", key)
  dir.create(figpath, recursive = TRUE, showWarnings = FALSE)
  pr <- getPredictionsResponsesForMethods(results_to_compare, cytosig_sig, L.set = L.set,
                                          seurat_object_oi, figpath)
  auc <- plotROCAndExtractAUC(pr$predictions, pr$responses, figpath, dataset_name = key)
  res <- map_dfr(names(auc), function(thr) map_dfr(names(auc[[thr]]), function(m) {
    v <- auc[[thr]][[m]]
    if (is.null(v$auc)) return(NULL)
    data.frame(dataset=key, method=m, threshold=as.numeric(thr), auroc=v$auc, n_true=v$n_true)
  }))
  write.table(res, out_csv, sep=",", row.names=FALSE,
              col.names=!file.exists(out_csv), append=file.exists(out_csv))
  message(sprintf("  thr2: %s", paste(sprintf("%s=%.3f", res$method[res$threshold==2], res$auroc[res$threshold==2]), collapse=" ")))
  rm(seurat_object_oi, dec_by_reg); gc(verbose=FALSE)
}
message("SCORE DEF ABLATION DONE")

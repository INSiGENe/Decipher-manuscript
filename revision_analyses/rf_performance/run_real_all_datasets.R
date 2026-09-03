# R2.4 extension (22 Aug 2026, overnight) — model-performance sweep over the
# remaining benchmark datasets: OOB / 5-fold CV / 50-permutation null /
# 10-refit stability for RF, plus lm / ridge / univariate-Spearman CV.
# NO donor hold-out (no donor lookups assembled for these datasets).
# 23 Aug: permutation null computed for RF only (null_methods = "rf") — the
# comparator nulls cost ~750 CV fits per model and are never reported.
# cord_pic / erp / tnbc were completed before this change (comparator nulls present, unused).
# Usage: Rscript run_real_all_datasets.R [n_cores]
# One CSV per dataset in results/sweep/<key>_per_tf_metrics.csv, resumable.
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source("R/regressors.R"); source("R/evaluate.R")
suppressMessages({library(randomForest); library(glmnet); library(parallel)})
args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, detectCores() - 2L)
dir.create("results/sweep", recursive = TRUE, showWarnings = FALSE)

drive_root <- "/Volumes/MegaEdgar/aws_pull_20260813/results"
local_root <- "../aws_data"
# priority: Fig 2e datasets not yet done, then the four without a threshold-2 slice
KEYS <- c("cord_pic" = "cord_pic", "erp" = "ERP", "tnbc" = "TNBC", "5yr_pic" = "5yr_pic",
          "bcg" = "BCG", "cz_influenza" = "cz_influenza", "MilCOVID_Azimuthl2" = "MilCOVID_Azimuthl2",
          "cz_cf_bronchial_biopsy" = "cz_cf_bronchial_biopsy", "cz_hnscc_hpv" = "cz_hnscc_hpv",
          "lupus" = "lupus", "sepsis" = "sepsis", "cz_hpap_t1d_islets" = "cz_hpap_t1d_islets",
          "cz_human_kidney_v1.5" = "cz_human_kidney_v1.5")
data_dir_for <- function(folder) {
  # 23 Aug: prefer the local archive subset (aws_data/results/<ds>/data), then the
  # older aws_data/<ds>_data stores, and only then the drive.
  cands <- c(file.path(local_root, "results", folder, "data"),
             file.path(local_root, paste0(folder, "_data")),
             file.path(drive_root, folder, "data"))
  cands[dir.exists(cands)][1]
}

# loader without donor structure (mirrors load_data.R minus membership/donor)
load_basic <- function(data_dir) {
  potentials_all <- readRDS(file.path(data_dir, "interaction_potentials_matrix_clusters_all_clusters.rds"))
  regulons_all   <- readRDS(file.path(data_dir, "regulon_scores_by_cluster.rds"))
  sig_all        <- readRDS(file.path(data_dir, "significant_regulons_by_cluster.rds"))
  pb             <- readRDS(file.path(data_dir, "pseudobulk_seurat.rds"))
  pb_meta        <- attr(pb, "meta.data"); rm(pb)
  out <- lapply(names(potentials_all), function(cl) {
    pot <- potentials_all[[cl]]; ids <- colnames(pot)
    if (is.null(sig_all[[cl]]) || is.null(regulons_all[[cl]])) return(NULL)
    sig_tfs <- sig_all[[cl]]$name[sig_all[[cl]]$class == "real"]
    sig_tfs <- sig_tfs[sig_tfs %in% rownames(regulons_all[[cl]])]
    if (!length(sig_tfs) || !all(ids %in% colnames(regulons_all[[cl]])) || !all(ids %in% rownames(pb_meta))) return(NULL)
    tf_activity <- regulons_all[[cl]][sig_tfs, ids, drop = FALSE]
    meta <- data.frame(metacell_id = ids, cluster = pb_meta[ids, "cluster"],
                       condition = pb_meta[ids, "condition"], stringsAsFactors = FALSE)
    list(cluster = cl, potentials = pot, tf_activity = tf_activity, meta = meta, membership = NULL)
  })
  names(out) <- names(potentials_all); Filter(Negate(is.null), out)
}

for (key in names(KEYS)) {
  out_csv <- sprintf("results/sweep/%s_per_tf_metrics.csv", key)
  dd <- data_dir_for(KEYS[[key]])
  if (!dir.exists(dd)) { message(key, ": data dir missing, skipping (", dd, ")"); next }
  dataset <- tryCatch(load_basic(dd), error = function(e) { message(key, ": LOAD FAIL ", conditionMessage(e)); NULL })
  if (is.null(dataset)) next
  done <- if (file.exists(out_csv)) unique(paste(read.csv(out_csv)$cluster, read.csv(out_csv)$tf)) else character(0)
  ord <- names(dataset)[order(sapply(dataset, function(d) nrow(d$tf_activity) * ncol(d$potentials)))]
  message(sprintf("[%s] === %s: %d clusters, %d models (data: %s)", format(Sys.time(), "%H:%M:%S"),
                  key, length(dataset), sum(sapply(dataset, function(d) nrow(d$tf_activity))), dd))
  for (cl in ord) {
    inp <- dataset[[cl]]
    tfs <- rownames(inp$tf_activity); tfs <- tfs[!paste(cl, tfs) %in% done]
    if (!length(tfs)) next
    message(sprintf("[%s] %s / %s: %d TFs (%d obs x %d preds)", format(Sys.time(), "%H:%M:%S"), key, cl,
                    length(tfs), ncol(inp$potentials), nrow(inp$potentials)))
    res <- mclapply(tfs, function(tf) tryCatch(evaluate_tf(tf, inp, null_methods = "rf"),
              error = function(e) { message("FAIL ", key, "/", cl, "/", tf, ": ", conditionMessage(e)); NULL }),
              mc.cores = n_cores, mc.preschedule = FALSE)
    res <- do.call(rbind, Filter(Negate(is.null), res))
    if (!is.null(res)) { res$dataset <- key
      # 23 Aug: retry the append — a late SIGCHLD from a worker once hit file() with EINTR and killed the run
      for (attempt in 1:10) {
        ok <- tryCatch({ write.table(res, out_csv, sep = ",", row.names = FALSE,
                                     col.names = !file.exists(out_csv), append = file.exists(out_csv)); TRUE },
                       error = function(e) { message("write attempt ", attempt, " failed: ", conditionMessage(e)); FALSE },
                       warning = function(w) { message("write attempt ", attempt, " warning: ", conditionMessage(w)); FALSE })
        if (ok) break; Sys.sleep(2)
      }
      if (!ok) stop("could not write ", out_csv, " after 10 attempts") }
  }
  message(sprintf("[%s] === %s DONE", format(Sys.time(), "%H:%M:%S"), key))
  rm(dataset); gc(verbose = FALSE)
}
message("SWEEP ALL DONE")

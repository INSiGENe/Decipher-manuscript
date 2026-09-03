# R2.4 ntree sensitivity sweep — importance stability + OOB R^2 at
# ntree = 100 (manuscript) / 500 / 1000, for every TF whose RF model beat the
# permutation null in the main run (those are the models whose edges get
# mechanistically interpreted; pilot flagged stability ~0.56 at ntree=100).
# Usage: Rscript run_ntree_sweep.R [n_cores]   (run AFTER run_real_covid.R)
# Writes results/real_covid_ntree_sweep.csv incrementally, resuming on restart.

setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source("R/regressors.R"); source("R/evaluate.R"); source("R/load_data.R")
suppressMessages({library(randomForest); library(parallel)})

args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, detectCores() - 2L)
main_csv <- "results/real_covid_per_tf_metrics.csv"
out_csv  <- "results/real_covid_ntree_sweep.csv"
NTREES <- c(100, 500, 1000)

stopifnot(file.exists(main_csv))
main <- read.csv(main_csv, stringsAsFactors = FALSE)
targets <- main[main$method == "rf" & main$beats_null, c("cluster", "tf")]
message(nrow(targets), " beat-null RF models to sweep")

done <- if (file.exists(out_csv)) {
  prev <- read.csv(out_csv, stringsAsFactors = FALSE)
  unique(paste(prev$cluster, prev$tf))
} else character(0)

dataset <- load_decipher_dataset(
  data_dir         = "../aws_data/covid_data",
  membership_csv   = "../donor_purity/metacell_membership.csv",
  donor_lookup_csv = "../donor_purity/cell_donor_lookup.csv")

sweep_tf <- function(inp, tf) {
  x <- t(inp$potentials)
  y <- as.numeric(inp$tf_activity[tf, ])
  do.call(rbind, lapply(NTREES, function(nt) {
    fit  <- fit_model(x, y, "rf", seed = 123, ntree = nt)
    stab <- importance_stability(x, y, n_rep = 10, seed = 123, ntree = nt)
    data.frame(cluster = inp$cluster, tf = tf, ntree = nt,
               oob_r2 = pooled_r2(y, fit$oob_pred),
               imp_stability_spearman = stab$imp_spearman,
               imp_top5_jaccard = stab$top5_jaccard,
               stringsAsFactors = FALSE)
  }))
}

for (cl in unique(targets$cluster)) {
  inp <- dataset[[cl]]
  tfs <- setdiff(targets$tf[targets$cluster == cl],
                 sub("^\\S+ ", "", grep(paste0("^", cl, " "), done, value = TRUE)))
  if (!length(tfs)) { message(cl, ": already complete"); next }
  message(sprintf("[%s] %s: %d TFs on %d cores", format(Sys.time(), "%H:%M:%S"), cl, length(tfs), n_cores))
  res_list <- mclapply(tfs, function(tf) {
    tryCatch(sweep_tf(inp, tf),
             error = function(e) { message("FAIL ", cl, "/", tf, ": ", conditionMessage(e)); NULL })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  res <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (!is.null(res))
    write.table(res, out_csv, sep = ",", row.names = FALSE,
                col.names = !file.exists(out_csv), append = file.exists(out_csv))
  message(sprintf("[%s] %s done", format(Sys.time(), "%H:%M:%S"), cl))
}
message("SWEEP DONE")

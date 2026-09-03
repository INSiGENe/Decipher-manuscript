# R2.4 real-data run — vaccination dataset (covid), archive RDS inputs.
# Usage: Rscript run_real_covid.R [n_cores]
# Writes results/real_covid_per_tf_metrics.csv incrementally (one row per
# TF x method), resuming past completed (cluster, tf) pairs on restart.

setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source("R/regressors.R"); source("R/evaluate.R"); source("R/load_data.R")
suppressMessages({library(randomForest); library(glmnet); library(parallel)})

args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, detectCores() - 2L)
out_csv <- "results/real_covid_per_tf_metrics.csv"
dir.create("results", showWarnings = FALSE)

dataset <- load_decipher_dataset(
  data_dir         = "../aws_data/covid_data",
  membership_csv   = "../donor_purity/metacell_membership.csv",
  donor_lookup_csv = "../donor_purity/cell_donor_lookup.csv")

done <- if (file.exists(out_csv)) {
  prev <- read.csv(out_csv, stringsAsFactors = FALSE)
  unique(paste(prev$cluster, prev$tf))
} else character(0)

# smallest clusters first so results land early
ord <- names(dataset)[order(sapply(dataset, function(d) nrow(d$tf_activity) * ncol(d$potentials)))]

for (cl in ord) {
  inp <- dataset[[cl]]
  tfs <- setdiff(rownames(inp$tf_activity), sub("^\\S+ ", "", grep(paste0("^", cl, " "), done, value = TRUE)))
  tfs <- rownames(inp$tf_activity)[rownames(inp$tf_activity) %in% tfs]
  if (!length(tfs)) { message(cl, ": already complete"); next }
  message(sprintf("[%s] %s: %d TFs on %d cores", format(Sys.time(), "%H:%M:%S"), cl, length(tfs), n_cores))

  res_list <- mclapply(tfs, function(tf) {
    tryCatch(evaluate_tf(tf, inp),
             error = function(e) { message("FAIL ", cl, "/", tf, ": ", conditionMessage(e)); NULL })
  }, mc.cores = n_cores, mc.preschedule = FALSE)

  res <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (!is.null(res))
    write.table(res, out_csv, sep = ",", row.names = FALSE,
                col.names = !file.exists(out_csv), append = file.exists(out_csv))
  message(sprintf("[%s] %s done (%d/%d TFs ok)", format(Sys.time(), "%H:%M:%S"), cl,
                  length(unique(res$tf)), length(tfs)))
}
message("ALL DONE")

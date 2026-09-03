# R2.4 follow-up (21 Aug 2026): does the importance RANKING at ntree=100 (the
# manuscript setting) match the ranking at ntree=1000, per model? The sweep
# measured within-ntree stability only; this is the direct cross-ntree test
# the "no conclusion depends on ntree" sentence needs. Both COVID datasets.
# Usage: Rscript run_ntree_ranking_concordance.R [n_cores]
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source("R/regressors.R"); source("R/evaluate.R"); source("R/load_data.R")
suppressMessages({library(randomForest); library(parallel)})
args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else max(1L, detectCores() - 2L)
out_csv <- "results/ntree_ranking_concordance.csv"

SETS <- list(
  vaccination = list(main = "results/real_covid_per_tf_metrics.csv",
                     data_dir = "../aws_data/covid_data",
                     membership_csv = "../donor_purity/metacell_membership.csv",
                     donor_lookup_csv = "../donor_purity/cell_donor_lookup.csv"),
  severity    = list(main = "results/real_severity_per_tf_metrics.csv",
                     data_dir = "../aws_data/SevCOVID_Azimuthl2_data",
                     membership_csv = "../donor_purity/severity/metacell_membership.csv",
                     donor_lookup_csv = "../donor_purity/severity/cell_donor_lookup.csv"))

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
one_tf <- function(inp, tf, ds) {
  x <- t(inp$potentials); y <- as.numeric(inp$tf_activity[tf, ])
  # two independent seeds at 100 (noise floor) and one fit at 1000
  i100a <- fit_model(x, y, "rf", seed = 123,  ntree = 100)$importance
  i100b <- fit_model(x, y, "rf", seed = 4567, ntree = 100)$importance
  i1000 <- fit_model(x, y, "rf", seed = 123,  ntree = 1000)$importance
  top <- function(v, n) order(v, decreasing = TRUE)[seq_len(n)]
  data.frame(dataset = ds, cluster = inp$cluster, tf = tf, n_pred = ncol(x),
    rho_100_vs_1000   = suppressWarnings(cor(i100a, i1000, method = "spearman")),
    rho_100_vs_100    = suppressWarnings(cor(i100a, i100b, method = "spearman")),
    top5_jac_100_1000 = jac(top(i100a, 5),  top(i1000, 5)),
    top5_jac_100_100  = jac(top(i100a, 5),  top(i100b, 5)),
    top10_jac_100_1000 = jac(top(i100a, 10), top(i1000, 10)),
    top1_same_100_1000 = top(i100a, 1) == top(i1000, 1),
    stringsAsFactors = FALSE)
}

done <- if (file.exists(out_csv)) {
  p <- read.csv(out_csv, stringsAsFactors = FALSE); unique(paste(p$dataset, p$cluster, p$tf))
} else character(0)

for (ds in names(SETS)) {
  S <- SETS[[ds]]
  main <- read.csv(S$main, stringsAsFactors = FALSE)
  targets <- main[main$method == "rf", c("cluster", "tf")]
  dataset <- load_decipher_dataset(S$data_dir, S$membership_csv, S$donor_lookup_csv)
  for (cl in unique(targets$cluster)) {
    inp <- dataset[[cl]]
    tfs <- targets$tf[targets$cluster == cl]
    tfs <- tfs[!paste(ds, cl, tfs) %in% done]
    if (!length(tfs)) next
    message(sprintf("[%s] %s / %s: %d TFs", format(Sys.time(), "%H:%M:%S"), ds, cl, length(tfs)))
    res <- mclapply(tfs, function(tf) tryCatch(one_tf(inp, tf, ds),
              error = function(e) { message("FAIL ", cl, "/", tf, ": ", conditionMessage(e)); NULL }),
              mc.cores = n_cores, mc.preschedule = FALSE)
    res <- do.call(rbind, Filter(Negate(is.null), res))
    if (!is.null(res)) write.table(res, out_csv, sep = ",", row.names = FALSE,
      col.names = !file.exists(out_csv), append = file.exists(out_csv))
  }
}
message("DONE")

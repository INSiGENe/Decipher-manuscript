# Regressor panel (24 Aug 2026): apples-to-apples comparison of importance-capable
# regressors on the LR->TF task. Methods: ranger (extratrees), xgboost, MARS (earth),
# RBF-SVM (e1071). ONE importance definition for all: out-of-fold permutation
# importance (deltaMSE on the held-out fold, averaged over 5 folds), i.e. the RF's own
# principle applied model-agnostically. Outputs per dataset:
#   results/regressor_panel/<key>_cv.csv     (per TF x method: cv_r2)
#   results/regressor_panel/<key>_edges.rds  (per method: interaction-level scores,
#     imp x sign(rho) x regulon_val, mean over TFs -- ready for AUROC injection,
#     which additionally needs the drive's cytosig z_scores; run separately)
# Usage: Rscript run_regressor_panel.R [n_cores] [keys...]
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
suppressMessages({library(dplyr); library(ranger); library(xgboost); library(earth); library(e1071); library(parallel)})
args <- commandArgs(TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 6L
local_root <- "../aws_data/results"
KEYS <- c("covid","SevCOVID_Azimuthl2","MilCOVID_Azimuthl2","cord_pic"="cord_pic","erp"="ERP","tnbc"="TNBC",
          "5yr_pic"="5yr_pic","bcg"="BCG","cz_influenza"="cz_influenza",
          "cz_cf_bronchial_biopsy"="cz_cf_bronchial_biopsy","cz_hnscc_hpv"="cz_hnscc_hpv")
KEYS <- setNames(ifelse(KEYS == "", names(KEYS), KEYS), ifelse(names(KEYS) == "", KEYS, names(KEYS)))
if (length(args) >= 2) KEYS <- KEYS[args[-1]]
dir.create("results/regressor_panel", recursive = TRUE, showWarnings = FALSE)

fitfun <- list(
  ranger_et = function(xtr, ytr) { m <- ranger(x = as.data.frame(xtr), y = ytr, num.trees = 300, splitrule = "extratrees", num.threads = 1, seed = 123); function(xx) predict(m, as.data.frame(xx), num.threads = 1)$predictions },
  xgboost   = function(xtr, ytr) { m <- xgboost(data = xtr, label = ytr, nrounds = 200, max_depth = 4, eta = 0.1, subsample = 0.8, nthread = 1, verbose = 0); function(xx) predict(m, xx) },
  mars      = function(xtr, ytr) { m <- earth(x = xtr, y = ytr, degree = 2); function(xx) as.numeric(predict(m, xx)) },
  svm_rbf   = function(xtr, ytr) { m <- svm(x = xtr, y = ytr, kernel = "radial"); function(xx) as.numeric(predict(m, xx)) })

eval_tf <- function(x, y, cond) {
  set.seed(123)
  fold <- integer(length(y)); for (cc in unique(cond)) { idx <- sample(which(cond == cc)); fold[idx] <- rep_len(1:5, length(idx)) }
  out <- list()
  for (mname in names(fitfun)) {
    pred <- numeric(length(y)); imp <- matrix(0, 5, ncol(x))
    okk <- TRUE
    for (k in 1:5) {
      tr <- fold != k; te <- fold == k
      f <- tryCatch(fitfun[[mname]](x[tr, , drop = FALSE], y[tr]), error = function(e) NULL)
      if (is.null(f)) { okk <- FALSE; break }
      pred[te] <- f(x[te, , drop = FALSE])
      base_mse <- mean((y[te] - pred[te])^2)
      set.seed(1000 + k)
      for (j in seq_len(ncol(x))) { xp <- x[te, , drop = FALSE]; xp[, j] <- sample(xp[, j]); imp[k, j] <- mean((y[te] - f(xp))^2) - base_mse }
    }
    if (!okk) { out[[mname]] <- list(cv = NA_real_, imp = rep(NA_real_, ncol(x))); next }
    out[[mname]] <- list(cv = 1 - mean((y - pred)^2) / mean((y - mean(y))^2), imp = colMeans(imp))
  }
  out
}

for (key in names(KEYS)) {
  cvf <- sprintf("results/regressor_panel/%s_cv.csv", key)
  edf <- sprintf("results/regressor_panel/%s_edges.rds", key)
  if (file.exists(edf)) { message(key, ": done, skipping"); next }
  dd <- file.path(local_root, KEYS[[key]], "data")
  if (!dir.exists(dd)) { message("SKIP ", key, ": no data dir"); next }
  message(sprintf("[%s] === %s", format(Sys.time(), "%H:%M:%S"), key))
  pot <- readRDS(file.path(dd, "interaction_potentials_matrix_clusters_all_clusters.rds"))
  reg <- readRDS(file.path(dd, "regulon_scores_by_cluster.rds"))
  sig <- readRDS(file.path(dd, "significant_regulons_by_cluster.rds"))
  der <- readRDS(file.path(dd, "decipher_scores_by_regulon_and_cluster.rds"))
  pb  <- readRDS(file.path(dd, "pseudobulk_seurat.rds")); pmeta <- attr(pb, "meta.data"); rm(pb)
  cv_rows <- list(); edge_scores <- list()
  for (cl in names(pot)) {
    if (is.null(sig[[cl]]) || is.null(reg[[cl]])) next
    tfs <- sig[[cl]]$name[sig[[cl]]$class == "real"]; tfs <- tfs[tfs %in% rownames(reg[[cl]])]
    if (!length(tfs)) next
    x <- t(pot[[cl]]); y_all <- reg[[cl]][tfs, colnames(pot[[cl]]), drop = FALSE]
    cond <- pmeta[colnames(pot[[cl]]), "condition"]
    message(sprintf("[%s] %s / %s: %d TFs (%d x %d)", format(Sys.time(), "%H:%M:%S"), key, cl, length(tfs), nrow(x), ncol(x)))
    res <- mclapply(tfs, function(tf) tryCatch(eval_tf(x, as.numeric(y_all[tf, ]), cond), error = function(e) NULL),
                    mc.cores = n_cores, mc.preschedule = FALSE)
    names(res) <- tfs; res <- Filter(Negate(is.null), res)
    d <- der[[cl]]; if (!is.data.frame(d)) d <- do.call(rbind, d)
    einfo <- d %>% select(interaction, regulon, spearman.cor, regulon.val) %>% distinct()
    for (mname in names(fitfun)) {
      cv_rows[[paste(cl, mname)]] <- data.frame(dataset = key, cluster = cl, tf = names(res), method = mname,
        cv_r2 = sapply(res, function(r) r[[mname]]$cv), stringsAsFactors = FALSE)
      imps <- bind_rows(lapply(names(res), function(tf) data.frame(interaction = colnames(x), regulon = tf, imp = res[[tf]][[mname]]$imp)))
      edge_scores[[paste(cl, mname)]] <- imps %>% inner_join(einfo, by = c("interaction", "regulon")) %>%
        mutate(contrib = imp * sign(spearman.cor) * regulon.val) %>%
        group_by(interaction) %>% summarize(score = mean(contrib), .groups = "drop") %>%
        mutate(cluster = cl, method = mname)
    }
  }
  write.csv(bind_rows(cv_rows), cvf, row.names = FALSE)
  saveRDS(bind_rows(edge_scores), edf)
  message(sprintf("[%s] === %s DONE", format(Sys.time(), "%H:%M:%S"), key))
}
message("PANEL ALL DONE")

# Core evaluation machinery for R2.4 (out-of-sample performance of LR->TF
# models, null-model comparison, donor-preserving CV, importance stability).
#
# Canonical input (one cluster) — see R/load_data.R for how to build it:
#   potentials : matrix, LR pairs (rows) x meta-cells (cols)
#   tf_activity: matrix, TFs (rows) x meta-cells (cols), same columns
#   meta       : data.frame(metacell_id, condition, seed_donor [optional])
#   membership : optional data.frame(metacell_id, cell_id, donor) for
#                donor-aware fold construction (drop contaminated meta-cells)

pooled_r2 <- function(y, pred) {
  ok <- is.finite(pred)
  if (sum(ok) < 3 || var(y[ok]) == 0) return(NA_real_)
  1 - mean((y[ok] - pred[ok])^2) / mean((y[ok] - mean(y[ok]))^2)
}

# k-fold assignment stratified by condition (so no fold is single-condition)
make_cv_folds <- function(meta, k = 5, seed = 123) {
  set.seed(seed)
  fold <- integer(nrow(meta))
  for (cond in unique(meta$condition)) {
    idx <- sample(which(meta$condition == cond))
    fold[idx] <- rep_len(seq_len(k), length(idx))
  }
  fold
}

# Leave-one-donor-out folds. Test = meta-cells seeded by the held-out donor.
# With a membership table, training additionally DROPS every meta-cell that
# contains any cell from the held-out donor (the conservative rule from the
# R2.3 feasibility notes) — report how many were dropped.
make_lodo_folds <- function(meta, membership = NULL) {
  stopifnot(!is.null(meta$seed_donor))
  lapply(unique(meta$seed_donor), function(d) {
    test <- which(meta$seed_donor == d)
    if (!is.null(membership)) {
      contaminated <- unique(membership$metacell_id[membership$donor == d])
      train <- which(!(meta$metacell_id %in% contaminated))
      n_dropped <- nrow(meta) - length(train) - length(test)
      # meta-cells seeded by another donor but containing d's cells are
      # excluded from BOTH sets; count them as dropped
      n_dropped <- sum(!(seq_len(nrow(meta)) %in% c(train, test)))
    } else {
      train <- which(meta$seed_donor != d)
      n_dropped <- 0L
    }
    list(donor = d, train = train, test = test, n_dropped = n_dropped)
  })
}

# Generic out-of-fold evaluation. folds: either an integer fold vector
# (k-fold) or a list of train/test index lists (LODO).
cv_r2 <- function(x, y, folds, method, seed = 123, ntree = 100) {
  pred <- rep(NA_real_, length(y))
  if (is.list(folds)) {
    for (f in folds) {
      if (length(f$train) < 10 || length(f$test) < 2) next
      fit <- fit_model(x[f$train, , drop = FALSE], y[f$train], method, seed, ntree)
      pred[f$test] <- fit$predict(x[f$test, , drop = FALSE])
    }
  } else {
    for (k in unique(folds)) {
      tr <- which(folds != k); te <- which(folds == k)
      fit <- fit_model(x[tr, , drop = FALSE], y[tr], method, seed, ntree)
      pred[te] <- fit$predict(x[te, , drop = FALSE])
    }
  }
  pooled_r2(y, pred)
}

# Permutation null: shuffle y, recompute the method's headline metric.
# rf uses OOB R^2 (cheap, and R2.4 explicitly asks for OOB); others use k-fold CV.
null_distribution <- function(x, y, meta, method, n_perm = 50, k = 5,
                              seed = 123, ntree = 100) {
  vapply(seq_len(n_perm), function(i) {
    set.seed(seed + i)
    yp <- sample(y)
    if (method == "rf") {
      fit <- fit_model(x, yp, "rf", seed = seed + i, ntree = ntree)
      pooled_r2(yp, fit$oob_pred)
    } else {
      cv_r2(x, yp, make_cv_folds(meta, k, seed + i), method, seed + i, ntree)
    }
  }, numeric(1))
}

# Importance stability (rf only): refit with different seeds, report mean
# pairwise Spearman correlation of importance vectors + mean pairwise
# Jaccard overlap of the top-5 predictors.
importance_stability <- function(x, y, n_rep = 10, seed = 123, ntree = 100, top_n = 5) {
  imps <- sapply(seq_len(n_rep), function(i)
    fit_model(x, y, "rf", seed = seed + 1000 * i, ntree = ntree)$importance)
  pairs <- combn(n_rep, 2)
  rho <- mean(apply(pairs, 2, function(p)
    suppressWarnings(cor(imps[, p[1]], imps[, p[2]], method = "spearman"))), na.rm = TRUE)
  tops <- apply(imps, 2, function(v) order(v, decreasing = TRUE)[seq_len(top_n)])
  jac <- mean(apply(pairs, 2, function(p) {
    a <- tops[, p[1]]; b <- tops[, p[2]]
    length(intersect(a, b)) / length(union(a, b))
  }))
  list(imp_spearman = rho, top5_jaccard = jac, mean_importance = rowMeans(imps))
}

# Evaluate one TF across methods. Returns a data.frame, one row per method.
evaluate_tf <- function(tf_name, input, methods = c("rf", "lm", "ridge", "spearman"),
                        n_perm = 50, n_stab = 10, k = 5, seed = 123, ntree = 100,
                        null_methods = methods) {   # 23 Aug: restrict nulls (e.g. "rf") — comparator nulls are never reported
  x <- t(input$potentials)                      # n_obs x n_pred
  y <- as.numeric(input$tf_activity[tf_name, ])
  meta <- input$meta
  has_donor <- !is.null(meta$seed_donor)
  kfolds <- make_cv_folds(meta, k, seed)
  lfolds <- if (has_donor) make_lodo_folds(meta, input$membership) else NULL

  do.call(rbind, lapply(methods, function(m) {
    full <- fit_model(x, y, m, seed, ntree)
    oob  <- if (m == "rf") pooled_r2(y, full$oob_pred) else NA_real_
    cv   <- cv_r2(x, y, kfolds, m, seed, ntree)
    lodo <- if (has_donor) cv_r2(x, y, lfolds, m, seed, ntree) else NA_real_
    nulls <- if (m %in% null_methods) null_distribution(x, y, meta, m, n_perm, k, seed, ntree) else NA_real_
    obs  <- if (m == "rf") oob else cv
    p_perm <- if (m %in% null_methods) (1 + sum(nulls >= obs, na.rm = TRUE)) / (1 + sum(is.finite(nulls))) else NA_real_
    stab <- if (m == "rf" && n_stab > 1) importance_stability(x, y, n_stab, seed, ntree)
            else list(imp_spearman = NA_real_, top5_jaccard = NA_real_)
    data.frame(
      cluster = input$cluster, tf = tf_name, method = m,
      n_obs = length(y), n_pred = ncol(x),
      oob_r2 = oob, cv_r2 = cv, lodo_r2 = lodo,
      null_metric = if (m == "rf") "oob_r2" else "cv_r2",
      null_q95 = unname(quantile(nulls, 0.95, na.rm = TRUE)),
      p_perm = p_perm, beats_null = p_perm < 0.05,
      imp_stability_spearman = stab$imp_spearman,
      imp_top5_jaccard = stab$top5_jaccard,
      stringsAsFactors = FALSE
    )
  }))
}

# Convenience: evaluate every TF in the input
evaluate_cluster <- function(input, ...) {
  res <- do.call(rbind, lapply(rownames(input$tf_activity), evaluate_tf, input = input, ...))
  rownames(res) <- NULL
  res
}

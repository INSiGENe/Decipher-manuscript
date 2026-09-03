# Pluggable regressor interface for the R2.4 evaluation harness.
#
# Every method exposes the same contract:
#   fit_model(x, y, method, seed, ntree) -> list(
#     predict    : function(newx) -> numeric predictions
#     oob_pred   : OOB predictions (rf only, else NULL)
#     importance : permutation importance vector (rf only, else NULL)
#   )
# x: numeric matrix, n_obs (meta-cells) x n_pred (LR pairs). y: numeric vector.
#
# "rf" mirrors the manuscript pipeline exactly (randomForest, ntree = 100,
# importance = TRUE, type-1 unscaled extraction — Decipher-manuscript
# R/decipher_scoring.R:136-140, :40).
# "lm" / "ridge" are the R2.2 ablation candidates ("linear/regularized
# regression instead of RF").
# "spearman" is the simplest ablation baseline: predict from the single
# training-set LR pair with the largest |Spearman correlation| to the TF
# (univariate linear fit) — the predictive analogue of "simple LR-TF
# Spearman correlation instead of RF regression".

fit_model <- function(x, y, method = c("rf", "lm", "ridge", "spearman"),
                      seed = 123, ntree = 100) {
  method <- match.arg(method)
  stopifnot(is.matrix(x), nrow(x) == length(y))
  # internal stable column names so formula interfaces never choke
  colnames(x) <- paste0("V", seq_len(ncol(x)))
  set.seed(seed)

  if (method == "rf") {
    m <- randomForest::randomForest(x = x, y = y, ntree = ntree, importance = TRUE)
    imp <- as.numeric(randomForest::importance(m, type = 1, scale = FALSE))
    names(imp) <- colnames(x)
    return(list(
      predict    = function(newx) {
        colnames(newx) <- paste0("V", seq_len(ncol(newx)))
        as.numeric(predict(m, newx))
      },
      oob_pred   = as.numeric(m$predicted),
      importance = imp
    ))
  }

  if (method == "lm") {
    df <- data.frame(y = y, x, check.names = FALSE)
    m <- suppressWarnings(lm(y ~ ., data = df))
    return(list(
      predict = function(newx) {
        colnames(newx) <- paste0("V", seq_len(ncol(newx)))
        suppressWarnings(as.numeric(predict(m, newdata = as.data.frame(newx))))
      },
      oob_pred = NULL, importance = NULL
    ))
  }

  if (method == "ridge") {
    m <- glmnet::cv.glmnet(x, y, alpha = 0, nfolds = 5)
    return(list(
      predict = function(newx) {
        colnames(newx) <- paste0("V", seq_len(ncol(newx)))
        as.numeric(predict(m, newx = newx, s = "lambda.min"))
      },
      oob_pred = NULL, importance = NULL
    ))
  }

  # spearman: univariate best-correlate baseline
  sp <- suppressWarnings(apply(x, 2, cor, y = y, method = "spearman"))
  sp[is.na(sp)] <- 0
  best <- which.max(abs(sp))
  xb <- x[, best]
  cf <- coef(lm(y ~ xb))
  list(
    predict = function(newx) {
      colnames(newx) <- paste0("V", seq_len(ncol(newx)))
      as.numeric(cf[1] + cf[2] * newx[, best])
    },
    oob_pred = NULL, importance = NULL,
    best_predictor = colnames(x)[best]
  )
}

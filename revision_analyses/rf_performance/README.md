# R2.4 evaluation harness (with R2.2 bridge)

Runnable skeleton for reviewer comment **R2.4** (out-of-sample performance of the
LR→TF random-forest models, null-model comparison, donor-preserving CV,
importance stability), designed so its pluggable regressor interface doubles as
the engine for the **R2.2** ablations (linear/ridge/Spearman instead of RF).

## Run the pilot

```sh
cd analysis_r2/rf_performance
Rscript run_pilot.R          # ~3 min; needs randomForest + glmnet (installed 13 Aug 2026)
```

Synthetic data with known ground truth (15 TFs: 6 linear-signal, 3
nonlinear-signal, 6 null; 40 LR pairs; 216 meta-cells; 6 donors with random
effects; 3 cells per meta-cell with ~20% donor contamination). Output:
`results/pilot_per_tf_metrics.csv` — one row per TF × method with `oob_r2`,
`cv_r2` (5-fold, condition-stratified), `lodo_r2` (leave-one-donor-out with
contaminated-meta-cell dropping), `null_q95`, `p_perm`, `beats_null`,
`imp_stability_spearman`, `imp_top5_jaccard`.

## Layout

| File | What |
|---|---|
| `R/regressors.R` | `fit_model(x, y, method)` — `rf` mirrors the manuscript pipeline exactly (randomForest, ntree=100, type-1 unscaled permutation importance); `lm` / `ridge` / `spearman` are the ablation candidates. **This interface is the R2.2 bridge** — swap variants into the full Decipher scoring later. |
| `R/evaluate.R` | Fold builders (stratified k-fold; LODO with the R2.3 contamination-drop rule), pooled R², permutation nulls, importance stability, `evaluate_cluster()`. |
| `R/synthetic_data.R` | Ground-truth generator shaped exactly like the real input. |
| `R/load_data.R` | **The plug-in point for real data** — documents which Decipher-manuscript pipeline objects to export (per-cluster interaction-potential matrix, real-regulon PAGODA2 scores, condition labels, and `analysis_r2/donor_purity/metacell_membership.csv` for donor folds). Loader is a stub until the RDS exports exist. |
| `run_pilot.R` | End-to-end pilot + console summary. |

## Pilot results (13 Aug 2026) — harness validation

- **Null test discriminates perfectly:** 100% of signal TFs beat the permutation
  null (rf/lm/ridge), 0% of null TFs do. The `beats_null` flag is usable as the
  R2.4 "drop/flag edges from models that don't beat null" criterion.
- **Donor machinery works and matters:** naive 5-fold CV is optimistic vs
  leave-one-donor-out (mean gap +0.044 R² on signal TFs at donor_sd=0.6);
  null TFs sink further below zero under LODO, as they should.
- **Preview of the R2.2 answer:** ridge matched/beat RF even on the
  "nonlinear" TFs (a product of log-normal potentials is ~log-linear), while
  the univariate Spearman baseline collapsed on multi-driver TFs. On real data
  this comparison IS the ablation result — whichever way it falls.
- **⚠ Importance stability is mediocre at the manuscript's ntree=100:** mean
  importance Spearman across 10 refits was only ~0.56 for signal TFs (~0.29
  for null), and only ~54% of true drivers appeared in the top-k importances
  (correlated predictors). Expect the reviewer's stability question to have
  teeth on real data; an ntree sensitivity sweep (100 vs 500 vs 1000) should
  probably be added when real data is wired in.

## Wiring real data (next steps)

1. Export per-cluster RDS from a `set.seed(123)` run of
   `6_decipher_pipeline_v1_modularized.R` (vaccination dataset first):
   `filtered_interaction_potentials_matrix_all_clusters` and the real-regulon
   PAGODA2 score matrices. Implement `load_decipher_cluster()`.
2. Donor folds: vaccination via `analysis_r2/donor_purity/metacell_membership.csv`
   (reconcile against revision-1 AWS outputs first — R2.3 caveat); severity
   after the AWS download.
3. Headline numbers for the response: proportion of TF models with p_perm<0.05,
   OOB/CV/LODO R² distributions, importance stability, RF-vs-ridge/lm/Spearman
   — then feed the surviving-TF filter back into the Decipher scores.
4. R2.2: reuse `fit_model()` variants inside the pipeline's
   `getRandomForestWeightsAllClusters()` to produce ablated Decipher scores for
   the benchmark AUROC comparison; add the one non-regressor ablation
   (LR-potential-only, no TF layer) separately.

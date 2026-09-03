# Revision analyses

Scripts and summary outputs for the analyses added during the second round of peer
review (Cell Communication and Signaling, 2026). Every table and figure in the
revised supplementary material (Supplementary Note 3, Supplementary Tables ST.1–ST.11,
Supplementary Figures SF.3–SF.6) is produced by a script in this folder from the
archived outputs of the manuscript pipeline (the `results/<comparison>/data/` stores
written by the pipeline scripts in `../scripts/`), without re-running the pipeline
itself, except where stated.

Paths are defined at the top of each script and point to the local archive layout
`results/<comparison>/{data, pre_processing, cellOracle/data}`; set them to your copy of
the archived results before running. Package versions: R 4.4–4.6 with randomForest
4.7-1.1, glmnet, ranger, xgboost, earth, e1071, pROC, dplyr/tidyr, Seurat 4/5;
CellOracle 0.19.0 in the pinned Docker image (`../DockerFiles`). Small result tables
(CSV) are included so every reported number can be checked without recomputation.

| Folder | Script(s) | Produces | Cited as |
|---|---|---|---|
| `rf_performance/` | `R/evaluate.R`, `R/regressors.R`, `R/load_data.R`; `run_real_covid.R`, `run_real_severity.R` (both COVID-19 datasets, donor-held-out folds); `run_real_all_datasets.R` (13 remaining comparisons); `run_ntree_sweep.R`, `run_ntree_ranking_concordance.R` | out-of-bag / cross-validated / donor-held-out R², 50-permutation nulls, importance stability, tree-count sensitivity; `results/*.csv`, `results/sweep/*.csv` | Supplementary Note 3; ST.3, ST.4; SF.3; Supplementary Data 1–2 |
| `ablations/` | `run_ablation_auroc.R`, `run_ablation_auroc_ols.R`, `run_score_def_ablation.R` (variants injected into the manuscript's CytoSig AUROC benchmark, `../R/vis_data_wrangling.R`); `run_regressor_panel.R` (extra-trees, XGBoost, MARS, RBF-SVM under one out-of-fold permutation-importance definition) | `results/ablation_auroc*.csv`, `results/ablation_score_defs.csv`, `results/regressor_panel/*_cv.csv` | ST.6; SF.5; Methods *Ablation of framework components* |
| `donor_purity/` | `vaccination_donor_purity.R`, `severity_donor_purity.R` (instrumented copies of the pipeline's meta-cell functions; memberships verified identical to the archived pseudobulk) | donor composition, seed-donor enrichment vs a donor-blind null, effective donor numbers; also the membership tables used for donor-held-out folds (not included: large) | ST.11(a) |
| `c8_alt_reference/` | `run_c8_altref.R <control_cluster>` (full pipeline rerun with CD16⁺ monocytes or cDC2 as the C8 reference), `compare_references.R` | `outputs/reference_comparison_summary.csv` | ST.7; Methods *Adapted Decipher analysis* |
| `celloracle_perturbation/` | `make_tf_ranking.R` (TF panels from the archived per-regulon scores); `ko_cd16mono.py`, `ko_cd16_vax.py` (CellOracle GRN reproduction and TF knockouts, severity and vaccination CD16⁺ monocytes); `ko_extra_panels.py` (activity-change-only and importance-only panels); `anchor_regulon_orientation.R`, `verify_target_anchor.R`, `anchor_all_datasets.R` (regulon-score orientation checks); `lr_ko_predict.R`, `lr_ko_simulate.py` (hybrid LR-knockout consistency check, not reported) | `outputs/grn_check.csv`, `outputs/ko_results_*.csv`, orientation tables | SF.4, SF.6; Methods *In silico perturbation* |
| `dataset_summary/` | `build_dataset_summary.R` | donors, cells, meta-cells, candidate LR pairs, retained TFs and predicted interactions per comparison × cell type × condition | ST.11(b); Supplementary Data 4 |
| `supplementary_assets/` | `build_tables.R` | assembles ST.1–ST.11 CSVs from the outputs above (`ST_*.csv`) | ST.1–ST.11 |
| `figures/` | `fig_rf_validation_supp*.R`, `fig_ablation_heatmap.R`, `fig_celloracle_ko*.R` | SF.3, SF.5, SF.4/SF.6 panels | SF.3–SF.6 |

Meta-cell sequencing depths (ST.1) are column sums of the archived `pseudobulk_seurat.rds`
count matrices per comparison and cell type.

# Builds the round-2 supplementary tables from analysis outputs (24 Aug 2026).
suppressMessages({library(dplyr); library(tidyr)})
out <- "supplementary_assets"
# (a) per-cluster LR–TF model validation summary — both COVID datasets
v <- read.csv("rf_performance/results/real_covid_per_tf_metrics.csv") %>% mutate(dataset="Vaccination")
s <- read.csv("rf_performance/results/real_severity_per_tf_metrics.csv") %>% mutate(dataset="Severity")
val <- bind_rows(v,s) %>% filter(method=="rf") %>% group_by(dataset, cluster) %>%
  summarize(n_models=n(), n_obs=first(n_obs), n_predictors=first(n_pred),
            median_oob_r2=round(median(oob_r2),2), median_cv_r2=round(median(cv_r2),2),
            median_lodo_r2=round(median(lodo_r2),2), pct_lodo_positive=round(100*mean(lodo_r2>0),0),
            median_imp_stability=round(median(imp_stability_spearman),2),
            n_beat_null=sum(beats_null), .groups="drop") %>% arrange(dataset, desc(median_lodo_r2))
write.csv(val, file.path(out,"ST_validation_per_cluster.csv"), row.names=FALSE)
# (b) ntree sweep + ranking concordance
sw <- read.csv("rf_performance/results/real_covid_ntree_sweep.csv") %>% group_by(ntree) %>%
  summarize(n_models=n(), median_oob_r2=round(median(oob_r2),3), median_imp_stability=round(median(imp_stability_spearman),3), .groups="drop")
write.csv(sw, file.path(out,"ST_ntree_sweep.csv"), row.names=FALSE)
cc <- read.csv("rf_performance/results/ntree_ranking_concordance.csv") %>% group_by(dataset) %>%
  summarize(n_models=n(), median_rho_100_vs_1000=round(median(rho_100_vs_1000),3), min_rho_100_vs_1000=round(min(rho_100_vs_1000),3),
            median_rho_100_vs_100=round(median(rho_100_vs_100),3), pct_same_top1=round(100*mean(top1_same_100_1000),1),
            pct_top5_ge4of5=round(100*mean(top5_jac_100_1000>=0.66),1), .groups="drop")
write.csv(cc, file.path(out,"ST_ntree_concordance.csv"), row.names=FALSE)
# (c) negative-importance summary
neg <- bind_rows(lapply(c(Vaccination="aws_data/covid_data", Severity="aws_data/SevCOVID_Azimuthl2_data"), function(dd) {
  x <- readRDS(file.path(dd,"decipher_scores_by_regulon_and_cluster.rds")); d <- bind_rows(lapply(x, function(cl) if (is.data.frame(cl)) cl else bind_rows(cl)))
  n <- d$imp.perm<0
  data.frame(n_values=nrow(d), n_negative=sum(n), pct_negative=round(100*mean(n),1),
             median_abs_negative=signif(median(abs(d$imp.perm[n])),2), median_positive=signif(median(d$imp.perm[!n & d$imp.perm>0]),2),
             negative_mass_pct=round(100*sum(abs(d$imp.perm[n]))/sum(d$imp.perm[d$imp.perm>0]),1),
             n_negative_in_top5=sum(d$perm.rank<5 & n)) }), .id="dataset")
write.csv(neg, file.path(out,"ST_negative_importance.csv"), row.names=FALSE)
# (d) ablation AUROC matrix (threshold 2)
a <- bind_rows(read.csv("ablations/results/ablation_auroc.csv"),
               read.csv("ablations/results/ablation_auroc_ols.csv") %>% filter(method=="OLSRegr"),
               read.csv("ablations/results/ablation_score_defs.csv") %>% filter(method!="Decipher")) %>% filter(threshold==2)
ord <- c("Decipher","SpearmanCorr","LROnly","RidgeRegr","OLSRegr","NoSpearmanSign","UnsignedMagnitude","ImpOnly","SumOverTFs")
w <- a %>% mutate(method=factor(method,ord)) %>% select(dataset, method, auroc) %>% pivot_wider(names_from=method, values_from=auroc)
rank_rows <- a %>% group_by(dataset) %>% mutate(rk=rank(-auroc)) %>% group_by(method) %>% summarize(v=round(mean(rk),2)) %>% mutate(method=factor(method,ord)) %>% arrange(method)
mean_rows <- a %>% group_by(method) %>% summarize(v=round(mean(auroc),3)) %>% mutate(method=factor(method,ord)) %>% arrange(method)
w2 <- bind_rows(w %>% mutate(across(-dataset, ~round(.,3))),
                tibble(dataset="Mean AUROC") %>% bind_cols(as_tibble(setNames(as.list(mean_rows$v), as.character(mean_rows$method)))),
                tibble(dataset="Mean rank") %>% bind_cols(as_tibble(setNames(as.list(rank_rows$v), as.character(rank_rows$method)))))
write.csv(w2, file.path(out,"ST_ablation_auroc.csv"), row.names=FALSE)
# (e) donor composition — vaccination + severity
sanitize_cl <- function(x) { x <- gsub("\\+", "_plus_", x); x <- gsub("-", "_minus_", x); x <- gsub(" ", "_", x); gsub("_+", "_", x) }
norm <- function(d) { names(d)[names(d)=="effective_n_donors"] <- "eff_n_donors"
  names(d)[names(d)=="n_donors_present"] <- "n_donors"; names(d)[names(d)=="null_frac_seed_donor"] <- "expected_frac_seed_donor"
  names(d)[names(d)=="out_cluster"] <- "cluster"; d }
don <- bind_rows(lapply(c(Vaccination="donor_purity", Severity="donor_purity/severity"), function(dd) {
  cnt <- norm(read.csv(file.path(dd,"counts_donors_cells_metacells.csv")))
  skw <- norm(read.csv(file.path(dd,"donor_composition_skew.csv"))) %>% mutate(cluster = sanitize_cl(cluster)) %>% select(cluster, condition, max_donor_share, eff_n_donors)
  pur <- norm(read.csv(file.path(dd,"observed_vs_null_purity.csv"))) %>% select(cluster, condition, mean_frac_seed_donor, expected_frac_seed_donor, enrichment)
  cnt %>% select(cluster, condition, n_donors, n_cells, n_metacells) %>%
    left_join(skw, by=c("cluster","condition")) %>% left_join(pur, by=c("cluster","condition")) %>%
    mutate(across(c(max_donor_share, eff_n_donors, mean_frac_seed_donor, expected_frac_seed_donor, enrichment), ~round(.,2))) }), .id="dataset")
write.csv(don, file.path(out,"ST_donor_composition.csv"), row.names=FALSE)
# (f) alt-reference
ar <- read.csv("c8_alt_reference/outputs/reference_comparison_summary.csv") %>%
  mutate(variant=recode(variant, CD16="CD16+ monocytes", cDC2="cDC2"), lig_spearman=round(lig_spearman,2))
write.csv(ar, file.path(out,"ST_c8_alt_reference.csv"), row.names=FALSE)
cat("built:", paste(list.files(out, pattern="^ST_"), collapse=", "), "\n")

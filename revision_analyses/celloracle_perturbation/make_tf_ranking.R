# Reconstructs inputs/tf_ranking_*.csv (the TF panels' ranking) from the archived
# per-regulon Decipher score tables. Rank = mean |decipher_score| over a TF's LR–TF
# edges (decipher_score = imp.perm × sign(spearman) × regulon.val). Verified 23 Aug 2026
# to reproduce the files used by ko_cd16mono.py / ko_cd16_vax.py to 1e-15.
suppressMessages(library(dplyr))
rank_tfs <- function(store, cluster, grn_csv) {
  x <- readRDS(file.path(store, "decipher_scores_by_regulon_and_cluster.rds"))[[cluster]]
  if (!is.data.frame(x)) x <- do.call(rbind, x)
  grn_tfs <- unique(read.csv(grn_csv, row.names = 1)$source)
  x %>% group_by(regulon) %>%
    summarize(regulon_val = first(regulon.val),
              mean_abs_score = mean(abs(decipher_score)),
              max_abs_score = max(abs(decipher_score)),
              top_edge = interaction[which.max(abs(decipher_score))], .groups = "drop") %>%
    mutate(in_grn = regulon %in% grn_tfs) %>% arrange(desc(mean_abs_score))
}
check <- function(new, file) {
  ref <- read.csv(file); m <- inner_join(ref, new, by = "regulon", suffix = c(".f", ".n"))
  cat(sprintf("%s: %d TFs, max |Δ mean_abs| = %.1e, top-10 order identical = %s\n", basename(file), nrow(m),
              max(abs(m$mean_abs_score.f - m$mean_abs_score.n)), identical(head(ref$regulon, 10), head(new$regulon, 10))))
}
sev <- rank_tfs("../aws_data/SevCOVID_Azimuthl2_data", "CD16_Mono", "inputs/CD16_Mono.csv")
vax <- rank_tfs("../aws_data/covid_data", "CD16_plus_monocytes", "inputs/CD16_plus_monocytes.csv")
check(sev, "inputs/tf_ranking_cd16mono.csv"); check(vax, "inputs/tf_ranking_cd16_vax.csv")

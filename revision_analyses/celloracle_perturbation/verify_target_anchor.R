# Empirical verification that shipped regulon scores follow the pipeline's
# target-anchored orientation convention (adjustScoresOrientation /
# intracellular_scoring.R:282): for each significant real regulon, correlate
# its score with the mean log1p expression of its (capped) target genes
# across the same meta-cells. Expect overwhelmingly positive; negatives are
# orientation anomalies relevant to the manuscript's directional claims.
suppressPackageStartupMessages({library(Seurat); library(dplyr)})
local_root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/aws_data"
out_csv <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/celloracle_perturbation/outputs/target_anchor_verification.csv"

rows <- list()
for (key in c("covid", "SevCOVID_Azimuthl2", "MilCOVID_Azimuthl2")) {
  ddir <- file.path(local_root, paste0(key, "_data"))
  reg_all <- readRDS(file.path(ddir, "regulon_scores_by_cluster.rds"))
  sig_all <- readRDS(file.path(ddir, "significant_regulons_by_cluster.rds"))
  cap_all <- readRDS(file.path(ddir, "capped_regulons_all_clusters.rds"))
  pb <- readRDS(file.path(ddir, "pseudobulk_seurat.rds"))
  expr <- GetAssayData(pb, assay = "RNA", layer = "data")
  for (cl in names(sig_all)) {
    rs <- reg_all[[cl]]; cap <- cap_all[[cl]]
    cols <- intersect(colnames(rs), colnames(expr))
    for (i in seq_len(nrow(sig_all[[cl]]))) {
      tf <- sig_all[[cl]]$name[i]
      targets <- intersect(cap$target[cap$source == tf], rownames(expr))
      if (length(targets) < 5 || !tf %in% rownames(rs)) next
      tmean <- Matrix::colMeans(log1p(expr[targets, cols, drop = FALSE]))
      r <- suppressWarnings(cor(as.numeric(rs[tf, cols]), tmean, method = "spearman"))
      rows[[length(rows)+1]] <- data.frame(dataset = key, cluster = cl, tf = tf,
        deltaPagoda = sig_all[[cl]]$deltaPagoda[i], n_targets = length(targets),
        cor_score_targets = r)
    }
  }
  message(key, " done")
  rm(pb, expr); gc(verbose = FALSE)
}
res <- bind_rows(rows)
write.csv(res, out_csv, row.names = FALSE)
cat("regulons checked:", nrow(res), "\n")
cat("positive target anchor:", sum(res$cor_score_targets > 0, na.rm=TRUE),
    sprintf("(%.1f%%)\n", 100*mean(res$cor_score_targets > 0, na.rm=TRUE)))
cat("anomalies (cor < 0):", sum(res$cor_score_targets < 0, na.rm=TRUE),
    "| strong anomalies (cor < -0.2):", sum(res$cor_score_targets < -0.2, na.rm=TRUE), "\n")
print(res %>% filter(cor_score_targets < -0.2) %>% arrange(cor_score_targets) %>% head(20))

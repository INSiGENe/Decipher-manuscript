# Orientation-anchoring of PAGODA regulon signs — ALL benchmark datasets, all
# clusters. For each significant regulon: Spearman(regulon score, log1p TF
# expression) across that cluster's meta-cells; anchored_delta = deltaPagoda *
# sign(cor). Output: one long CSV for the manuscript-direction audit.

suppressPackageStartupMessages({library(Seurat); library(dplyr)})
local_root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/aws_data"
drive_root <- "/Volumes/MegaEdgar/aws_pull_20260813/results"
out_csv <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/celloracle_perturbation/outputs/regulon_orientation_anchor_all.csv"

DATASETS <- c("covid"="covid","SevCOVID_Azimuthl2"="SevCOVID_Azimuthl2",
  "MilCOVID_Azimuthl2"="MilCOVID_Azimuthl2","5yr_pic"="5yr_pic","bcg"="BCG",
  "cord_pic"="cord_pic","erp"="ERP","lupus"="lupus","sepsis"="sepsis","tnbc"="TNBC",
  "cz_influenza"="cz_influenza","cz_hpap_t1d_islets"="cz_hpap_t1d_islets",
  "cz_hnscc_hpv"="cz_hnscc_hpv","cz_human_kidney_v1.5"="cz_human_kidney_v1.5",
  "cz_cf_bronchial_biopsy"="cz_cf_bronchial_biopsy")

all_rows <- list()
for (key in names(DATASETS)) {
  ddir <- file.path(local_root, paste0(DATASETS[[key]], "_data"))
  if (!dir.exists(ddir)) ddir <- file.path(drive_root, DATASETS[[key]], "data")
  ok <- tryCatch({
    reg_all <- readRDS(file.path(ddir, "regulon_scores_by_cluster.rds"))
    sig_all <- readRDS(file.path(ddir, "significant_regulons_by_cluster.rds"))
    pb <- readRDS(file.path(ddir, "pseudobulk_seurat.rds"))
    expr_all <- GetAssayData(pb, assay = "RNA", layer = "data")
    TRUE
  }, error = function(e) { message("SKIP ", key, ": ", conditionMessage(e)); FALSE })
  if (!ok) next
  for (cl in names(sig_all)) {
    rs <- reg_all[[cl]]; sg <- sig_all[[cl]]
    cols <- intersect(colnames(rs), colnames(expr_all))
    for (i in seq_len(nrow(sg))) {
      tf <- sg$name[i]
      if (!tf %in% rownames(expr_all) || !tf %in% rownames(rs)) {
        all_rows[[length(all_rows)+1]] <- data.frame(dataset=key, cluster=cl, tf=tf,
          deltaPagoda=sg$deltaPagoda[i], cor_score_expr=NA, anchor_sign=NA,
          anchored_delta=NA, flipped=NA, weak_anchor=NA); next
      }
      r <- suppressWarnings(cor(as.numeric(rs[tf, cols]),
                                log1p(as.numeric(expr_all[tf, cols])), method="spearman"))
      s <- sign(r)
      all_rows[[length(all_rows)+1]] <- data.frame(dataset=key, cluster=cl, tf=tf,
        deltaPagoda=sg$deltaPagoda[i], cor_score_expr=r, anchor_sign=s,
        anchored_delta=sg$deltaPagoda[i]*s, flipped=!is.na(s) && s < 0,
        weak_anchor=!is.na(r) && abs(r) < 0.2)
    }
  }
  message(key, " done (", sum(sapply(all_rows, function(x) x$dataset[1]==key)), " regulons)")
  rm(pb, expr_all); gc(verbose = FALSE)
}
res <- bind_rows(all_rows)
write.csv(res, out_csv, row.names = FALSE)
message("TOTAL: ", nrow(res), " regulon x cluster rows | flipped: ",
        sum(res$flipped, na.rm=TRUE), " (", round(100*mean(res$flipped, na.rm=TRUE),1),
        "%) | weak anchors: ", sum(res$weak_anchor, na.rm=TRUE))
message("ANCHOR ALL DONE")

# Orientation-anchoring of PAGODA regulon signs — CD16_Mono, SevCOVID_Azimuthl2.
# PAGODA2 regulon-score orientation is arbitrary per regulon; anchor it by
# correlating each regulon's score with its TF's own expression across the
# same meta-cells. anchored_delta = deltaPagoda * sign(cor). Used to (a) give
# the KO panel a principled per-TF direction prediction and (b) audit any
# manuscript statement of TF activity direction in this cluster.

suppressPackageStartupMessages({library(Seurat); library(dplyr)})
base <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/aws_data/SevCOVID_Azimuthl2_data"
out  <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/celloracle_perturbation/outputs"

reg_scores <- readRDS(file.path(base, "regulon_scores_by_cluster.rds"))[["CD16_Mono"]]
sig <- readRDS(file.path(base, "significant_regulons_by_cluster.rds"))[["CD16_Mono"]]
pb <- readRDS(file.path(base, "pseudobulk_seurat.rds"))
expr <- GetAssayData(pb, assay = "RNA", layer = "data")[, colnames(reg_scores)]

anchor <- do.call(rbind, lapply(sig$name, function(tf) {
  if (!tf %in% rownames(expr) || !tf %in% rownames(reg_scores)) {
    return(data.frame(tf = tf, cor_score_expr = NA, anchor_sign = NA))
  }
  e <- as.numeric(expr[tf, ]); s <- as.numeric(reg_scores[tf, ])
  r <- suppressWarnings(cor(s, log1p(e), method = "spearman"))
  data.frame(tf = tf, cor_score_expr = r, anchor_sign = sign(r))
}))
anchor <- anchor %>%
  left_join(sig %>% select(name, deltaPagoda), by = c("tf" = "name")) %>%
  mutate(anchored_delta = deltaPagoda * anchor_sign,
         flipped = !is.na(anchor_sign) & anchor_sign < 0,
         weak_anchor = !is.na(cor_score_expr) & abs(cor_score_expr) < 0.2)
write.csv(anchor, file.path(out, "regulon_orientation_anchor_cd16mono.csv"), row.names = FALSE)

cat("anchored:", sum(!is.na(anchor$anchor_sign)), "/", nrow(anchor),
    "| flipped orientation:", sum(anchor$flipped),
    "| weak anchors (|rho|<0.2):", sum(anchor$weak_anchor), "\n\n")

# re-test KO panel signs
ko <- read.csv(file.path(out, "ko_results_cd16mono.csv")) %>%
  left_join(anchor %>% select(tf, cor_score_expr, anchored_delta, weak_anchor), by = "tf") %>%
  mutate(predicted_sign_anchored = -sign(anchored_delta),
         sign_match_anchored = predicted_sign_anchored == observed_sign)
write.csv(ko, file.path(out, "ko_results_cd16mono_anchored.csv"), row.names = FALSE)

for (g in c("top", "bottom")) {
  k <- ko %>% filter(group == g, !is.na(sign_match_anchored))
  cat(g, "panel anchored sign match:", sum(k$sign_match_anchored), "/", nrow(k),
      " (naive was:", sum(k$sign_match == "True"), ")\n")
  ks <- k %>% filter(!weak_anchor)
  cat("   strong anchors only:", sum(ks$sign_match_anchored), "/", nrow(ks), "\n")
}
print(ko %>% filter(group == "top") %>%
  select(tf, regulon_val, cor_score_expr, anchored_delta, cosine_mean_delta_vs_axis,
         predicted_sign_anchored, observed_sign, sign_match_anchored) %>%
  as.data.frame(), digits = 3)

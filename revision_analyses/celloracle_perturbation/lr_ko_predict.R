# Hybrid LR perturbation, step 1 (23 Aug 2026): in the LR–TF random forests for severity
# CD16+ monocytes, set the signalling potential of a ligand's / receptor's pairs to zero
# and read off each TF's predicted activity change; map to an expression direction via
# the orientation anchors. Output: outputs/lr_ko_tf_deltas.csv (one row per KO x TF).
suppressMessages({library(randomForest); library(dplyr)})
set.seed(123)
store <- "../aws_data/SevCOVID_Azimuthl2_data"; cl <- "CD16_Mono"
pot <- readRDS(file.path(store, "interaction_potentials_matrix_clusters_all_clusters.rds"))[[cl]]
reg <- readRDS(file.path(store, "regulon_scores_by_cluster.rds"))[[cl]]
sig <- readRDS(file.path(store, "significant_regulons_by_cluster.rds"))[[cl]]; tfs <- sig$name[sig$class == "real"]
pb <- readRDS(file.path(store, "pseudobulk_seurat.rds")); meta <- attr(pb, "meta.data")[colnames(pot), ]; rm(pb)
sc <- readRDS(file.path(store, "decipher_scores_by_cluster.rds"))[[cl]] %>% arrange(desc(abs(decipher_score)))
anchor <- read.csv("outputs/regulon_orientation_anchor_cd16mono.csv")
x <- t(pot); colnames(x) <- rownames(pot); y <- reg[tfs, colnames(pot), drop = FALSE]
cat("meta-cells:", nrow(x), " pairs:", ncol(x), " TFs:", nrow(y), " conditions:", paste(names(table(meta$condition)), table(meta$condition)), "\n")
# fit one RF per TF (manuscript settings) and cache predictions at baseline
fits <- lapply(tfs, function(tf) { set.seed(123); randomForest(x = x, y = as.numeric(y[tf, ]), ntree = 100) }); names(fits) <- tfs
base <- sapply(fits, function(f) predict(f, x))            # meta-cells x TFs
sdy  <- apply(y, 1, sd)
# KO panels: top 5 / bottom 5 / random 5 LR pairs by |Decipher score|
n <- nrow(sc); pairs <- bind_rows(sc[1:5, ] %>% mutate(panel = "top"), sc[(n-4):n, ] %>% mutate(panel = "bottom"),
                                  sc[sample(6:(n-5), 5), ] %>% mutate(panel = "random"))
ko_one <- function(cols, label, panel, kind, anchor_pair) {
  xk <- x; xk[, cols] <- 0
  pk <- sapply(fits, function(f) predict(f, xk))
  dz <- colMeans(pk - base) / sdy                                     # SD units, all meta-cells
  dz_sev <- colMeans((pk - base)[meta$condition != "control", , drop = FALSE]) / sdy
  data.frame(ko = label, panel = panel, kind = kind, anchor_pair = anchor_pair, n_pairs_zeroed = length(cols),
             tf = tfs, dz = as.numeric(dz), dz_severe = as.numeric(dz_sev), stringsAsFactors = FALSE)
}
out <- list()
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i, ]
  rc <- sc$interaction[sc$receptor == p$receptor]; lc <- sc$interaction[sc$ligand == p$ligand]
  out[[length(out)+1]] <- ko_one(rc, paste0("receptorKO:", p$receptor), p$panel, "receptor", p$interaction)
  out[[length(out)+1]] <- ko_one(lc, paste0("ligandKO:", p$ligand), p$panel, "ligand", p$interaction)
  out[[length(out)+1]] <- ko_one(p$interaction, paste0("pairKO:", p$interaction), p$panel, "pair", p$interaction)   # 23 Aug: single-interaction KO
}
res <- bind_rows(out) %>% left_join(anchor %>% select(tf, anchor_sign, weak_anchor), by = "tf") %>%
  mutate(expr_dir = sign(dz) * anchor_sign)
write.csv(res, "outputs/lr_ko_tf_deltas.csv", row.names = FALSE)
write.csv(pairs %>% select(panel, interaction, ligand, receptor, decipher_score), "outputs/lr_ko_pairs.csv", row.names = FALSE)
cat("written; KOs:", length(unique(res$ko)), "| per-KO max |dz| range:", paste(round(range(tapply(abs(res$dz), res$ko, max)), 2), collapse="-"), "\n")
print(res %>% group_by(ko, panel, kind) %>% summarize(n_pairs = first(n_pairs_zeroed), top_tf = tf[which.max(abs(dz))], max_abs_dz = round(max(abs(dz)), 2), n_dz_gt_0.25 = sum(abs(dz) > 0.25), .groups = "drop") %>% arrange(panel, kind), n = 40, width = 120)

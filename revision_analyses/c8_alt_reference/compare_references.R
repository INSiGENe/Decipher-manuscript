# Compare C8 Decipher results across control-arm references (R2.5 robustness).
# Baseline = CD14+ monocytes (manuscript choice), variants = CD16+ mono, cDC2.
# Ligand-level comparisons are used where representative-LR selection makes
# interaction-level comparison unstable.
suppressMessages(library(dplyr))
C8 <- "CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells"
out_root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/c8_alt_reference/outputs"

load_variant <- function(v) {
  list(scores = readRDS(file.path(out_root, v, "decipher_scores_by_cluster.rds"))[[C8]],
       sig    = readRDS(file.path(out_root, v, "significant_regulons_by_cluster.rds"))[[C8]])
}
base <- load_variant("CD14_plus_monocytes")
variants <- list(CD16 = load_variant("CD16_plus_monocytes"), cDC2 = load_variant("cDC2"))

lig_scores <- function(sc) sc %>% group_by(ligand) %>%
  filter(abs(decipher_score) == max(abs(decipher_score))) %>% ungroup() %>%
  select(ligand, decipher_score) %>% arrange(desc(decipher_score))

base_lig <- lig_scores(base$scores)
narrative <- c("SLAMF7","CCL2","SERPING1","IFNG","MIF","ICAM4","HGF")

res <- list()
for (v in names(variants)) {
  vd <- variants[[v]]
  vl <- lig_scores(vd$scores)
  m_int <- inner_join(base$scores %>% select(interaction, b = decipher_score),
                      vd$scores %>% select(interaction, s = decipher_score), by = "interaction")
  m_lig <- inner_join(base_lig %>% rename(b = decipher_score),
                      vl %>% rename(s = decipher_score), by = "ligand")
  top10_b <- head(base_lig$ligand, 10); top10_v <- head(vl$ligand, 10)
  sig_overlap <- length(intersect(base$sig$name, vd$sig$name))
  cat(sprintf("\n=== %s vs baseline ===\n", v))
  cat(sprintf("interactions: %d shared | score spearman %.3f\n", nrow(m_int), cor(m_int$b, m_int$s, method="spearman")))
  cat(sprintf("ligands: %d shared | score spearman %.3f | top-10 ligand overlap %d/10\n",
      nrow(m_lig), cor(m_lig$b, m_lig$s, method="spearman"), length(intersect(top10_b, top10_v))))
  cat(sprintf("significant TFs: %d (base %d) | overlap %d\n", nrow(vd$sig), nrow(base$sig), sig_overlap))
  cat("variant top-10 ligands:", paste(top10_v, collapse=", "), "\n")
  for (lg in narrative) {
    rb <- which(base_lig$ligand == lg); rv <- which(vl$ligand == lg)
    cat(sprintf("  %-9s base rank %s -> %s rank %s\n", lg,
        ifelse(length(rb), rb, "absent"), v, ifelse(length(rv), rv, "absent")))
  }
  res[[v]] <- data.frame(variant = v, n_shared_lig = nrow(m_lig),
    lig_spearman = cor(m_lig$b, m_lig$s, method="spearman"),
    top10_overlap = length(intersect(top10_b, top10_v)),
    sig_tf_overlap = sig_overlap, n_sig_tf = nrow(vd$sig))
}
cat("\nbaseline top-10 ligands:", paste(head(base_lig$ligand,10), collapse=", "), "\n")
write.csv(bind_rows(res), file.path(out_root, "reference_comparison_summary.csv"), row.names = FALSE)

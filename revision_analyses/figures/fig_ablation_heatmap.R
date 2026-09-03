# Fig 2 new panel (numbering TBD) — ablation AUROC heatmap on the CytoSig benchmark.
# Standalone; does not touch figure_2.R. Inputs: analysis_r2/ablations/results/*.csv
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(viridisLite)})
root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2"
a1 <- read.csv(file.path(root, "ablations/results/ablation_auroc.csv"))
a2 <- read.csv(file.path(root, "ablations/results/ablation_auroc_ols.csv")) %>% filter(method == "OLSRegr")
a3 <- read.csv(file.path(root, "ablations/results/ablation_score_defs.csv")) %>% filter(method != "Decipher")
d <- bind_rows(a1, a2, a3) %>% filter(threshold == 2)

variant_levels <- c("Decipher","SpearmanCorr","LROnly","RidgeRegr","OLSRegr",
                    "NoSpearmanSign","UnsignedMagnitude","ImpOnly","SumOverTFs")
variant_labels <- c("Full\nframework","Spearman\ncorrelation","LR potential\nonly","Ridge\nregression",
                    "OLS\nregression","No sign","Unsigned\nmagnitude","Importance\nonly","Sum over\nTFs")
ds_levels <- c("covid","5yr_pic","bcg","cord_pic","erp","tnbc","cz_influenza","cz_hnscc_hpv",
               "cz_cf_bronchial_biopsy","SevCOVID_Azimuthl2","MilCOVID_Azimuthl2")
ds_labels <- c("COVID-19 vaccination","Poly-IC (adult)","BCG vaccination","Poly-IC (cord)",
               "Breast cancer ICB (ER+)","Breast cancer ICB (TNBC)","Influenza","HNSCC",
               "Cystic fibrosis","Severe COVID-19","Mild COVID-19")

d <- d %>% filter(method %in% variant_levels, dataset %in% ds_levels) %>%
  mutate(method = factor(method, variant_levels, variant_labels),
         dataset = factor(dataset, ds_levels, ds_labels))
ntrue <- d %>% filter(method == "Full\nframework") %>% select(dataset, n_true)

# summary rows
ranks <- d %>% group_by(dataset) %>% mutate(rank = rank(-auroc, ties.method = "average")) %>% ungroup()
summ <- bind_rows(
  d %>% group_by(method) %>% summarize(value = mean(auroc), .groups = "drop") %>% mutate(row = "Mean AUROC"),
  ranks %>% group_by(method) %>% summarize(value = mean(rank), .groups = "drop") %>% mutate(row = "Mean rank"))
summ$row <- factor(summ$row, c("Mean AUROC","Mean rank"))

star <- d %>% left_join(ntrue, by = "dataset", suffix = c("", ".ref")) %>%
  mutate(lab = sprintf("%.2f%s", auroc, ifelse(n_true.ref < 10, "*", "")))

p_main <- ggplot(star, aes(method, dataset, fill = auroc)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab, colour = auroc > 0.60), size = 3.1) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  scale_fill_viridis_c(option = "D", name = "AUROC") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 9.5, lineheight = 0.85),
        legend.position = "right", legend.key.height = unit(1.2, "cm"), plot.margin = margin(4, 4, 0, 4))

p_sum <- ggplot(summ, aes(method, row)) +
  geom_tile(fill = "grey93", colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(row == "Mean rank", sprintf("%.2f", value), sprintf("%.3f", value)),
                fontface = ifelse(method == "Full\nframework", "bold", "plain")), size = 3.1) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0, 4, 4, 4))

library(patchwork)
p <- (p_main / p_sum) + plot_layout(heights = c(11, 2.2)) +
  plot_annotation(caption = "* fewer than 10 CytoSig-active ligands; AUROC at median z-score > 2")
ggsave("fig_ablation_heatmap.png", p, width = 10.5, height = 6.6, dpi = 200, bg = "white")
ggsave("fig_ablation_heatmap.pdf", p, width = 10.5, height = 6.6)
cat("written\n")

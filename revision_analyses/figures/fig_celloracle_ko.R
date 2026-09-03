# R1.1 — CellOracle in-silico TF knockout panels (severity CD16+ mono = main Fig 5 panel;
# vaccination CD16+ mono = supplementary). Standalone.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})
root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/celloracle_perturbation/outputs"
th <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
grp_lab <- c(top = "Decipher top-10 TFs", bottom = "Lowest-ranked TFs", random = "Non-prioritized TFs")
make_panels <- function(csv, axis_lab, toward_lab, away_lab, tag_a, tag_b) {
  d <- read.csv(file.path(root, csv)) %>%
    mutate(group = factor(grp_lab[group], grp_lab),
           cos = cosine_mean_delta_vs_axis, mag = shift_magnitude,
           dir = ifelse(cos < 0, toward_lab, away_lab)) %>%
    arrange(group, cos) %>% mutate(tf = factor(tf, tf))
  pa <- ggplot(d, aes(cos, tf, colour = dir)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_segment(aes(x = 0, xend = cos, yend = tf), linewidth = 0.6) +
    geom_point(aes(size = mag)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = setNames(c("#21918c", "#d95f02"), c(toward_lab, away_lab)), name = "Predicted shift") +
    scale_size_continuous(range = c(1.5, 6), name = "Shift magnitude") +
    labs(x = paste0("Cosine alignment with ", axis_lab), y = NULL, title = tag_a) + th +
    theme(strip.text.y = element_text(angle = 0, hjust = 0), legend.position = "bottom", legend.box = "vertical")
  top <- d %>% filter(group == grp_lab["top"]); ctrl <- d %>% filter(group != grp_lab["top"])
  p <- wilcox.test(abs(top$cos), abs(ctrl$cos))$p.value
  dd <- d %>% mutate(g2 = ifelse(group == grp_lab["top"], "Decipher\ntop-10", "Control\npanels"))
  pb <- ggplot(dd, aes(g2, abs(cos), fill = g2)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.8) +
    scale_fill_manual(values = c("#fde725", "#bbbbbb"), guide = "none") +
    annotate("text", x = 1.5, y = max(abs(dd$cos)) * 1.08, label = sprintf("Wilcoxon p = %.3g", p), size = 3.6) +
    labs(x = NULL, y = "|cosine alignment|", title = tag_b) + th + theme(panel.grid.major.y = element_line())
  list(pa = pa, pb = pb, p = p, d = d)
}
sev <- make_panels("ko_results_cd16mono_anchored.csv", "healthy→severe axis", "toward healthy", "toward severe", "a  CD16+ monocytes, severity", "b")
vax <- make_panels("ko_results_cd16_vax.csv", "baseline→post-vaccination axis", "toward baseline", "toward post-vaccination", "a  CD16+ monocytes, vaccination", "b")
ggsave("fig_celloracle_ko_severity.png", (sev$pa | sev$pb) + plot_layout(widths = c(2.2, 1)), width = 10.5, height = 6.2, dpi = 200, bg = "white")
ggsave("fig_celloracle_ko_severity.pdf", (sev$pa | sev$pb) + plot_layout(widths = c(2.2, 1)), width = 10.5, height = 6.2)
ggsave("fig_celloracle_ko_vaccination_supp.png", (vax$pa | vax$pb) + plot_layout(widths = c(2.2, 1)), width = 10.5, height = 6.2, dpi = 200, bg = "white")
ggsave("fig_celloracle_ko_vaccination_supp.pdf", (vax$pa | vax$pb) + plot_layout(widths = c(2.2, 1)), width = 10.5, height = 6.2)
cat(sprintf("severity p = %.4f | vaccination p = %.3f\n", sev$p, vax$p))

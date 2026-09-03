# Supplementary figure — LR–TF model validation (R2.4). Standalone.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(viridisLite)})
root <- "/Users/edgarbasto/Documents/Decipher_revision_2/analysis_r2/rf_performance/results"
v <- read.csv(file.path(root, "real_covid_per_tf_metrics.csv")) %>% mutate(dataset = "Vaccination")
s <- read.csv(file.path(root, "real_severity_per_tf_metrics.csv")) %>% mutate(dataset = "Severity")
d <- bind_rows(v, s) %>% mutate(dataset = factor(dataset, c("Vaccination","Severity")))
conc <- read.csv(file.path(root, "ntree_ranking_concordance.csv")) %>%
  mutate(dataset = factor(ifelse(dataset == "vaccination", "Vaccination", "Severity"), c("Vaccination","Severity")))
rf <- d %>% filter(method == "rf")
pal <- c(Vaccination = "#440154", Severity = "#21918c")
th <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())

# (a) observed OOB vs null 95th percentile
pa <- ggplot(rf, aes(null_q95, oob_r2, colour = dataset)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.55, size = 1.3) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = expression("Null 95th percentile OOB"~R^2), y = expression("Observed OOB"~R^2),
       title = "a  Observed vs permutation null (one point per model)") + th +
  theme(legend.position = c(0.8, 0.2), legend.background = element_rect(fill = "white", colour = NA))

# (b) LODO by cell type
clean <- function(x) x %>% gsub("_plus_", "+ ", .) %>% gsub("_minus_", "- ", .) %>% gsub("_", " ", .) %>%
  gsub("CD14\\+ BDCA1\\+ PD- L1\\+ cells", "C8 (CD14+BDCA1+PD-L1+)", .)
lodo <- rf %>% mutate(cl = clean(cluster)) %>% group_by(dataset, cl) %>% mutate(med = median(lodo_r2)) %>% ungroup()
pb <- ggplot(lodo, aes(reorder(cl, med), lodo_r2, fill = dataset)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85, linewidth = 0.3) +
  scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~dataset, scales = "free_x") + coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, y = expression("Donor-held-out"~R^2)) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# (c) CV R² by regressor
meth_lab <- c(rf = "Random forest", ridge = "Ridge", lm = "OLS", spearman = "Univariate\n(best LR pair)")
cv <- d %>% mutate(method = factor(meth_lab[method], meth_lab))
pc <- ggplot(cv, aes(method, cv_r2, fill = dataset)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = pal, name = NULL) + coord_cartesian(ylim = c(-0.25, 1)) +
  labs(x = NULL, y = expression("5-fold CV"~R^2)) + th +
  theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = "white", colour = NA))

# (d) tree-count concordance vs noise floor
cc <- conc %>% select(dataset, rho_100_vs_1000, rho_100_vs_100) %>%
  pivot_longer(-dataset, names_to = "comparison", values_to = "rho") %>%
  mutate(comparison = factor(comparison, c("rho_100_vs_1000","rho_100_vs_100"),
                             c("100 vs 1,000 trees", "100 vs 100 trees (different seed)")))
pd <- ggplot(cc, aes(rho, fill = comparison)) +
  geom_histogram(binwidth = 0.02, alpha = 0.65, position = "identity", colour = NA) +
  scale_fill_manual(values = c("#fde725", "#31688e"), name = NULL) +
  facet_wrap(~dataset, ncol = 1, scales = "free_y") +
  labs(x = "Spearman correlation of importance rankings", y = "Models",
       title = "d  Ranking concordance across tree counts") + th +
  theme(legend.position = "bottom")

p <- (pa | pc) / (pb | pd) + plot_layout(heights = c(1, 1.25))
ggsave("fig_rf_validation_supp_bare.png", p, width = 13, height = 10, dpi = 170, bg = "white")
ggsave("fig_rf_validation_supp_bare.pdf", p, width = 13, height = 10)
cat("written\n")

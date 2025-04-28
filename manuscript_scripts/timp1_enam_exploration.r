library(Seurat)
library(ggplot2)
library(patchwork)

# Define genes of interest
genes_of_interest <- c("TIMP1", "ENAM")

# Optionally, create a new metadata column combining severity and celltype
test$severity_celltype <- paste(test$condition, test$cluster, sep = "_")

# Plot each gene as a violin plot
p_list <- lapply(genes_of_interest, function(gene) {
  VlnPlot(
    test,
    features = gene,
    group.by = "severity_celltype",
    pt.size = 0,
    cols = NULL
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(gene)
})

# Combine plots
p_combined <- wrap_plots(p_list, ncol = 1)
# Save the combined violin plot
ggsave(
  filename = "figures/timp1_enam_violin_by_severity_and_celltype.png",
  plot = p_combined,
  width = 40,
  height = 30,
  dpi = 150
)


library(Seurat)
library(ggplot2)

# Load the object
test <- readRDS("results/SevCOVID_Azimuthl2/pre_processing/seurat_object_oi.rds")

# Generate the plot
p <- VlnPlot(
  object = test,
  features = "SLAMF7",
  group.by = "predicted.celltype.l2",
  split.by = "severity_group",
  pt.size = 0.1
) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("SLAMF7_violin_by_cluster_and_severity.png", plot = p, width = 10, height = 6, dpi = 300)

test <- readRDS("results/SevCOVID_Azimuthl2/data/seurat_object_oi.rds")

test <- readRDS("results/covid/pre_processing/seurat_object_oi.rds")

# Generate the plot
p <- VlnPlot(
  object = test,
  features = "SLAMF7",
  group.by = "cluster",
  split.by = "severity_group",
  pt.size = 0.1
) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("SLAMF7_violin_by_cluster_and_severity.png", plot = p, width = 10, height = 6, dpi = 300)

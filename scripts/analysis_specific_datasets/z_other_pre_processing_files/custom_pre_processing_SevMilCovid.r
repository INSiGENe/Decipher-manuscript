#restore environment
library(renv)
renv::restore()

############
#libraries ----
############

library(Seurat)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

cleanSymbols <- function(string) {
  # Remove or replace various symbols with safer alternatives
  string <- stringr::str_remove_all(string, "/")
  string <- stringr::str_replace_all(string, " ", "_")
  string <- stringr::str_replace_all(string, "\\+", "_plus_")
  string <- stringr::str_replace_all(string, "\\-", "_minus_")
  string <- stringr::str_replace_all(string, "\\(", "")
  string <- stringr::str_replace_all(string, "\\)", "")
  string <- stringr::str_replace_all(string, "%", "_percent_")
  string <- stringr::str_replace_all(string, "\\.", "_dot_")
  string <- stringr::str_replace_all(string, ",", "_comma_")
  string <- stringr::str_replace_all(string, ":", "_colon_")
  string <- stringr::str_replace_all(string, ";", "_semicolon_")
  string <- stringr::str_replace_all(string, "&", "_and_")
  string <- stringr::str_replace_all(string, "\\?", "_question_")
  string <- stringr::str_replace_all(string, "!", "_exclamation_")
  string <- stringr::str_replace_all(string, "\"", "_quote_")
  string <- stringr::str_replace_all(string, "'", "_apostrophe_")
  string <- stringr::str_replace_all(string, "=", "_equals_")
  string <- stringr::str_replace_all(string, "\\*", "_asterisk_")
  string <- stringr::str_replace_all(string, "#", "_hash_")
  string <- stringr::str_replace_all(string, "@", "_at_")
  string <- stringr::str_replace_all(string, "\\$", "_dollar_")
  string <- stringr::str_replace_all(string, "\\^", "_caret_")
  string <- stringr::str_replace_all(string, "<", "_less_than_")
  string <- stringr::str_replace_all(string, ">", "_greater_than_")
  string <- stringr::str_replace_all(string, "\\[", "_lbracket_")
  string <- stringr::str_replace_all(string, "\\]", "_rbracket_")
  string <- stringr::str_replace_all(string, "\\{", "_lbrace_")
  string <- stringr::str_replace_all(string, "\\}", "_rbrace_")
  string <- stringr::str_replace_all(string, "\\|", "_pipe_")
  string <- stringr::str_replace_all(string, "\\\\", "_backslash_")
  string <- stringr::str_replace_all(string, "/", "_slash_")
  string <- stringr::str_replace_all(string, "__", "_")  # Remove double underscores
  return(string)
}


############
#parameters ----
############

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a dataset key, e.g. 'cz_rcc'")
}
dataset_key <- args[1]

config <- jsonlite::fromJSON(txt = "scripts/config.json")

if (!dataset_key %in% names(config)) {
  stop(paste("Dataset key not found in config:", dataset_key))
}

cfg <- config[[dataset_key]][["custom_pre_processing"]]

if (is.null(cfg)) {
  stop(paste("No 'custom_pre_processing' section found for dataset:", dataset_key))
}

############
#analysis ----
############

set.seed(123)

pre_processing_path <- cfg$input_path
seurat_object_rds_path <- file.path(pre_processing_path, "combined_seurat_for_processing.rds")
combined_seurat_for_processing <- readRDS(seurat_object_rds_path)

# Check if gene names are Ensembl IDs (rough heuristic: Ensembl IDs start with ENS)
gene_names <- rownames(combined_seurat_for_processing)
has_ensembl <- all(grepl("^ENS", gene_names))

if (has_ensembl) {
  message("Gene names appear to be Ensembl IDs, mapping to HGNC symbols...")

  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_names,
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  ) %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  # Subset matrix and rename
  common_genes <- intersect(mapping$ENSEMBL, rownames(combined_seurat_for_processing))
  combined_seurat_for_processing <- subset(combined_seurat_for_processing, features = common_genes)
  rownames(combined_seurat_for_processing) <- mapping$SYMBOL[match(rownames(combined_seurat_for_processing), mapping$ENSEMBL)]
}

combined_seurat_for_processing$severity_group <- ifelse(
  is.na(combined_seurat_for_processing$disease_severity) | combined_seurat_for_processing$disease_severity == "",
  "Healthy",
  combined_seurat_for_processing$disease_severity
)

# Extract and clean metadata
# Clean and extract metadata, marking 'case', 'control', and 'other'
combined_seurat_for_processing@meta.data <- combined_seurat_for_processing@meta.data %>%
  mutate(
    cluster = cleanSymbols(.data[[cfg$cluster_meta_field]]),
    condition = .data[[cfg$condition_meta_field]],
    sample_key = paste(.data[[cfg$individual_meta_field]], condition, sep = "_"),
    condition_original = case_when(
      condition == cfg$case_condition ~ "case",
      condition == cfg$control_condition ~ "control",
      TRUE ~ "other"  # optional: label unclassified conditions
    )
  )

  # Keep only case and control
combined_seurat_for_processing <- subset(
  combined_seurat_for_processing,
  subset = condition_original %in% c("case", "control")
)

# Save to output path
dir.create(cfg$output_path,recursive = TRUE)
saveRDS(combined_seurat_for_processing, file = file.path(cfg$output_path, "seurat_object_oi.rds"))

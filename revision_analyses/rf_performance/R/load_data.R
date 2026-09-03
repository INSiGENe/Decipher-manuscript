# Data interface: hand REAL Decipher pipeline outputs to the harness.
#
# The harness consumes, per cluster, a list with:
#   cluster    : character label
#   potentials : matrix, representative LR pairs (rows, "LIGAND-RECEPTOR")
#                x meta-cells (cols) — the filtered interaction-potential
#                matrix the RF models are trained on
#   tf_activity: matrix, significant real regulons (rows) x meta-cells
#                (cols, same order as potentials) — exactly the TF models the
#                manuscript pipeline fits (decipher_scoring.R:131-140 loops
#                significant_regulon_deltas$name; all are class "real")
#   meta       : data.frame(metacell_id, cluster, condition, seed_donor)
#   membership : data.frame(metacell_id, cell_id, donor) for the LODO
#                contamination-drop rule
#
# Source of the real data (verified 18 Aug 2026):
#   - `analysis_r2/aws_data/covid_data/` — local copy of the AWS archive's
#     `results/covid/data/` store, i.e. the objects the manuscript run produced:
#     `interaction_potentials_matrix_clusters_all_clusters.rds`,
#     `regulon_scores_by_cluster.rds`, `significant_regulons_by_cluster.rds`,
#     `pseudobulk_seurat.rds` (condition labels; meta-cell IDs are seed-cell
#     barcodes, some carrying native ".1" suffixes from sample-merge barcode
#     collisions — NEVER strip suffixes, match IDs verbatim).
#   - `analysis_r2/donor_purity/metacell_membership.csv` — regenerated
#     meta-cell membership. RECONCILED 18 Aug 2026: its 9832
#     (cluster, condition, seed) keys match the archive pseudobulk exactly,
#     so the R2.3 "reconcile before quoting" caveat is closed for vaccination.
#   - `analysis_r2/donor_purity/cell_donor_lookup.csv` — barcode -> pt_id for
#     all 59,918 cells (from the archived seurat_object_oi.rds), giving every
#     meta-cell a seed donor.
#
# Cluster quirk: CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells (C8) has no
# control arm; its matrix's 600 control columns are CD14_plus_monocytes
# control meta-cells (the pipeline's paramPairings substitution). meta$cluster
# reports each column's true cluster of origin.

DEFAULT_COVID_DATA <- file.path(dirname(dirname(getwd())), "aws_data", "covid_data")

load_decipher_dataset <- function(data_dir,
                                  membership_csv,
                                  donor_lookup_csv,
                                  clusters = NULL) {
  potentials_all <- readRDS(file.path(data_dir, "interaction_potentials_matrix_clusters_all_clusters.rds"))
  regulons_all   <- readRDS(file.path(data_dir, "regulon_scores_by_cluster.rds"))
  sig_all        <- readRDS(file.path(data_dir, "significant_regulons_by_cluster.rds"))
  pb             <- readRDS(file.path(data_dir, "pseudobulk_seurat.rds"))
  pb_meta        <- attr(pb, "meta.data")   # avoids needing Seurat loaded
  rm(pb)

  mem_raw <- read.csv(membership_csv, stringsAsFactors = FALSE)
  lookup  <- read.csv(donor_lookup_csv, stringsAsFactors = FALSE)
  # donor column name differs by dataset (pt_id for vaccination, sample_id/donor
  # for severity) — take whichever is present
  lk_donor_col  <- intersect(c("pt_id", "sample_id", "donor"), colnames(lookup))[1]
  mem_donor_col <- intersect(c("pt_id", "sample_id", "donor"), colnames(mem_raw))[1]
  stopifnot(!is.na(lk_donor_col), !is.na(mem_donor_col))
  donor_of <- setNames(as.character(lookup[[lk_donor_col]]), lookup$barcode)

  if (is.null(clusters)) clusters <- names(potentials_all)
  out <- lapply(clusters, function(cl) {
    pot <- potentials_all[[cl]]
    ids <- colnames(pot)

    stopifnot(all(ids %in% rownames(pb_meta)))
    col_meta <- pb_meta[ids, c("cluster", "condition")]

    sig_tfs <- sig_all[[cl]]$name
    stopifnot(all(sig_all[[cl]]$class == "real"),
              all(sig_tfs %in% rownames(regulons_all[[cl]])),
              identical(colnames(regulons_all[[cl]]), ids))
    tf_activity <- regulons_all[[cl]][sig_tfs, ids, drop = FALSE]

    seed_donor <- unname(donor_of[ids])
    stopifnot(!anyNA(seed_donor))

    meta <- data.frame(
      metacell_id = ids,
      cluster     = col_meta$cluster,
      condition   = col_meta$condition,
      seed_donor  = seed_donor,
      stringsAsFactors = FALSE
    )

    # membership join on (true cluster of origin, condition, seed barcode);
    # verbatim seed IDs — verified unique per stratum and 100% matched
    key_col <- paste(meta$cluster, meta$condition, meta$metacell_id)
    key_mem <- paste(mem_raw$out_cluster, mem_raw$condition, mem_raw$metacell_seed)
    mem <- mem_raw[key_mem %in% key_col, ]
    membership <- data.frame(
      metacell_id = mem$metacell_seed,
      cell_id     = mem$member_barcode,
      donor       = as.character(mem[[mem_donor_col]]),
      stringsAsFactors = FALSE
    )
    stopifnot(setequal(unique(membership$metacell_id), ids))

    list(cluster = cl, potentials = pot, tf_activity = tf_activity,
         meta = meta, membership = membership)
  })
  names(out) <- clusters
  out
}

load_decipher_cluster <- function(cluster,
                                  data_dir,
                                  membership_csv,
                                  donor_lookup_csv) {
  load_decipher_dataset(data_dir, membership_csv, donor_lookup_csv,
                        clusters = cluster)[[1]]
}

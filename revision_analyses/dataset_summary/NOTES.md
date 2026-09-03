# Dataset summary table (round-2 reviewer request) — provenance notes

Built 2 Sep 2026 by `build_dataset_summary.R` from the read-only archive copies under
`analysis_r2/aws_data/results/<comparison>/`. Outputs:

- `ST_dataset_summary.csv` — one row per comparison x analysed cell type x condition (256 rows, 15 comparisons, 128 cell types).
- `ST_dataset_summary_by_comparison.csv` — roll-up per comparison.
- `cache/<comparison>.rds` — metadata-only cache of the two Seurat objects plus the small count lists (40 MB total; git-ignored). Delete to force a full reload.

## Where each count comes from

| Column | Object | Rule |
|---|---|---|
| `cell_type` | `data/decipher_scores_by_cluster.rds` | `names()` of the list = clusters that reached the final scoring step ("analysed"). In all 15 comparisons this is identical to the clusters present in `pseudobulk_seurat.rds` and to the names of `L_set_relevant_features_all_clusters.rds`, i.e. every cluster that survived meta-cell generation was scored. |
| `condition` | `pseudobulk_seurat.rds@meta.data$condition` | Always `case` / `control` (pipeline labels). |
| `condition_column_sc` | `pre_processing/seurat_object_oi.rds@meta.data` | The single-cell object keeps the dataset-native labels in `condition` and the pipeline's `case`/`control` in `condition_original` (14 datasets) or `original_condition` (BCG). The script picks whichever column holds exactly `case`/`control`. Mapping to native labels is listed per dataset below. |
| `n_cells` | `seurat_object_oi.rds@meta.data` | Cells with `cluster == cell_type` and pipeline condition == `condition`. This is the post-QC single-cell object handed to the meta-cell module (all cells in the object, before the meta-cell input subsample). |
| `n_cells_input_to_metacells` | derived | `min(n_cells, 1200 * (k + 1))`: `generateMetaCellMatrices()` randomly subsamples each arm to at most `1200*(k+1)` cells before building meta-cells (`R/pre_processing_meta_cells.R`, Decipher-manuscript repo). |
| `n_metacells` | `pseudobulk_seurat.rds@meta.data` | Meta-cells (columns) with that `cluster` and `condition`. Rule in code: a cluster is kept iff `floor(min(arm cells)/(k+1)) >= 100` (`min_meta_cells`), and each arm then gets `min(600, floor(min(arm cells)/(k+1)))` meta-cells (both arms equal). `expected_metacells_per_arm` / `metacell_rule_ok` in the CSV check this; it holds everywhere except the covid pairing quirk (below). |
| `n_donors` | `seurat_object_oi.rds@meta.data[[donor_column]]` | Unique donor IDs among the cells counted in `n_cells` (so per cell type x condition). Roll-up gives per-arm and total unique donors across the whole object. `donor_column` names the column; NA when none exists (BCG). |
| `n_candidate_lr_pairs` | `L_set_relevant_features_all_clusters.rds[[cluster]]` | Unique `interaction` strings. This is the L.set (2,293 pairs in the database, identical for all 15 runs) filtered to ligands expressed anywhere (>=10 % of meta-cells in some cluster) and receptors expressed in the receiving cluster (`getRelevantFeaturesForEachCluster`). Per cell type; repeated on both condition rows. |
| `n_lr_pairs_diff_potential` | `interaction_deltas_by_cluster.rds[[cluster]]` | Rows = candidate pairs whose interaction potential differs between arms (`calculateInteractionDeltasAllClusters`); this is the set that enters the random-forest step. Extra column, not requested. |
| `n_retained_tfs` | `significant_regulons_by_cluster.rds[[cluster]]` | Unique `name` (TF) among regulons with a significant PAGODA delta between arms. |
| `n_predicted_interactions` | `decipher_scores_by_cluster.rds[[cluster]]` | `nrow()` of the per-cluster tibble (one row per ligand-receptor interaction with a Decipher score). `n_unique_predicted_lr` (unique `interaction`) is identical in every cluster, i.e. no duplicate interaction rows. |
| `k`, `min_meta_cells` | `data/parameter_record.csv` | k = 1 (8 runs), 2 (6 runs: BCG, covid, kidney, lupus, TNBC), 3 (sepsis). `min_meta_cells` = 100 everywhere. |

`dataset_label` is a hand-written human-readable label (see config block in the script); the
dataset attributions in parentheses are my best reading of the metadata (GEO/CELLxGENE IDs,
patient-ID prefixes) and should be checked against the manuscript's Methods before use.

## Per-dataset donor / condition situation

| comparison | donor column | donors (case / control / total) | native labels (case = ...) | notes |
|---|---|---|---|---|
| 5yr_pic | `Sample` (Donor1, Donor2) | 2 / 2 / 2 | `condition`: PIC = case, CTRL = control | Paired design: both donors contribute to both arms. `orig.ident`/`Group` = donor x stimulus (4 samples). |
| BCG | none | NA | `condition`: D21 = case, D0 = control (`original_condition` holds case/control) | Only `orig.ident` = `GSE232186_RAW` and `timepoint`. No donor, sample or animal ID survives in the object. Lung tissue cell types (AT2, alveolar macrophages, ciliated). |
| cord_pic | `Sample` (Donor1, Donor2) | 2 / 2 / 2 | `condition`: PIC = case, CTRL = control | Paired, as 5yr_pic. |
| covid | `pt_id` | 6 / 6 / 6 | `condition`: day 22 = case, day 0 = control | Paired (`sample_id` = 12 = 6 donors x 2 days). Reproduces the verified anchor (cluster B: 2,373 case / 2,901 control cells, 600 + 600 meta-cells, k = 2). |
| cz_cf_bronchial_biopsy | `donor_id` | 3 / 19 / 22 | `condition`: "cystic fibrosis" = case, "normal" = control | Strongly unbalanced: CF arm is 3 donors from one lab/publication (Berg 2025, 10x v3), controls are 19 donors from three other publications (10x v2). `sample_key` = donor_condition (22, same as donors). |
| cz_hnscc_hpv | `donor_id` | 10 / 7 / 10 | `condition`: "oropharynx squamous cell carcinoma" = case, "normal" = control | Partially paired: 7 donors have both tumour and adjacent-normal samples, 3 donors tumour only. 21 `sample_id`/`library_id` (some donors have 2 tumour libraries). |
| cz_hpap_t1d_islets | `donor_id` | 5 / 19 / 24 | `condition`: "type 1 diabetes mellitus" = case, "normal" = control | Control arm pools 19 donors of which 9-ish are autoantibody-positive non-diabetic (`disease_state` = AAB, 20,650 cells) and the rest `Control` (26,594 cells). Cell-type cluster "unknown" (4,424 cells) was analysed as a cell type. |
| cz_human_kidney_v1.5 | `donor_id` (= `SampleID`) | 37 / 26 / 63 | `condition`: "chronic kidney disease" = case, "normal" = control | 67 `SpecimenID`/`LibraryID` vs 63 donors (a few donors with 2 specimens). Controls mix living donors (36,961 cells) and stone donors (6,419); CKD mixes DKD and HCKD. |
| cz_influenza | `donor_id` (= `Sample.ID`) | 5 / 4 / 9 | `condition`: "influenza" = case, "normal" = control | One sample per donor. |
| ERP | `patient_id` | 3 / 12 / 15 | `condition`: E = case, NE = control (T-cell clonotype expansion) | Unbalanced (3 expanders). Same BIOKEY cohort as TNBC, disjoint patients. |
| lupus | `donor_id` (= `ind_cov`, `sample_uuid`) | 26 / 22 / 48 | `condition`: "systemic lupus erythematosus" = case, "normal" = control | One sample per donor; 24 `library_uuid` (multiplexed pools spanning both arms). |
| MilCOVID_Azimuthl2 | `sample_id` (one sample per patient; = `orig.ident`, `gsm`, `sample_name`) | 3 / 5 / 8 | `condition`: Moderate = case, Healthy = control | The 5 healthy controls (COV07/08/09/17/18) are the same samples used as controls in SevCOVID_Azimuthl2 — the two comparisons are not independent. |
| sepsis | `donor_id` | 8 / 7 / 15 | `condition`: ICU-SEP = case, ICU-NoSEP = control | Clusters are cell states (e.g. Mono_MS1..MS4), not cell types. Cells come from two sorting gates (`biosample_id` CD45 / DC). |
| SevCOVID_Azimuthl2 | `sample_id` | 4 / 5 / 9 | `condition`: Severe = case, Healthy = control | Reproduces the verified anchor (9 donors: 4 severe, 5 healthy; k = 1). Controls shared with MilCOVID. |
| TNBC | `patient_id` | 5 / 7 / 12 | `condition`: E = case, NE = control | |

"Donors per condition" for a paired design (5yr_pic, cord_pic, covid, partially hnscc) counts
the same individuals in both arms; the roll-up's `n_donors_total` is the unique count.

## Things that looked inconsistent or need a footnote

1. **covid, cluster `CD14_plus_BDCA1_plus_PD_minus_L1_plus_cells` (C8)**: control arm has only
   20 cells and 0 meta-cells in `pseudobulk_seurat.rds`, yet the cluster was scored (600 case
   meta-cells). This is the documented paramPairings substitution (HANDOVER "C8 quirk"): its
   reference arm is the CD14+ monocyte control meta-cells. Correspondingly
   `CD14_plus_monocytes` control has 600 meta-cells where the rule would give 438 (= its case
   arm), because 600 control meta-cells were generated to serve as the C8 reference. These are
   the only rows with `metacell_rule_ok = FALSE`. Under the standard rule C8 would have been
   dropped (floor(20/3) = 6 < 100). The table reports the object as archived; a footnote
   should state the substituted reference for C8.
2. **covid meta-cell totals are unequal between arms** (5,135 case vs 4,697 control) purely
   because of item 1; every other comparison is exactly balanced.
3. **Cells "retained"**: `n_cells` counts all cells of the cluster/arm in the post-QC object.
   Because of the `1200*(k+1)` input cap and the k-neighbour aggregation, the number of cells
   actually absorbed into meta-cells is smaller for large clusters (e.g. lupus CD4 T:
   32,085 cells, 3,600 sampled, 600 meta-cells x 3 cells). `n_cells_input_to_metacells`
   gives the sampled count; exact absorbed-cell counts would need the meta-cell membership,
   which is not stored in the archive.
4. **Unbalanced or pooled control arms** worth flagging in the table legend:
   cz_cf_bronchial_biopsy (3 vs 19 donors, different labs/chemistries), cz_hpap_t1d_islets
   (AAB+ donors inside the control arm), cz_human_kidney_v1.5 (living + stone donors),
   ERP (3 vs 12).
5. **MilCOVID and SevCOVID share their 5 healthy control samples.**
6. **Cell types named "unknown" (islets) and "blood_cell" (influenza)** were analysed as
   clusters; they are CELLxGENE ontology labels carried through unchanged.
7. **BCG** has no donor/sample identifier at all in the archived object, so donor counts must
   come from the source publication (GSE232186) if required.
8. `n_candidate_lr_pairs` is identical on the case and control rows of a cell type (it is a
   per-cluster quantity); likewise `n_retained_tfs` and `n_predicted_interactions`. The
   roll-up sums them once per cell type (case rows only).

## Rerunning

```
cd analysis_r2
caffeinate -i Rscript dataset_summary/build_dataset_summary.R        # uses cache if present
rm -r dataset_summary/cache && caffeinate -i Rscript dataset_summary/build_dataset_summary.R   # full reload (~1.5 min)
```
Full reload of all 15 comparisons took ~80 s on the Mac mini (objects are read serially,
metadata extracted, object freed).

# In-silico TF knockout in CD16+ monocytes (SevCOVID_Azimuthl2) with CellOracle.
# Round-2 response to the reviewer's request for perturbational validation
# (round-1 R1.5 follow-up).
#
# Preprocessing through knn_imputation mirrors the manuscript's
# scripts/analysis_cellxgene_datasets/5_cell_oracle.py VERBATIM (same seeds,
# same parameters), so the fitted GRN corresponds to the manuscript's
# CD16_Mono GRN; the regenerated links table is cross-checked against the
# archived CSV before simulating.
#
# New (this script only): links.filter_links (manuscript params) ->
# fit_GRN_for_simulation -> simulate_shift KO per TF -> score each KO's
# predicted expression shift against the observed Severe-vs-Healthy axis.
#
# TF panels (from Decipher's own CD16_Mono scores, inputs/tf_ranking_cd16mono.csv):
#   top      — 10 highest mean |decipher_score| significant TFs
#   bottom   — 8 lowest mean |decipher_score| significant TFs (weak-signal control)
#   random   — 10 GRN TFs that are NOT significant regulons (negative control)
# Prediction: KO of a Severe-elevated TF (regulon_val > 0) shifts cells toward
# Healthy (negative cosine with the Healthy->Severe axis); KO of a
# Severe-suppressed TF (regulon_val < 0) shifts toward Severe (positive cosine).
# Controls should show weaker, sign-incoherent alignment.

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import random
import celloracle as co

print("celloracle version:", co.__version__, flush=True)

random.seed(123)
np.random.seed(123)

BASE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(BASE, "inputs")
OUT = os.path.join(BASE, "outputs")
os.makedirs(OUT, exist_ok=True)

N_PROPAGATION = 3
CONDITION_CANDIDATES = ["condition", "severity_group", "orig.condition", "condition_original"]

# ---------------- preprocessing: verbatim 5_cell_oracle.py ----------------
adata = sc.read_h5ad(os.path.join(INP, "CD16_Mono.h5ad"))
cluster_name = adata.obs["cluster"].unique()[0]
print("cluster:", cluster_name, "| cells:", adata.shape[0], flush=True)
print("obs columns:", list(adata.obs.columns), flush=True)

cond_col = next(c for c in CONDITION_CANDIDATES
                if c in adata.obs.columns and
                set(map(str, adata.obs[c].unique())) >= {"Severe", "Healthy"})
print("condition column:", cond_col, adata.obs[cond_col].value_counts().to_dict(), flush=True)
condition = adata.obs[cond_col].astype(str).copy()

base_GRN = co.data.load_human_promoter_base_GRN()

sc.pp.filter_genes(adata, min_counts=1)
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')
filter_result = sc.pp.filter_genes_dispersion(adata.X, n_top_genes=3000, log=False)
adata = adata[:, filter_result.gene_subset]
sc.pp.normalize_per_cell(adata)
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

oracle = co.Oracle()
adata.X = adata.layers["raw_count"].copy()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="cluster",
                                   embedding_name="X_pca")
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
k = int(0.025 * n_cell)
print(f"n_comps={n_comps} n_cell={n_cell} k={k}", flush=True)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8,
                      b_maxl=k * 4, n_jobs=6)

links = oracle.get_links(cluster_name_for_GRN_unit="cluster", alpha=10, verbose_level=10)

# ---------------- cross-check regenerated GRN vs archived CSV ----------------
regen = links.links_dict[cluster_name]
regen.to_csv(os.path.join(OUT, "regenerated_GRN_CD16_Mono.csv"))
arch = pd.read_csv(os.path.join(INP, "CD16_Mono.csv"), index_col=0)
m = regen.merge(arch, on=["source", "target"], suffixes=("_new", "_old"))
cor = np.corrcoef(m["coef_mean_new"], m["coef_mean_old"])[0, 1]
print(f"GRN check: regen {len(regen)} edges, archived {len(arch)}, "
      f"overlap {len(m)}, coef_mean r={cor:.4f}", flush=True)
pd.DataFrame([{"n_regen": len(regen), "n_arch": len(arch), "n_overlap": len(m),
               "coef_mean_r": cor}]).to_csv(os.path.join(OUT, "grn_check.csv"), index=False)

# ---------------- simulation setup (new) ----------------
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
oracle.to_hdf5(os.path.join(OUT, "cd16mono.celloracle.oracle"))

# ---------------- TF panels ----------------
rank = pd.read_csv(os.path.join(INP, "tf_ranking_cd16mono.csv"))
rank = rank[rank["in_grn"]].reset_index(drop=True)
sim_genes = set(oracle.adata.var_names)
active = set(oracle.active_regulatory_genes)

rank_ok = rank[rank["regulon"].isin(sim_genes & active)]
top_tfs = rank_ok.head(10)["regulon"].tolist()
bottom_tfs = rank_ok.tail(8)["regulon"].tolist()

grn_csv = pd.read_csv(os.path.join(INP, "CD16_Mono.csv"), index_col=0)
nonsig_pool = sorted((set(grn_csv["source"]) - set(rank["regulon"])) & sim_genes & active)
rng = np.random.RandomState(123)
random_tfs = list(rng.choice(nonsig_pool, size=min(10, len(nonsig_pool)), replace=False))

panels = ([(t, "top") for t in top_tfs] + [(t, "bottom") for t in bottom_tfs]
          + [(t, "random") for t in random_tfs])
print("panels:", panels, flush=True)

reg_val = dict(zip(rank["regulon"], rank["regulon_val"]))

# ---------------- severity axis ----------------
imputed = oracle.adata.layers["imputed_count"]
imputed = np.asarray(imputed.todense()) if hasattr(imputed, "todense") else np.asarray(imputed)
cond = condition.loc[oracle.adata.obs_names].values
sev_mask, hea_mask = cond == "Severe", cond == "Healthy"
axis = imputed[sev_mask].mean(axis=0) - imputed[hea_mask].mean(axis=0)   # Healthy -> Severe
axis_norm = axis / np.linalg.norm(axis)
print(f"axis built: {sev_mask.sum()} severe vs {hea_mask.sum()} healthy cells", flush=True)

def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan

# ---------------- KO sweep ----------------
rows = []
for tf, group in panels:
    try:
        oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=N_PROPAGATION)
        delta = oracle.adata.layers["delta_X"]
        delta = np.asarray(delta.todense()) if hasattr(delta, "todense") else np.asarray(delta)
        mean_delta = delta.mean(axis=0)
        per_cell = (delta @ axis_norm) / np.maximum(np.linalg.norm(delta, axis=1), 1e-12)
        rv = reg_val.get(tf, np.nan)
        rows.append({
            "tf": tf, "group": group, "regulon_val": rv,
            "cosine_mean_delta_vs_axis": cosine(mean_delta, axis),
            "mean_per_cell_cosine": float(np.nanmean(per_cell)),
            "frac_cells_negative_cosine": float(np.nanmean(per_cell < 0)),
            "shift_magnitude": float(np.linalg.norm(mean_delta)),
            "predicted_sign": -np.sign(rv) if not np.isnan(rv) else np.nan,
            "observed_sign": np.sign(cosine(mean_delta, axis)),
        })
        print(f"KO {tf} ({group}): cos={rows[-1]['cosine_mean_delta_vs_axis']:.3f} "
              f"mag={rows[-1]['shift_magnitude']:.3f}", flush=True)
    except Exception as e:
        print(f"KO {tf} FAILED: {e}", flush=True)
        rows.append({"tf": tf, "group": group, "regulon_val": reg_val.get(tf, np.nan),
                     "error": str(e)})

res = pd.DataFrame(rows)
res["sign_match"] = res.get("predicted_sign") == res.get("observed_sign")
res.to_csv(os.path.join(OUT, "ko_results_cd16mono.csv"), index=False)
print(res.to_string(), flush=True)
print("DONE", flush=True)

# In-silico TF KO — vaccination dataset (covid), CD16+ monocytes. Robustness
# companion to ko_cd16mono.py (severity); same design, axis = control(day0) ->
# case(day22). Preprocessing verbatim from 5_cell_oracle.py.
import os
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd, scanpy as sc, random
import celloracle as co
print("celloracle:", co.__version__, flush=True)
random.seed(123); np.random.seed(123)

BASE = os.path.dirname(os.path.abspath(__file__))
INP, OUT = os.path.join(BASE, "inputs"), os.path.join(BASE, "outputs")
H5AD, GRNCSV = "CD16_plus_monocytes.h5ad", "CD16_plus_monocytes.csv"
RANKCSV, TAG = "tf_ranking_cd16_vax.csv", "cd16_vax"

adata = sc.read_h5ad(os.path.join(INP, H5AD))
cluster_name = adata.obs["cluster"].unique()[0]
print("cluster:", cluster_name, "cells:", adata.shape[0], flush=True)
print("obs cols:", list(adata.obs.columns), flush=True)
cond_col = None
for c in ["condition", "orig.condition", "condition_original"]:
    if c in adata.obs.columns:
        vals = set(map(str, adata.obs[c].unique()))
        print(c, "->", vals, flush=True)
        if vals >= {"case", "control"}: cond_col, case_v, ctrl_v = c, "case", "control"; break
        if vals >= {"22", "0"}: cond_col, case_v, ctrl_v = c, "22", "0"; break
assert cond_col, "no condition column resolved"
print("using:", cond_col, case_v, ctrl_v, flush=True)
condition = adata.obs[cond_col].astype(str).copy()

base_GRN = co.data.load_human_promoter_base_GRN()
sc.pp.filter_genes(adata, min_counts=1)
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')
fr = sc.pp.filter_genes_dispersion(adata.X, n_top_genes=3000, log=False)
adata = adata[:, fr.gene_subset]
sc.pp.normalize_per_cell(adata)
adata.raw = adata; adata.layers["raw_count"] = adata.raw.X.copy()
sc.pp.log1p(adata); sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
oracle = co.Oracle()
adata.X = adata.layers["raw_count"].copy()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="cluster", embedding_name="X_pca")
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0]
n_comps = min(n_comps, 50)
k = int(0.025 * oracle.adata.shape[0])
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=6)
links = oracle.get_links(cluster_name_for_GRN_unit="cluster", alpha=10, verbose_level=10)

regen = links.links_dict[cluster_name]
arch = pd.read_csv(os.path.join(INP, GRNCSV), index_col=0)
m = regen.merge(arch, on=["source","target"], suffixes=("_new","_old"))
cor = np.corrcoef(m["coef_mean_new"], m["coef_mean_old"])[0,1]
print(f"GRN check: {len(regen)} vs {len(arch)}, overlap {len(m)}, r={cor:.4f}", flush=True)

links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

rank = pd.read_csv(os.path.join(INP, RANKCSV))
rank = rank[rank["in_grn"]].reset_index(drop=True)
sim_genes, active = set(oracle.adata.var_names), set(oracle.active_regulatory_genes)
rank_ok = rank[rank["regulon"].isin(sim_genes & active)]
top_tfs = rank_ok.head(10)["regulon"].tolist()
bottom_tfs = rank_ok.tail(8)["regulon"].tolist()
nonsig = sorted((set(arch["source"]) - set(rank["regulon"])) & sim_genes & active)
rng = np.random.RandomState(123)
random_tfs = list(rng.choice(nonsig, size=min(10, len(nonsig)), replace=False))
panels = [(t,"top") for t in top_tfs] + [(t,"bottom") for t in bottom_tfs] + [(t,"random") for t in random_tfs]
print("panels:", panels, flush=True)
reg_val = dict(zip(rank["regulon"], rank["regulon_val"]))

imp = oracle.adata.layers["imputed_count"]
imp = np.asarray(imp.todense()) if hasattr(imp, "todense") else np.asarray(imp)
cond = condition.loc[oracle.adata.obs_names].values
case_m, ctrl_m = cond == case_v, cond == ctrl_v
axis = imp[case_m].mean(axis=0) - imp[ctrl_m].mean(axis=0)
axis_n = axis / np.linalg.norm(axis)
print(f"axis: {case_m.sum()} case vs {ctrl_m.sum()} control", flush=True)

rows = []
for tf, group in panels:
    try:
        oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
        d = oracle.adata.layers["delta_X"]
        d = np.asarray(d.todense()) if hasattr(d, "todense") else np.asarray(d)
        md = d.mean(axis=0)
        cosv = float(md @ axis_n / max(np.linalg.norm(md), 1e-12) / np.linalg.norm(axis_n) * np.linalg.norm(axis_n)) if np.linalg.norm(md) > 0 else np.nan
        cosv = float(md @ axis_n / np.linalg.norm(md)) if np.linalg.norm(md) > 0 else np.nan
        rv = reg_val.get(tf, np.nan)
        rows.append({"tf": tf, "group": group, "regulon_val": rv,
                     "cosine_mean_delta_vs_axis": cosv,
                     "shift_magnitude": float(np.linalg.norm(md)),
                     "observed_sign": float(np.sign(cosv))})
        print(f"KO {tf} ({group}): cos={cosv:.3f} mag={rows[-1]['shift_magnitude']:.3f}", flush=True)
    except Exception as e:
        print(f"KO {tf} FAILED: {e}", flush=True)
        rows.append({"tf": tf, "group": group, "error": str(e)})
pd.DataFrame(rows).to_csv(os.path.join(OUT, f"ko_results_{TAG}.csv"), index=False)
print("DONE", flush=True)

# Hybrid LR perturbation, step 2 (23 Aug 2026): translate each LR knockout's predicted TF
# activity changes (lr_ko_tf_deltas.csv, from Decipher's own RF models) into a CellOracle
# multi-TF perturbation (expression direction from the orientation anchors; down -> 0,
# up -> 95th percentile), simulate, and score the shift against the severity axis.
# NOT orthogonal validation: first hop is Decipher predicting Decipher. Consistency check.
import os, numpy as np, pandas as pd, celloracle as co
BASE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(BASE, "outputs")
oracle = co.load_hdf5(os.path.join(OUT, "cd16mono.celloracle.oracle"))
imp = oracle.adata.layers["imputed_count"]; imp = np.asarray(imp.todense()) if hasattr(imp, "todense") else np.asarray(imp)
cond = oracle.adata.obs["disease_status"].astype(str).values
axis = imp[cond == "COVID-19"].mean(axis=0) - imp[cond == "Healthy"].mean(axis=0)
genes = list(oracle.adata.var_names); gidx = {g: i for i, g in enumerate(genes)}
q95 = {g: float(np.quantile(imp[:, gidx[g]], 0.95)) for g in genes}
ok = set(genes) & set(oracle.active_regulatory_genes)
def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b); return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan
def simulate(pert):
    oracle.simulate_shift(perturb_condition=pert, n_propagation=3)
    d = oracle.adata.layers["delta_X"]; d = np.asarray(d.todense()) if hasattr(d, "todense") else np.asarray(d)
    md = d.mean(axis=0); return cosine(md, axis), float(np.linalg.norm(md))
dz = pd.read_csv(os.path.join(OUT, "lr_ko_tf_deltas.csv"))
dz = dz[dz.tf.isin(ok) & dz.expr_dir.notna()]
KINDS = os.environ.get("KINDS"); SUFFIX = os.environ.get("OUT_SUFFIX", ""); SKIP_NULL = os.environ.get("SKIP_NULL") == "1"
if KINDS: dz = dz[dz.kind.isin(KINDS.split(","))]
rows = []
# 23 Aug fix: a TF enters the perturbation only if its predicted activity change is at least
# THR standard deviations; a knockout with no TF passing produces no shift (the right answer
# for a weak LR pair). Two thresholds; at most 10 TFs.
for ko, g in dz.groupby("ko"):
    for thr in (0.1, 0.25):
        top = g[g.dz.abs() >= thr].reindex(g[g.dz.abs() >= thr].dz.abs().sort_values(ascending=False).index).head(10)
        base = {"ko": ko, "panel": g.panel.iloc[0], "kind": g.kind.iloc[0], "anchor_pair": g.anchor_pair.iloc[0], "thr": thr,
                "n_tfs": int(len(top)), "tfs": ";".join(top.tf), "dirs": ";".join(str(int(e)) for e in top.expr_dir),
                "max_abs_dz": float(g.dz.abs().max())}
        if len(top) == 0:
            rows.append({**base, "cosine": np.nan, "magnitude": 0.0}); print(f"{ko} thr={thr}: no TF passes -> no shift", flush=True); continue
        pert = {t: (0.0 if e < 0 else q95[t]) for t, e in zip(top.tf, top.expr_dir)}
        try:
            c, m = simulate(pert)
        except Exception as ex:
            print(f"{ko} thr={thr} FAILED: {ex}", flush=True); continue
        rows.append({**base, "cosine": c, "magnitude": m})
        print(f"{ko} [{g.panel.iloc[0]}/{g.kind.iloc[0]}] thr={thr} n={len(top)}: cos={c:.3f} mag={m:.3f} tfs={','.join(top.tf)}", flush=True)
# null: random TF sets, sizes matched to what top-panel knockouts perturb, random directions
rng = np.random.RandomState(123); pool = sorted(ok)
for k in ((3, 5, 7) if not SKIP_NULL else ()):
    for r in range(12):
        tfs = list(rng.choice(pool, size=k, replace=False)); dirs = rng.choice([-1, 1], size=k)
        pert = {t: (0.0 if e < 0 else q95[t]) for t, e in zip(tfs, dirs)}
        try: c, m = simulate(pert)
        except Exception as ex: print(f"null k={k} r={r} FAILED: {ex}", flush=True); continue
        rows.append({"ko": f"null_{k}_{r}", "panel": "null", "kind": "random_tf_set", "anchor_pair": "", "thr": np.nan, "n_tfs": k,
                     "tfs": ";".join(tfs), "dirs": ";".join(str(int(e)) for e in dirs), "max_abs_dz": np.nan, "cosine": c, "magnitude": m})
        print(f"null k={k} r={r}: cos={c:.3f} mag={m:.3f}", flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, f"lr_ko_sim_results{SUFFIX}.csv"), index=False); print("LRKO DONE", flush=True)

# Extra KO panels (23 Aug 2026): do importance-only vs dTF-only rankings pick TFs whose
# knockout aligns with the severity axis? Reuses the saved oracle (same GRN/fit as the
# main run). ZBTB7A re-run as a reproduction check of the reloaded object.
import os, numpy as np, pandas as pd, celloracle as co
BASE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(BASE, "outputs")
oracle = co.load_hdf5(os.path.join(OUT, "cd16mono.celloracle.oracle"))
imp = oracle.adata.layers["imputed_count"]; imp = np.asarray(imp.todense()) if hasattr(imp, "todense") else np.asarray(imp)
cond = oracle.adata.obs["disease_status"].astype(str).values
axis = imp[cond == "COVID-19"].mean(axis=0) - imp[cond == "Healthy"].mean(axis=0)
def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b); return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan
rank = pd.read_csv(os.path.join(BASE, "inputs", "tf_ranking_cd16mono.csv")); reg_val = dict(zip(rank.regulon, rank.regulon_val))
panels = [("ZBTB7A", "check")] + [(t, "imp_only") for t in ["IRF1","STAT2","NFATC2","JUND","ETS2","JUN","EGR2"]] \
       + [(t, "dtf_only") for t in ["HES4","NR4A2","E2F5","FOXL2","EGR3"]]
ok = set(oracle.adata.var_names) & set(oracle.active_regulatory_genes)
rows = []
for tf, g in panels:
    if tf not in ok: print(f"SKIP {tf}: not simulable", flush=True); rows.append({"tf": tf, "group": g, "error": "not simulable"}); continue
    oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
    d = oracle.adata.layers["delta_X"]; d = np.asarray(d.todense()) if hasattr(d, "todense") else np.asarray(d)
    md = d.mean(axis=0)
    rows.append({"tf": tf, "group": g, "regulon_val": reg_val.get(tf, np.nan), "cosine_mean_delta_vs_axis": cosine(md, axis),
                 "shift_magnitude": float(np.linalg.norm(md))})
    print(f"KO {tf} ({g}): cos={rows[-1]['cosine_mean_delta_vs_axis']:.3f} mag={rows[-1]['shift_magnitude']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "ko_results_cd16mono_extra.csv"), index=False); print("EXTRA DONE", flush=True)

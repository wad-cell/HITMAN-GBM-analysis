"""
B3 - step02: PRIMARY DESeq2 ~MAN*Field on all 20 samples (B3.4 + base for B3.5)
Full-gene results for C1-C5; saves per-contrast TSVs + model artifacts.
"""
import os, time, json, pickle, traceback
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")
LOG  = os.path.join(TEMP, "b3_step02_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)

try:
    t_start = time.time()
    meta = pd.read_pickle(os.path.join(TEMP, "b3_meta.pkl"))
    genes = pd.read_csv(os.path.join(TEMP, "b3_genes.tsv"), sep="\t")["gene_id"].tolist()
    rounded = np.load(os.path.join(TEMP, "b3_counts_rounded.npy"))
    counts_df = pd.DataFrame(rounded.T, index=meta["sample"].values, columns=genes, dtype=float)
    counts_df = counts_df.loc[meta.index]
    log("prepared", counts_df.shape)

    dds = DeseqDataSet(counts=counts_df, metadata=meta, design="~MAN*Field",
                       refit_cooks=True, quiet=False)
    log("dds built", round(time.time()-t_start,1), "design cols:", list(dds.obsm["design_matrix"].columns))
    dds.deseq2()
    log("deseq2 done", round(time.time()-t_start,1))

    cols = list(dds.obsm["design_matrix"].columns)  # Intercept, MAN[T.MAN], Field[T.F], MAN[T.MAN]:Field[T.F]
    idx = {c: i for i, c in enumerate(cols)}
    contrasts = {
        "C1": "ND-F - ND-NF (Field without MAN)",
        "C2": "MAN-NF - ND-NF (MAN without Field)",
        "C3": "MAN-F - MAN-NF (Field with MAN)",
        "C4": "MAN-F - ND-NF (Total HITMAN-associated)",
        "C5": "(MAN-F-MAN-NF)-(ND-F-ND-NF) = MANxField interaction [PRIMARY]",
    }
    vecs = {
        "C1": np.array([0,0,1,0], float),
        "C2": np.array([0,1,0,0], float),
        "C3": np.array([0,0,1,1], float),
        "C4": np.array([0,1,1,1], float),
        "C5": np.array([0,0,0,1], float),
    }

    summaries = {}
    for cname in ["C1","C2","C3","C4","C5"]:
        t0 = time.time()
        st = DeseqStats(dds, contrast=vecs[cname], quiet=True)
        st.summary()
        res = st.results_df.copy()
        res = res.reset_index().rename(columns={"index": "gene"})
        res.columns = ["gene","baseMean","log2FC","lfcSE","stat","pvalue","padj"]
        res["contrast"] = cname
        res.to_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{cname}.tsv"), sep="\t", index=False)
        sig = res.dropna(subset=["padj"])
        n_sig = int((sig["padj"] < 0.05).sum())
        sig05 = sig[sig["padj"] < 0.05]
        n_up = int((sig05["log2FC"] > 0).sum()); n_down = int((sig05["log2FC"] < 0).sum())
        summaries[cname] = {
            "tested_genes_with_padj": int(len(sig)),
            "n_FDR005": n_sig, "up": n_up, "down": n_down,
            "median_log2FC_all": float(res["log2FC"].median()),
            "mad_log2FC_all": float((res["log2FC"] - res["log2FC"].median()).abs().median()),
            "q1_log2FC_sig": float(sig05["log2FC"].quantile(0.25)) if n_sig else None,
            "q3_log2FC_sig": float(sig05["log2FC"].quantile(0.75)) if n_sig else None,
            "seconds": round(time.time()-t0,1),
        }
        log("contrast", cname, "sig", n_sig, "up", n_up, "down", n_down, "t", round(time.time()-t0,1))

    # save artifacts
    with open(os.path.join(TEMP, "b3_primary_artifacts.pkl"), "wb") as f:
        pickle.dump({
            "dds": dds, "design_cols": cols, "idx": idx,
            "contrasts": contrasts, "vecs": vecs,
            "summaries": summaries,
            "elapsed": round(time.time()-t_start,1),
            "sample_names": list(meta["sample"].values),
            "design_matrix": dds.obsm["design_matrix"],
            "cooks": dds.obsm.get("cooks_outliers", None),
        }, f)
    with open(os.path.join(TEMP, "b3_primary_summaries.json"), "w") as f:
        json.dump(summaries, f, indent=2)
    log("ARTIFACTS SAVED elapsed", round(time.time()-t_start,1))
    log("STEP02 OK")
except Exception as e:
    log("STEP02 EXC", traceback.format_exc())

"""
B3 - step03: SENSITIVITY DESeq2 without T.Cells_Q (19 samples) + primary comparison (B3.6)
"""
import os, time, json, pickle, traceback
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")
LOG  = os.path.join(TEMP, "b3_step03_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)

try:
    t_start = time.time()
    meta = pd.read_pickle(os.path.join(TEMP, "b3_meta.pkl"))
    genes = pd.read_csv(os.path.join(TEMP, "b3_genes.tsv"), sep="\t")["gene_id"].tolist()
    rounded = np.load(os.path.join(TEMP, "b3_counts_rounded.npy"))
    # drop T.Cells_Q
    keep = [s for s in meta["sample"].values if s != "T.Cells_Q"]
    log("sensitivity samples:", len(keep), "dropped T.Cells_Q")
    meta_s = meta.loc[keep].copy()
    rounded_s = rounded[:, [meta.index.get_loc(s) for s in keep]]

    counts_df = pd.DataFrame(rounded_s.T, index=meta_s["sample"].values, columns=genes, dtype=float)
    counts_df = counts_df.loc[meta_s.index]
    dds_s = DeseqDataSet(counts=counts_df, metadata=meta_s, design="~MAN*Field",
                         refit_cooks=True, quiet=False)
    dds_s.deseq2()
    log("sens deseq2 done", round(time.time()-t_start,1))
    cols_s = list(dds_s.obsm["design_matrix"].columns)
    log("sens design cols", cols_s)

    vecs = {
        "C1": np.array([0,0,1,0], float),
        "C2": np.array([0,1,0,0], float),
        "C3": np.array([0,0,1,1], float),
        "C4": np.array([0,1,1,1], float),
        "C5": np.array([0,0,0,1], float),
    }
    sens_res = {}
    for cname in ["C1","C2","C3","C4","C5"]:
        st = DeseqStats(dds_s, contrast=vecs[cname], quiet=True)
        st.summary()
        res = st.results_df.copy().reset_index().rename(columns={"index": "gene"})
        res.columns = ["gene","baseMean","log2FC","lfcSE","stat","pvalue","padj"]
        sens_res[cname] = res
    with open(os.path.join(TEMP, "b3_sensitivity_results.pkl"), "wb") as f:
        pickle.dump({"results": sens_res, "meta_s": meta_s,
                     "dropped": "T.Cells_Q", "n": len(keep),
                     "design_matrix": dds_s.obsm["design_matrix"],
                     "cooks": dds_s.obsm.get("cooks_outliers", None)}, f)
    log("sens artifacts saved", round(time.time()-t_start,1))
    log("STEP03 OK")
except Exception as e:
    log("STEP03 EXC", traceback.format_exc())

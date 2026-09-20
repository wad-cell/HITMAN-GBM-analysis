"""B3 - step07: ranked GSEA C1-C5 x Hallmark/Reactome/GOBP (B3.7)"""
import os, json, time, traceback
import numpy as np, pandas as pd
import gseapy as gp
BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")
LOG  = os.path.join(TEMP, "b3_step07_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)
try:
    mapf = pd.read_csv(os.path.join(TEMP, "b3_ensg_symbol_map.tsv"), sep="\t", dtype=str)
    sym = dict(zip(mapf["gene_id"], mapf["symbol"]))
    log("symbol map", len(sym))

    rnks = {}
    for c in ["C1","C2","C3","C4","C5"]:
        df = pd.read_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{c}.tsv"), sep="\t").set_index("gene")
        df["core"] = [g.split(".")[0] for g in df.index]
        df["symbol"] = df["core"].map(sym)
        d2 = df[df["symbol"].notna() & df["stat"].notna() & np.isfinite(df["stat"])].copy()
        # deduplicate symbol by max |stat|
        d2["astat"] = d2["stat"].abs()
        keep = d2.sort_values("astat", ascending=False).groupby("symbol").head(1)
        rnk = pd.DataFrame({"gene": keep["symbol"].values, "score": keep["stat"].values})
        rnk = rnk.dropna().drop_duplicates(subset="gene")
        rnk.to_csv(os.path.join(TEMP, f"b3_rnk_{c}.rnk"), sep="\t", index=False, header=False)
        log(c, "rnk genes", len(rnk), "pos", int((rnk['score']>0).sum()), "neg", int((rnk['score']<0).sum()))

    libs = {"HALLMARK": "b3_gslib_HALLMARK.gmt",
            "REACTOME": "b3_gslib_REACTOME.gmt",
            "GOBP": "b3_gslib_GOBP.gmt"}
    all_rows = []
    for c in ["C1","C2","C3","C4","C5"]:
        rnk_path = os.path.join(TEMP, f"b3_rnk_{c}.rnk")
        for lib, gmt in libs.items():
            t0 = time.time()
            res = gp.prerank(rnk=rnk_path, gene_sets=os.path.join(TEMP, gmt),
                             min_size=15, max_size=500, permutation_num=1000,
                             threads=8, seed=123, no_plot=True, processes=1,
                             ascending=False, outdir=None)
            rr = res.res2d
            rr.insert(0, "contrast", c); rr.insert(1, "library", lib)
            all_rows.append(rr)
            log(lib, c, "rows", len(rr), "sec", round(time.time()-t0,1))
    big = pd.concat(all_rows, ignore_index=True)
    big.to_csv(os.path.join(OUT, "B3_GSEA_ALL.tsv"), sep="\t", index=False)
    for lib in libs:
        sub = big[big["library"] == lib]
        sub.to_csv(os.path.join(OUT, f"B3_GSEA_{lib}.tsv"), sep="\t", index=False)
    log("STEP07 OK total rows", len(big))
except Exception:
    log("STEP07 EXC", traceback.format_exc())

"""B3.5 supplementary deep-dive: extreme stat genes, C5 set differences, Py filtering audit."""
import os, json
import numpy as np
import pandas as pd

BASE = r"${PROJECT_ROOT}"
OUT = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")

def load_py(cn):
    df = pd.read_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{cn}.tsv"), sep="\t")
    for c in ["baseMean","log2FC","lfcSE","stat","pvalue","padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("gene")

def load_r(cn):
    df = pd.read_csv(os.path.join(TEMP, f"B35_R_{cn}.tsv"), sep="\t")
    for c in ["baseMean","log2FC","lfcSE","stat","pvalue","padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("gene")

# ---- 1. extreme stat values C2 / C5 ----
for cn in ["C2","C5"]:
    py = load_py(cn); r = load_r(cn)
    j = pd.concat([py["stat"], r["stat"]], axis=1, keys=["py","r"]).dropna()
    j["abs_diff"] = (j["py"]-j["r"]).abs()
    print("=== ", cn, " stat extremes ===")
    print("py stat max", j["py"].max(), "min", j["py"].min())
    print("r  stat max", j["r"].max(), "min", j["r"].min())
    print(j.sort_values("abs_diff", ascending=False).head(10))

# ---- 2. log2FC largest diffs for C4/C5 ----
for cn in ["C4","C5"]:
    py = load_py(cn); r = load_r(cn)
    j = pd.concat([py["log2FC"], r["log2FC"]], axis=1, keys=["py","r"]).dropna()
    j["abs_diff"] = (j["py"]-j["r"]).abs()
    print("=== ", cn, " log2FC largest diffs ===")
    print(j.sort_values("abs_diff", ascending=False).head(10))

# ---- 3. Py-side filtering audit ----
rows = []
for cn in ["C1","C2","C3","C4","C5"]:
    py = load_py(cn)
    input_universe = len(py)
    stat_avail = py["stat"].notna().sum()
    pval_avail = py["pvalue"].notna().sum()
    padj_eval = py["padj"].notna().sum()
    padj_na_with_pval = ((py["padj"].isna()) & (py["pvalue"].notna())).sum()
    removed = pval_avail - padj_eval
    fdr = (py["padj"] < 0.05).sum()
    rows.append(dict(contrast=cn, input_universe=input_universe, statistic_available=int(stat_avail),
                     pvalue_available=int(pval_avail), padj_evaluable=int(padj_eval),
                     padj_NA=int(padj_na_with_pval), removed_by_independent_filtering=int(removed),
                     FDR_005=int(fdr)))
    print(rows[-1])

# R-side audit reload
r_rows = []
for cn in ["C1","C2","C3","C4","C5"]:
    r = load_r(cn)
    input_universe = len(r)
    stat_avail = r["stat"].notna().sum()
    pval_avail = r["pvalue"].notna().sum()
    padj_eval = r["padj"].notna().sum()
    padj_na_with_pval = ((r["padj"].isna()) & (r["pvalue"].notna())).sum()
    removed = pval_avail - padj_eval
    fdr = (r["padj"] < 0.05).sum()
    r_rows.append(dict(contrast=cn, input_universe=input_universe, statistic_available=int(stat_avail),
                       pvalue_available=int(pval_avail), padj_evaluable=int(padj_eval),
                       padj_NA=int(padj_na_with_pval), removed_by_independent_filtering=int(removed),
                       FDR_005=int(fdr)))
    print("R", r_rows[-1])

pd.DataFrame(rows).to_csv(os.path.join(TEMP, "B35_filtering_audit_py.tsv"), sep="\t", index=False)
pd.DataFrame(r_rows).to_csv(os.path.join(TEMP, "B35_filtering_audit_r_full.tsv"), sep="\t", index=False)

# ---- 4. C5 discordant genes (sig in one engine not other) ----
py = load_py("C5"); r = load_r("C5")
py_sig = set(py.index[py["padj"] < 0.05]); r_sig = set(r.index[r["padj"] < 0.05])
only_py = py_sig - r_sig; only_r = r_sig - py_sig
print("C5 sig only in py:", len(only_py), "only in r:", len(only_r))
# characterize only_r genes by py padj proximity
r_disc = r.loc[list(only_r), ["baseMean","log2FC","stat","padj"]].copy()
py_padj_for = py["padj"].reindex(list(only_r))
r_disc["py_padj"] = py_padj_for.values
print("only_r median baseMean", r_disc["baseMean"].median(), "median py_padj", r_disc["py_padj"].median())
print(r_disc.sort_values("py_padj").head(8))

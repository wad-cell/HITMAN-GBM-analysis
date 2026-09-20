"""B3.5: PyDESeq2 vs R/DESeq2 engine comparison + filtering audit + gene annotation."""
import os, json, math
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = r"${PROJECT_ROOT}"
OUT = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")

# ---------------- load helpers ----------------
def load_py(contrast):
    df = pd.read_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{contrast}.tsv"), sep="\t")
    df = df.rename(columns={"gene": "gene"}).copy()
    for c in ["baseMean","log2FC","lfcSE","stat","pvalue","padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("gene")

def load_r(contrast):
    df = pd.read_csv(os.path.join(TEMP, f"B35_R_{contrast}.tsv"), sep="\t")
    df = df.rename(columns={"gene": "gene"}).copy()
    for c in ["baseMean","log2FC","lfcSE","stat","pvalue","padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("gene")

def pearson_spearman(a, b):
    m = pd.concat([a, b], axis=1).dropna()
    if len(m) < 3:
        return None, None, 0
    pr = sps.pearsonr(m.iloc[:,0], m.iloc[:,1])[0]
    sr = sps.spearmanr(m.iloc[:,0], m.iloc[:,1])[0]
    return float(pr), float(sr), len(m)

def jaccard(s1, s2):
    u = len(s1 | s2)
    return float(len(s1 & s2) / u) if u else None

def concordance_metrics(py, r, contrast, stable=None):
    """Full comparison for one contrast."""
    both = pd.concat([py["log2FC"], r["log2FC"]], axis=1, keys=["py","r"]).dropna()
    out = {}
    # log2FC correlations on genes with both values
    pr, sr, n_lfc = pearson_spearman(py["log2FC"], r["log2FC"])
    out["n_both_log2FC"] = n_lfc
    out["log2FC_pearson"] = pr
    out["log2FC_spearman"] = sr
    # statistic correlation
    prs, srs, n_stat = pearson_spearman(py["stat"], r["stat"])
    out["n_both_stat"] = n_stat
    out["stat_pearson"] = prs
    out["stat_spearman"] = srs
    # sign concordance on genes with both stat non-NA (non-zero FC)
    fc = pd.concat([py["log2FC"], r["log2FC"]], axis=1, keys=["py","r"]).dropna()
    fc = fc[(fc["py"] != 0) & (fc["r"] != 0)]
    if len(fc):
        out["sign_concordance"] = float((np.sign(fc["py"]) == np.sign(fc["r"])).mean())
    else:
        out["sign_concordance"] = None
    out["n_sign"] = len(fc)
    # FDR counts
    py_sig = set(py.index[py["padj"] < 0.05])
    r_sig = set(r.index[r["padj"] < 0.05])
    out["py_n_FDR005"] = len(py_sig)
    out["r_n_FDR005"] = len(r_sig)
    out["sig_jaccard"] = jaccard(py_sig, r_sig)
    out["py_sig_in_r"] = float(len(py_sig & r_sig) / len(py_sig)) if py_sig else None
    out["r_sig_in_py"] = float(len(py_sig & r_sig) / len(r_sig)) if r_sig else None
    # top N overlap by |log2FC| (signed ranking by stat)
    def top_overlap(pyfc, rfc, N):
        a = pyfc.dropna().sort_values("py", key=lambda s: -s.abs()) if False else None
        return None
    # use stat rank (like DESeq2 default ordering), report overlap of union-top
    def topN_overlap(py_series, r_series, N, key="stat"):
        a = pd.concat([py_series, r_series], axis=1, keys=["py","r"]).dropna()
        py_top = set(a["py"].sort_values(ascending=False).head(N).index)
        r_top = set(a["r"].sort_values(ascending=False).head(N).index)
        return float(len(py_top & r_top) / N)
    out["top50_overlap"] = topN_overlap(py["stat"], r["stat"], 50)
    out["top100_overlap"] = topN_overlap(py["stat"], r["stat"], 100)
    out["top250_overlap"] = topN_overlap(py["stat"], r["stat"], 250)
    # C5-specific: original 1891 same-direction replication
    if contrast == "C5":
        py_sig_genes = py_sig
        r_fc = r["log2FC"]
        repl = pd.DataFrame({"py_fc": py.loc[list(py_sig_genes), "log2FC"],
                             "r_fc": r_fc.reindex(list(py_sig_genes))}).dropna()
        same = float((np.sign(repl["py_fc"]) == np.sign(repl["r_fc"])).mean())
        out["py1891_replicated_same_dir"] = same
        out["py1891_replicated_n"] = len(repl)
        # stable 1602
        if stable is not None:
            st = set(stable) & py_sig_genes
            repl2 = pd.DataFrame({"py_fc": py.loc[list(st), "log2FC"],
                                  "r_fc": r_fc.reindex(list(st))}).dropna()
            same2 = float((np.sign(repl2["py_fc"]) == np.sign(repl2["r_fc"])).mean())
            out["stable1602_same_dir"] = same2
            out["stable1602_n"] = len(repl2)
            out["stable1602_in_py_sig"] = len(st)
            # stable in R sig too
            st_r = set(st) & r_sig
            out["stable1602_in_r_sig"] = len(st_r)
    return out

# ---------------- main ----------------
contrasts = ["C1","C2","C3","C4","C5"]
stable_df = pd.read_csv(os.path.join(TEMP, "b3_stable_genes_c5.tsv"), sep="\t")
stable_genes = set(stable_df["gene"].dropna())

results = {}
for cn in contrasts:
    py = load_py(cn); r = load_r(cn)
    results[cn] = concordance_metrics(py, r, cn, stable=stable_genes)
    print(cn, json.dumps(results[cn], default=str)[:800])

# per-contrast summary table
comp_rows = []
for cn in contrasts:
    row = {"contrast": cn}
    row.update(results[cn])
    comp_rows.append(row)
comp = pd.DataFrame(comp_rows)
comp.to_csv(os.path.join(TEMP, "B35_engine_comparison_all.tsv"), sep="\t", index=False)
print("\nFull comparison saved.")

# C4/C5 dedicated tables with gene-level concordance export
for cn in ["C4","C5"]:
    py = load_py(cn); r = load_r(cn)
    gene_df = pd.DataFrame({
        "gene": py.index,
        "py_log2FC": py["log2FC"].values,
        "py_stat": py["stat"].values,
        "py_pvalue": py["pvalue"].values,
        "py_padj": py["padj"].values,
        "r_log2FC": r["log2FC"].reindex(py.index).values,
        "r_stat": r["stat"].reindex(py.index).values,
        "r_pvalue": r["pvalue"].reindex(py.index).values,
        "r_padj": r["padj"].reindex(py.index).values,
    })
    gene_df["log2FC_diff"] = gene_df["py_log2FC"] - gene_df["r_log2FC"]
    gene_df["stat_diff"] = gene_df["py_stat"] - gene_df["r_stat"]
    gene_df["py_sig005"] = gene_df["py_padj"] < 0.05
    gene_df["r_sig005"] = gene_df["r_padj"] < 0.05
    gene_df["sign_same"] = np.sign(gene_df["py_log2FC"]) == np.sign(gene_df["r_log2FC"])
    gene_df.to_csv(os.path.join(TEMP, f"B35_{cn}_engine_comparison_full.tsv"), sep="\t", index=False)
    print(cn, "gene-level rows", len(gene_df))
print("DONE")

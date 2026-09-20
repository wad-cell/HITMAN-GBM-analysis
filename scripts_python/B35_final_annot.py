"""B3.5 final: robust correlations (extreme-gene trimmed), gene annotation, deliverable TSVs."""
import os, json
import numpy as np
import pandas as pd
from scipy import stats as sps

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

# ---- robust Pearson after excluding |stat|>100 in either engine (Py divergent low-count genes) ----
robust = {}
for cn in ["C1","C2","C3","C4","C5"]:
    py = load_py(cn); r = load_r(cn)
    j = pd.concat([py["stat"], r["stat"]], axis=1, keys=["py","r"]).dropna()
    keep = (j["py"].abs() <= 100) & (j["r"].abs() <= 100)
    jk = j[keep]
    pr = float(sps.pearsonr(jk["py"], jk["r"])[0])
    # log2FC trimmed similarly by |stat|<=100 (already the divergent ones)
    fc = pd.concat([py["log2FC"], r["log2FC"]], axis=1, keys=["py","r"]).dropna()
    fc = fc[(j["py"].abs() <= 100) & (j["r"].abs() <= 100)] if len(fc)==len(j) else fc
    pr_fc = float(sps.pearsonr(fc["py"], fc["r"])[0])
    robust[cn] = dict(n_kept=len(jk), n_removed_extreme=len(j)-len(jk),
                      stat_pearson_trimmed=pr, log2FC_pearson_trimmed=pr_fc)
    print(cn, robust[cn])

pd.DataFrame(robust).T.to_csv(os.path.join(TEMP, "B35_robust_correlations.tsv"), sep="\t")

# ---- gene annotation (from Ensembl BioMart file) ----
ann = pd.read_csv(os.path.join(TEMP, "B35_ensembl_annotation.tsv"), sep="\t", dtype=str)
print("annotation rows", len(ann), ann.columns.tolist())
# ann columns: gene_id(unversioned ENSG), gene_id_version(ENSG.ver), symbol, biotype
ann_map = ann.set_index("gene_id_version")

genes = pd.read_csv(os.path.join(TEMP, "b3_genes.tsv"), sep="\t")
print("gene list rows", len(genes), genes.columns.tolist())

g = genes["gene_id"].astype(str)
unver = g.str.replace(r"\.[0-9]+$", "", regex=True)
out_df = pd.DataFrame({"ensg_versioned": g, "ensg_unversioned": unver})
sub = ann.set_index("gene_id")
# versioned match
v_match = ann_map.reindex(g)
out_df["hgnc_symbol"] = v_match["symbol"].values
out_df["gene_biotype"] = v_match["biotype"].values
# fill missing via unversioned
miss = out_df["hgnc_symbol"].isna()
if miss.any():
    u2 = sub.reindex(unver[miss])
    out_df.loc[miss, "hgnc_symbol"] = u2["symbol"].values
    out_df.loc[miss, "gene_biotype"] = u2["biotype"].values
out_df["hgnc_symbol"] = out_df["hgnc_symbol"].fillna("")
out_df["gene_biotype"] = out_df["gene_biotype"].fillna("")
print("symbol coverage:", (out_df["hgnc_symbol"]!="").mean().round(4),
      "biotype coverage:", (out_df["gene_biotype"]!="").mean().round(4))
out_df.to_csv(os.path.join(TEMP, "B35_gene_annotation_full.tsv"), sep="\t", index=False)
print(out_df.head(3))
print("annotation OK")

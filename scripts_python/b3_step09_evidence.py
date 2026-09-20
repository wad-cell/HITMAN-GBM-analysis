# -*- coding: utf-8 -*-
"""B3.10 evidence: C5 stable genes + C5 sensitivity GSEA (DIAGNOSTIC)."""
import pickle, json, gzip, os
import numpy as np, pandas as pd
import gseapy as gp

TMP = r"${PROJECT_ROOT}\temp"
OUT = r"${PROJECT_ROOT}\output"

# ---- A. stable genes C5 (primary vs sensitivity, both FDR<0.05 same sign)
sens = pickle.load(open(os.path.join(TMP, "b3_sensitivity_results.pkl"), "rb"))
c5s = sens["results"]["C5"].copy()
c5s["gene"] = c5s["gene"].astype(str)
prim = pd.read_csv(os.path.join(OUT, "B3_FULL_RESULTS_C5_INTERACTION.tsv"), sep="\t")
prim["gene"] = prim["gene"].astype(str)
m = c5s[["gene", "log2FC", "padj", "stat"]].merge(
    prim[["gene", "log2FC", "padj", "stat"]], on="gene", suffixes=("_sens", "_prim"))
stable = m[(m["padj_sens"] < 0.05) & (m["padj_prim"] < 0.05) &
           (np.sign(m["log2FC_sens"]) == np.sign(m["log2FC_prim"]))].copy()
stable["min_padj"] = np.minimum(stable["padj_sens"], stable["padj_prim"])
stable["max_padj"] = np.maximum(stable["padj_sens"], stable["padj_prim"])
stable = stable.sort_values("min_padj")
sym = pd.read_csv(os.path.join(TMP, "b3_ensg_symbol_map.tsv"), sep="\t", dtype=str)
sym["g_strip"] = sym["gene_id"].str.split(".").str[0]
c5s["g_strip"] = c5s["gene"].str.split(".").str[0]
smap = dict(zip(sym["g_strip"], sym["symbol"]))
stable["symbol"] = stable["gene"].str.split(".").str[0].map(smap)
stable.to_csv(os.path.join(TMP, "b3_stable_genes_c5.tsv"), sep="\t", index=False)
print("stable C5 genes both FDR<0.05 same sign:", len(stable))
print(stable.head(15)[["gene", "symbol", "log2FC_prim", "log2FC_sens", "min_padj"]].to_string(index=False))

# ---- B. C5 sensitivity ranked GSEA (signed stat, 19-sample)
rnk = c5s.merge(sym, on="g_strip", how="inner")
rnk = rnk.drop_duplicates("symbol")[["symbol", "stat"]].dropna()
rnk = rnk[np.isfinite(rnk["stat"])]
rnk.to_csv(os.path.join(TMP, "b3_rnk_C5_sens.rnk"), sep="\t", index=False, header=False)
print("C5 sens rnk size:", len(rnk))

libs = {
    "HALLMARK": os.path.join(TMP, "b3_gslib_HALLMARK.gmt"),
    "REACTOME": os.path.join(TMP, "b3_gslib_REACTOME.gmt"),
    "GOBP": os.path.join(TMP, "b3_gslib_GOBP.gmt"),
}
prim_g = {}
for k in libs:
    g = pd.read_csv(os.path.join(OUT, f"B3_GSEA_{k}.tsv"), sep="\t")
    g = g[g["contrast"] == "C5"] if "contrast" in g.columns else g
    prim_g[k] = g

rows = []
for lib, gmt in libs.items():
    pre = gp.prerank(rnk=os.path.join(TMP, "b3_rnk_C5_sens.rnk"),
                     gene_sets=gmt, min_size=5, max_size=5000,
                     permutation_num=1000, seed=1, outdir=None, no_plot=True)
    res = pre.res2d.copy()
    res["lib"] = lib
    res["contrast"] = "C5_sens"
    rows.append(res)
sens_gsea = pd.concat(rows, ignore_index=True)
sens_gsea.to_csv(os.path.join(TMP, "b3_gsea_C5_sensitivity_raw.tsv"), sep="\t", index=False)

# intersect FDR<0.05 primary & sensitivity for C5
out_rows = []
for lib in libs:
    p = prim_g[lib].copy()
    s = sens_gsea[sens_gsea["lib"] == lib].copy()
    pp = p.rename(columns=lambda c: c + "_prim") if False else p
    s = s.rename(columns={"Term": "Term", "NES": "NES_sens", "FDR q-val": "FDR_sens",
                          "FWER p-val": "FWER_sens", "Leading edge": "LE_sens", "lib": "lib_s"})
    p = p.rename(columns={"NES": "NES_prim", "FDR q-val": "FDR_prim",
                          "FWER p-val": "FWER_prim", "Leading edge": "LE_prim", "lib": "lib_p"})
    j = p.merge(s, on=["lib", "Term"], how="inner")
    j = j[(j["FDR_prim"] < 0.05) & (j["FDR_sens"] < 0.05)]
    j["dir_agree"] = np.sign(j["NES_prim"]) == np.sign(j["NES_sens"])
    out_rows.append(j)
stab = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame()
stab.to_csv(os.path.join(TMP, "b3_gsea_C5_stable.tsv"), sep="\t", index=False)
summary = {
    "n_prim_sig": {lib: int((prim_g[lib]["FDR q-val"] < 0.05).sum()) for lib in libs},
    "n_sens_sig": {lib: int((sens_gsea[sens_gsea["lib"] == lib]["FDR q-val"] < 0.05).sum()) for lib in libs},
    "n_stable_both_FDR005": {lib: int((stab["lib"] == lib).sum()) for lib in libs} if len(stab) else {k: 0 for k in libs},
    "stable_direction_conflict": int((stab["dir_agree"] == False).sum()) if len(stab) else 0,
}
json.dump(summary, open(os.path.join(TMP, "b3_c5_sens_gsea_summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
if len(stab):
    print(stab.sort_values(["lib", "FDR_prim"]).head(30)[["lib", "Term", "NES_prim", "NES_sens", "FDR_prim", "FDR_sens"]].to_string(index=False))

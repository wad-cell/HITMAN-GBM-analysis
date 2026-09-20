# -*- coding: utf-8 -*-
"""Extract precise GSEA points + leading-edge/stable-gene coverage for reports."""
import json, os
import numpy as np, pandas as pd

OUT = r"${PROJECT_ROOT}\output"
TMP = r"${PROJECT_ROOT}\temp"

libs = ["HALLMARK", "REACTOME", "GOBP"]
allg = {}
for lib in libs:
    df = pd.read_csv(os.path.join(OUT, f"B3_GSEA_{lib}.tsv"), sep="\t")
    df = df[df["contrast"] == "C5"].copy()
    df["FDR"] = df["FDR q-val"]
    allg[lib] = df

out = {}
for lib in libs:
    df = allg[lib]
    sig = df[df["FDR"] < 0.05].copy()
    sig = sig.sort_values("NES", ascending=False)
    out[lib] = {
        "n_sig": int(len(sig)),
        "top_up": [{"t": r["Term"], "nes": float(r["NES"]), "fdr": float(r["FDR"])}
                   for _, r in sig.head(12).iterrows()],
        "top_down": [{"t": r["Term"], "nes": float(r["NES"]), "fdr": float(r["FDR"])}
                     for _, r in sig.tail(12).sort_values("NES").iterrows()],
    }

# C5-only: significant in C5 but not significant in any other contrast (same term, both FDR<0.05 counted; opposite-direction C-sig also disqualifies "C5-only" only if same-direction? Use strict: padj>=0.05 in C1-C4)
other = {}
for lib in libs:
    d = pd.read_csv(os.path.join(OUT, f"B3_GSEA_{lib}.tsv"), sep="\t")
    other[lib] = {c: d[(d["contrast"] == c)][["Term", "NES", "FDR q-val"]]
                  for c in ["C1", "C2", "C3", "C4"]}

c5only = {}
for lib in libs:
    sig = allg[lib][allg[lib]["FDR"] < 0.05]
    rows = []
    for _, r in sig.iterrows():
        sig_else = False
        for c, od in other[lib].items():
            m = od[od["Term"] == r["Term"]]
            if len(m) and (m["FDR q-val"].iloc[0] < 0.05) and (np.sign(m["NES"].iloc[0]) == np.sign(r["NES"])):
                sig_else = True
                break
        if not sig_else:
            rows.append({"t": r["Term"], "nes": float(r["NES"]), "fdr": float(r["FDR"])})
    rows.sort(key=lambda x: x["fdr"])
    c5only[lib] = rows
out["C5_only"] = c5only

# ---- stable genes coverage of C5 significant pathways (descriptive)
st = pd.read_csv(os.path.join(TMP, "b3_stable_genes_c5.tsv"), sep="\t")
stab_sym = set(st["symbol"].dropna())
cover = {}
for lib in libs:
    df = allg[lib][allg[lib]["FDR"] < 0.05]
    rows = []
    for _, r in df.iterrows():
        le = str(r.get("Leading edge", ""))
        genes = [x for x in str(r.get("lead_genes", le)).split(",")] if r.get("lead_genes") is not None else []
        if not genes:
            genes = [x for x in le.replace(";", ",").split(",") if x]
        # fallback: name tokens not usable; coverage on leading edge gene list only if present
        if genes:
            inter = len(stab_sym.intersection(set(genes)))
            rows.append((r["Term"], inter, len(genes)))
    cover[lib] = rows[:10]
out["stable_leading_coverage_sample"] = cover

json.dump(out, open(os.path.join(TMP, "b3_final_gsea_points.json"), "w"), indent=1, default=str)
print("n_sig:", {lib: out[lib]["n_sig"] for lib in libs})
print("C5_only counts:", {lib: len(c5only[lib]) for lib in libs})
for lib in libs:
    print("--", lib, "top_up:", [(x["t"], round(x["nes"], 2)) for x in out[lib]["top_up"][:5]])
    print("   top_down:", [(x["t"], round(x["nes"], 2)) for x in out[lib]["top_down"][:5]])
    print("   only:", [x["t"] for x in c5only[lib]][:12])

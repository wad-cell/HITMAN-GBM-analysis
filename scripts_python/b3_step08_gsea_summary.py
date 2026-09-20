"""B3 - step08: GSEA summary for report (replication anchors, C5-specific)"""
import os, json, traceback
import numpy as np, pandas as pd
BASE = r"${PROJECT_ROOT}"
OUT = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")
LOG = os.path.join(TEMP, "b3_step08_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)
try:
    df = pd.read_csv(os.path.join(OUT, "B3_GSEA_ALL.tsv"), sep="\t")
    df["TermU"] = df["Term"].str.upper()
    summ = {}
    for lib in ["HALLMARK", "REACTOME", "GOBP"]:
        sub = df[df["library"] == lib]
        cnt = sub[sub["FDR q-val"] < 0.05].groupby("contrast").size()
        summ[f"n_sig_{lib}"] = {c: int(cnt.get(c, 0)) for c in ["C1","C2","C3","C4","C5"]}
    # C5 significant term lists per library
    c5 = df[(df["contrast"] == "C5") & (df["FDR q-val"] < 0.05)]
    c5terms = {}
    for lib in ["HALLMARK", "REACTOME", "GOBP"]:
        sub = c5[c5["library"] == lib].sort_values("NES", ascending=False)
        c5terms[lib] = [{"term": r.Term, "NES": round(float(r.NES), 3),
                         "p": round(float(r[5]), 4), "FDR": round(float(r[6]), 4)}
                        for r in sub.itertuples()]
    # overlap of significant term sets across contrasts (per library, Hallmark+Reactome)
    def sig_sets(lib):
        return {c: set(df[(df.library == lib) & (df.contrast == c) & (df["FDR q-val"] < 0.05)]["Term"]) for c in ["C1","C2","C3","C4","C5"]}
    overlap = {}
    for lib in ["HALLMARK", "REACTOME", "GOBP"]:
        ss = sig_sets(lib)
        c5s = ss["C5"]
        overlap[lib] = {
            "C5_only": sorted(c5s - ss["C4"] - ss["C1"] - ss["C2"]),
            "C5_shared_C4": sorted(c5s & ss["C4"]),
            "C5_shared_C1": sorted(c5s & ss["C1"]),
            "C5_shared_C2": sorted(c5s & ss["C2"]),
            "n_C5": len(c5s), "n_C5_only": len(c5s - ss["C4"] - ss["C1"] - ss["C2"]),
        }
    # replication anchors pattern table
    patterns = {
        "UPR/ER_stress": ["UNFOLDED PROTEIN", "ENDOPLASMIC RETICULUM STRESS", "ER-ASSOCIATED", "ERAD", "ATF6", "IRE1", "XBP1"],
        "autophagy": ["AUTOPHAG", "MTORC1", "LYSOSOM"],
        "cell_cycle": ["CELL CYCLE", "G2M", "G2/M", "E2F", "MITOTIC", "S PHASE", "M PHASE"],
        "apoptosis": ["APOPTOS", "P53", "PROGRAMMED CELL DEATH"],
        "cytoskeleton_adhesion": ["CYTOSKELET", "FOCAL ADHESION", "CELL ADHESION", "CADHERIN", "INTEGRIN", "ACTIN"],
        "mitochondrial": ["MITOCHONDRI", "OXIDATIVE PHOSPHORYLATION", "RESPIRATORY ELECTRON"],
        "proliferation": ["MYC TARGETS", "E2F TARGETS", "G2M CHECKPOINT"],
        "inflammation": ["TNFA", "TNF-ALPHA", "NF-KB", "INFLAMMATORY", "INTERFERON", "IL6", "JAK/STAT"],
        "hypoxia": ["HYPOXIA"],
        "glycolysis": ["GLYCOLYSIS"],
    }
    anchor = {}
    for key, pats in patterns.items():
        rows = []
        for _, r in df[df["TermU"].str.contains("|".join(pats), regex=True)].iterrows():
            rows.append({"contrast": r.contrast, "library": r.library, "term": r.Term,
                         "NES": round(float(r.NES), 3), "FDR": round(float(r["FDR q-val"]), 4)})
        anchor[key] = sorted(rows, key=lambda x: (x["contrast"], x["NES"]))
    with open(os.path.join(TEMP, "b3_gsea_summary.json"), "w") as f:
        json.dump({"counts": summ, "c5_terms": c5terms, "overlap": overlap, "anchors": anchor}, f, indent=1, default=str)
    log("counts", summ)
    log("HALLMARK C5 overlap:", overlap["HALLMARK"])
    log("REACTOME nC5_only", len(overlap["REACTOME"]["C5_only"]), overlap["REACTOME"]["C5_only"][:15])
    log("STEP08 OK")
except Exception:
    log("STEP08 EXC", traceback.format_exc())

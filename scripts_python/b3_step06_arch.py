"""B3 - step06: DE summary + effect architecture (B3.8) -> B3_DE_SUMMARY.tsv, B3_EFFECT_ARCHITECTURE.md"""
import os, json, traceback
import numpy as np, pandas as pd
from scipy import stats as ss
BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp")
LOG  = os.path.join(TEMP, "b3_step06_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)
try:
    # load results into single wide frame
    frames = {}
    for c in ["C1","C2","C3","C4","C5"]:
        df = pd.read_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{c}.tsv"), sep="\t").set_index("gene")
        frames[c] = df
    wide = pd.DataFrame(index=frames["C1"].index)
    for c in ["C1","C2","C3","C4","C5"]:
        df = frames[c]
        wide[f"{c}_log2FC"] = df["log2FC"]
        wide[f"{c}_padj"] = df["padj"]
        wide[f"{c}_baseMean"] = df["baseMean"]
    def sig(c):
        return (wide[f"{c}_padj"] < 0.05).fillna(False)
    # ---- DE summary
    de_rows = []
    for c in ["C1","C2","C3","C4","C5"]:
        s = sig(c)
        fc = wide[f"{c}_log2FC"][s]
        de_rows.append({
            "contrast": c,
            "n_FDR005": int(s.sum()),
            "n_up": int((fc > 0).sum()),
            "n_down": int((fc < 0).sum()),
            "median_log2FC_sig": round(float(fc.median()), 4) if len(fc) else "",
            "q25_log2FC_sig": round(float(fc.quantile(0.25)), 4) if len(fc) else "",
            "q75_log2FC_sig": round(float(fc.quantile(0.75)), 4) if len(fc) else "",
            "max_abs_log2FC_sig": round(float(fc.abs().max()), 4) if len(fc) else "",
        })
    de = pd.DataFrame(de_rows)
    de.to_csv(os.path.join(OUT, "B3_DE_SUMMARY.tsv"), sep="\t", index=False)
    log("DE summary saved")

    # ---- effect architecture
    s1, s2, s4, s5 = sig("C1"), sig("C2"), sig("C4"), sig("C5")
    uni = s4 | s5
    n_uni = int(uni.sum())
    log("C4|C5 universe", n_uni)

    class1 = s5 & ~s1 & ~s2                      # interaction only, no marginal
    class2 = s5 & (s1 | s2)                      # interaction with marginal
    class3 = s4 & ~s5 & s1 & ~s2                 # additive field-driven
    class4 = s4 & ~s5 & ~s1 & s2                 # additive man-driven
    class5 = s4 & ~s5 & s1 & s2                  # additive both
    class6 = s4 & ~s5 & ~s1 & ~s2                # C4-only
    labels = {
        "interaction-only (C5 sig, no marginal C1/C2)": class1,
        "interaction with marginal (C5 sig & C1/C2 sig)": class2,
        "additive field-driven (C4 sig, C1 only)": class3,
        "additive MAN-driven (C4 sig, C2 only)": class4,
        "additive both (C4 sig, C1+C2)": class5,
        "C4-only, no marginal & no C5 (low-power flag)": class6,
    }
    counts = {k: int(v.sum()) for k, v in labels.items()}
    log("architecture counts", counts)

    # numeric additive vs interaction evidence
    g = wide.loc[uni]
    c4 = g["C4_log2FC"]; c1 = g["C1_log2FC"]; c2 = g["C2_log2FC"]; c5 = g["C5_log2FC"]
    r_add = ss.spearmanr((c1 + c2), c4)[0]
    r_c5_c4 = ss.spearmanr(c5, c4)[0]
    frac_c5_abs = (c5.abs() / c4.abs().replace(0, np.nan)).dropna()
    median_ratio = float(frac_c5_abs.median())
    n_c5gtc4 = int((c5.abs() > c4.abs()).sum())
    # for C5 sig genes: direction relative to C4
    g_sig5 = wide.loc[class1 | class2]
    same_dir_c4 = (np.sign(g_sig5["C5_log2FC"]) == np.sign(g_sig5["C4_log2FC"])).mean()
    log("numerics", {"spearman_add": r_add, "spearman_c5_c4": r_c5_c4,
                     "median_|C5/C4|": median_ratio, "n_|C5|>|C4|": n_c5gtc4,
                     "C5sig_same_dir_as_C4_frac": float(same_dir_c4)})

    # gene-level table for universe
    arch = pd.DataFrame(index=wide.index)
    for k, v in labels.items():
        arch[k] = v
    arch["C4_log2FC"] = wide["C4_log2FC"]; arch["C5_log2FC"] = wide["C5_log2FC"]
    arch["C1_log2FC"] = wide["C1_log2FC"]; arch["C2_log2FC"] = wide["C2_log2FC"]
    arch = arch[uni]
    arch.to_csv(os.path.join(TEMP, "b3_architecture_genes.tsv"), sep="\t")

    # save numeric summary json for report generation
    summ = {
        "DE": de.to_dict("records"),
        "architecture_counts": counts,
        "spearman_additive_pred_C4": round(float(r_add), 4),
        "spearman_C5_C4": round(float(r_c5_c4), 4),
        "median_abs_C5_over_C4": round(float(median_ratio), 4),
        "n_absC5_gt_absC4": n_c5gtc4,
        "n_universe_C4orC5": n_uni,
        "C5sig_same_dir_as_C4_frac": round(float(same_dir_c4), 4),
        "C4_sig_frac_with_C5sig": round(float((s4 & s5).sum() / max(s4.sum(), 1)), 4),
    }
    with open(os.path.join(TEMP, "b3_arch_summary.json"), "w") as f:
        json.dump(summ, f, indent=2)
    log("STEP06 OK")
except Exception:
    log("STEP06 EXC", traceback.format_exc())

"""
B3 - step04: C5 deep diagnostics (B3.5) + sensitivity comparison (B3.6) + QC figures (DIAGNOSTIC_ONLY)
"""
import os, json, pickle, traceback
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sstats

BASE = r"${PROJECT_ROOT}"
OUT  = os.path.join(BASE, "output"); TEMP = os.path.join(BASE, "temp"); QC = os.path.join(OUT, "B3_QC")
LOG  = os.path.join(TEMP, "b3_step04_log.txt")
def log(*a):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(" ".join(str(x) for x in a) + "\n")
    print(*a, flush=True)

try:
    # ---- load primary per-contrast results
    prim = {}
    for cname in ["C1","C2","C3","C4","C5"]:
        df = pd.read_csv(os.path.join(OUT, f"B3_FULL_RESULTS_{cname}.tsv"), sep="\t")
        df = df.set_index("gene")
        prim[cname] = df
    with open(os.path.join(TEMP, "b3_primary_artifacts.pkl"), "rb") as f:
        art = pickle.load(f)
    cooks_prim = art["cooks"]
    with open(os.path.join(TEMP, "b3_sensitivity_results.pkl"), "rb") as f:
        sens_pkg = pickle.load(f)
    sens = sens_pkg["results"]
    cooks_sens = sens_pkg["cooks"]
    g = prim["C5"]
    gs = sens["C5"]

    out_diag = {}

    # ================= B3.5 C5 interaction deep check =================
    d = g["log2FC"].dropna()
    d_sig = g.loc[g["padj"] < 0.05, "log2FC"].dropna()
    diag35 = {
        "n_genes_with_stat": int(g["pvalue"].notna().sum()),
        "n_FDR005": int((g["padj"] < 0.05).sum()),
        "log2FC_median": float(d.median()),
        "log2FC_mad": float((d - d.median()).abs().median()),
        "log2FC_pct25": float(d.quantile(0.25)),
        "log2FC_pct75": float(d.quantile(0.75)),
        "log2FC_sig_median": float(d_sig.median()) if len(d_sig) else None,
        "log2FC_sig_pct25": float(d_sig.quantile(0.25)) if len(d_sig) else None,
        "log2FC_sig_pct75": float(d_sig.quantile(0.75)) if len(d_sig) else None,
        "max_abs_log2FC": float(d.abs().max()),
        "pct_sig_abs_ge_1": float((d_sig.abs() >= 1).mean()) if len(d_sig) else None,
    }
    log("C5 diag", diag35)

    # strongest genes
    sig = g[g["padj"] < 0.05].copy()
    top_neg = sig.nsmallest(30, "log2FC")
    top_pos = sig.nlargest(30, "log2FC")
    top_by_p = sig.nsmallest(30, "padj")
    top_genes = pd.concat([
        top_by_p.assign(tier="top30_by_padj"),
        top_pos.assign(tier="top30_pos"),
        top_neg.assign(tier="top30_neg"),
    ])
    top_genes = top_genes[["tier","baseMean","log2FC","lfcSE","stat","pvalue","padj"]]
    top_genes.to_csv(os.path.join(OUT, "B3_INTERACTION_TOP_GENES.tsv"), sep="\t")
    log("top genes saved", top_genes.shape)

    # MA / volcano / hist / C5-vs-C4
    fig, ax = plt.subplots(1, 3, figsize=(18,5))
    ax[0].scatter(np.log10(g["baseMean"]+1), g["log2FC"], s=2, alpha=0.3, color="grey")
    m = g["padj"] < 0.05
    ax[0].scatter(np.log10(g.loc[m,"baseMean"]+1), g.loc[m,"log2FC"], s=2, alpha=0.6, color="crimson")
    ax[0].set_title("C5 MA plot (DIAGNOSTIC_ONLY)")
    ax[0].set_xlabel("log10(baseMean+1)"); ax[0].set_ylabel("C5 log2FC")
    ax[1].scatter(g["log2FC"], -np.log10(g["pvalue"]+1e-300), s=2, alpha=0.3, color="grey")
    ax[1].scatter(g.loc[m,"log2FC"], -np.log10(g.loc[m,"pvalue"]+1e-300), s=2, alpha=0.6, color="crimson")
    ax[1].set_title("C5 volcano (DIAGNOSTIC_ONLY)"); ax[1].set_xlabel("C5 log2FC"); ax[1].set_ylabel("-log10 P")
    ax[2].hist(d.clip(-4,4), bins=80, color="steelblue", alpha=0.7)
    ax[2].axvline(0, color="black", lw=0.8); ax[2].set_title("C5 interaction log2FC dist (DIAGNOSTIC_ONLY)")
    plt.tight_layout(); plt.savefig(os.path.join(QC, "B3_QC_C5_diagnostics.png"), dpi=130); plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(12,5))
    c4 = prim["C4"]["log2FC"]
    common = g.index.intersection(c4.dropna().index)
    ax[0].scatter(c4.loc[common], g.loc[common,"log2FC"], s=2, alpha=0.25, color="grey")
    msig = g.loc[common, "padj"] < 0.05
    ax[0].scatter(c4.loc[common][msig], g.loc[common,"log2FC"][msig], s=2, alpha=0.5, color="crimson")
    ax[0].set_title("C5 log2FC vs C4 log2FC (DIAGNOSTIC_ONLY)")
    ax[0].set_xlabel("C4 log2FC (total)"); ax[0].set_ylabel("C5 log2FC (interaction)")
    # additive expectation: if additive, C4 gene effect approx C1+C2 (when no interaction). show C4 vs (C1+C2)
    c1 = prim["C1"]["log2FC"]; c2 = prim["C2"]["log2FC"]
    add_pred = (c1 + c2).reindex(common)
    ax[1].scatter(add_pred, c4.loc[common], s=2, alpha=0.25, color="grey")
    ax[1].plot([-6,6],[-6,6], color="black", lw=0.8)
    ax[1].set_title("C4 vs C1+C2 (additive prediction) (DIAGNOSTIC_ONLY)")
    ax[1].set_xlabel("C1+C2 log2FC"); ax[1].set_ylabel("C4 log2FC")
    plt.tight_layout(); plt.savefig(os.path.join(QC, "B3_QC_C5_vs_C4_additivity.png"), dpi=130); plt.close()

    # Cook's distance presence (primary model)
    try:
        cooks_mat = None
        if "_mu_LFC" in art["dds"].obsm:
            art["dds"].calculate_cooks()
            cooks_mat = art["dds"].layers.get("cooks", None)
        if cooks_mat is not None:
            sample_names = pd.read_pickle(os.path.join(TEMP,"b3_meta.pkl"))["sample"].values
            cooks_df = pd.DataFrame(cooks_mat, index=sample_names)
            # transpose to genes x samples for intersect with g index (keep if shape matches)
            if cooks_mat.shape[1] == g.shape[0]:
                cooks_df = pd.DataFrame(cooks_mat.T, index=g.index)
            # F-distribution threshold approximation
            import scipy.stats as fs
            k = 4  # design rank
            n = cooks_df.shape[1]
            thr = fs.f.ppf(0.99, k, n - k)
            maxc = cooks_df.max(axis=1)
            n_cooks_any = int((maxc > thr).sum())
            sig_genes = g.index[m]
            n_sig_cook = int(maxc.loc[sig_genes].gt(thr).sum())
            cooks_summary = {"cooks_matrix_shape": list(cooks_df.shape),
                             "f_threshold_99": round(float(thr),4),
                             "n_genes_cooks_above_thr": n_cooks_any,
                             "n_C5sig_cooks_above_thr": n_sig_cook}
            # per-sample contributions among C5 sig genes exceeding threshold
            sample_hits = (cooks_df.loc[sig_genes] > thr).sum(axis=0).sort_values(ascending=False)
            cooks_summary["sample_hits_C5sig"] = {str(s): int(v) for s, v in sample_hits.items()}
            # refit/replaced genes info
            v = art["dds"].var
            cooks_summary["n_genes_replaced_refit"] = int(v["replaced"].sum()) if "replaced" in v else None
        else:
            cooks_summary = {"note": "cooks unavailable (low_memory or missing mu_LFC)"}
    except Exception as e:
        cooks_summary = {"note": f"cooks computation failed: {e}"}
    log("cooks summary", cooks_summary)
    diag35["cooks"] = cooks_summary

    # low-count dependence: fraction of C5 sig by baseMean quartile
    qs = pd.qcut(g["baseMean"].rank(method="first"), 4, labels=["Q1_low","Q2","Q3","Q4_high"])
    lowdep = pd.DataFrame({"q": qs, "sig": g["padj"] < 0.05}).groupby("q", observed=True)["sig"].agg(["sum","count"])
    lowdep["frac"] = lowdep["sum"]/lowdep["count"]
    diag35["lowcount_dependence_by_baseMean_quartile"] = lowdep.round(4).to_dict("index")

    # relation C5 sig vs C1/C2/C3/C4 sig & direction
    rel = {}
    for cc, name in [("C1","C1_sig"),("C2","C2_sig"),("C3","C3_sig"),("C4","C4_sig")]:
        cs = prim[cc]["padj"] < 0.05
        cs = cs.reindex(g.index).fillna(False)
        both = cs & m
        rel[name] = {
            "n_C5sig_also_that": int(both.sum()),
            "frac_of_C5sig": float(both.sum()/max(m.sum(),1)),
            "sign_concord_C5_same_dir": float((np.sign(g.loc[both,"log2FC"])==np.sign(prim[cc].loc[both,"log2FC"])).mean()) if both.sum() else None,
        }
    diag35["relation_C5sig_vs_other_contrasts"] = rel

    # ================= B3.6 sensitivity comparison (C5 focus + all) =================
    def compare(pri, sen, name):
        if "gene" in sen.columns:
            sen = sen.set_index("gene")
        j = pri.index.intersection(sen.index)
        x = pri.loc[j,"log2FC"].dropna(); y = sen.loc[j,"log2FC"].dropna()
        j2 = x.index.intersection(y.index)
        r_pearson = sstats.pearsonr(x.loc[j2], y.loc[j2])[0]
        r_spearman = sstats.spearmanr(x.loc[j2], y.loc[j2])[0]
        sc = (np.sign(x.loc[j2]) == np.sign(y.loc[j2])).mean()
        xp = pri.loc[j,"pvalue"]; yp = sen.loc[j,"pvalue"]
        # stat correlation (Spearman primary; Pearson is unstable due to extreme |stat| genes)
        xs = pri.loc[j2,"stat"]; ys = sen.loc[j2,"stat"]
        r_stat_spearman = sstats.spearmanr(xs, ys)[0]
        # Pearson trimmed to central 98% for reference
        q = xs.quantile([0.01, 0.99]); mask = (xs >= q.iloc[0]) & (xs <= q.iloc[1]) & ys.notna()
        r_stat_pearson_trimmed = sstats.pearsonr(xs[mask], ys[mask])[0] if mask.sum() > 2 else None
        # top overlap by padj
        def topk(padj_series, k):
            return set(padj_series.dropna().nsmallest(k).index)
        topres = {}
        for k in [50,100,250]:
            a = topk(pri.loc[j,"padj"], k); b = topk(sen.loc[j,"padj"], k)
            topres[f"top{k}"] = {"n_union": len(a|b), "n_intersection": len(a&b),
                                 "jaccard": round(len(a&b)/max(len(a|b),1),4)}
        # FDR sets
        a_sig = set(pri.index[pri["padj"] < 0.05]); b_sig = set(sen.index[sen["padj"] < 0.05])
        fdr = {"n_primary": len(a_sig), "n_sensitivity": len(b_sig),
               "inter": len(a_sig & b_sig), "union": len(a_sig | b_sig),
               "jaccard": round(len(a_sig & b_sig)/max(len(a_sig | b_sig),1),4)}
        return {"name": name, "n_joint": len(j2), "pearson_log2FC": round(float(r_pearson),4),
                "spearman_log2FC": round(float(r_spearman),4),
                "spearman_stat": round(float(r_stat_spearman),4),
                "pearson_stat_trimmed_98pct": (round(float(r_stat_pearson_trimmed),4) if r_stat_pearson_trimmed is not None else None),
                "sign_concordance": round(float(sc),4),
                "top_overlap": topres, "fdr_set": fdr}

    comp_all = {cname: compare(prim[cname], sens[cname], cname) for cname in ["C1","C2","C3","C4","C5"]}
    c5comp = comp_all["C5"]
    # verdict thresholds
    sc = c5comp["sign_concordance"]; jac = c5comp["fdr_set"]["jaccard"]
    t250 = c5comp["top_overlap"]["top250"]["jaccard"]
    if sc >= 0.95 and jac >= 0.6 and t250 >= 0.7:
        verdict = "ROBUST"
    elif sc >= 0.85 and jac >= 0.3:
        verdict = "PARTIALLY ROBUST"
    else:
        verdict = "OUTLIER-SENSITIVE"
    log("sensitivity C5", c5comp, "verdict", verdict)

    with open(os.path.join(TEMP, "b3_diag_json.json"), "w") as f:
        json.dump({"B3.5": diag35, "B3.6": comp_all, "verdict": verdict}, f, indent=2, default=str)
    log("STEP04 OK")
except Exception as e:
    log("STEP04 EXC", traceback.format_exc())

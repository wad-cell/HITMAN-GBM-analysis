# -*- coding: utf-8 -*-
"""C3 main analysis: Q1 correlation architecture, Q2/Q3 subtype context, C5 exploratory.
Frozen inputs only: c3_data/tcga_ssgsea_wide.tsv, cgga_ssgsea_wide.tsv,
cgga_clinical_iv.tsv, IDH1-Sample-Status.txt (TCGA official Verhaak subtype labels, n=44 overlap).
No signature/direction/scoring change.
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C3D = os.path.join(BASE, "temp", "c3_data")
OUT = os.path.join(BASE, "output")

A8 = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
      "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
B2 = ["M_shared_ER_UPR", "M_autophagy_lysosome"]  # grade B at module-score level (component caveat)
C5S = ["C5_UP", "C5_DOWN", "C5_STABLE_UP", "C5_STABLE_DOWN"]

# Pre-registered expected signs among A-level module pairs (Q1), frozen before results.
EXPECTED = {
 ("M_cholesterol_SREBP", "M_TNF_NFkB"): 1,
 ("M_cholesterol_SREBP", "M_p53_apoptosis"): 1,
 ("M_TNF_NFkB", "M_p53_apoptosis"): 1,
 ("C4_UP", "M_cholesterol_SREBP"): 1,
 ("C4_UP", "M_TNF_NFkB"): 1,
 ("C4_UP", "M_p53_apoptosis"): 1,
 ("M_shared_antiproliferative", "M_shared_ER_UPR"): 1,
 ("C4_DOWN", "M_shared_antiproliferative"): 1,
 ("C4_DOWN", "M_shared_ER_UPR"): 1,
 ("C4_UP", "C4_DOWN"): -1,
 ("M_shared_antiproliferative", "M_cholesterol_SREBP"): -1,
 ("M_shared_antiproliferative", "M_TNF_NFkB"): -1,
 ("M_shared_antiproliferative", "M_p53_apoptosis"): -1,
 ("M_shared_ER_UPR", "M_cholesterol_SREBP"): -1,
 ("M_shared_ER_UPR", "M_TNF_NFkB"): -1,
 ("M_shared_ER_UPR", "M_p53_apoptosis"): -1,
 ("C4_DOWN", "M_cholesterol_SREBP"): -1,
 ("C4_DOWN", "M_TNF_NFkB"): -1,
 ("C4_DOWN", "M_p53_apoptosis"): -1,
 ("C4_UP", "M_shared_antiproliferative"): -1,
 ("C4_UP", "M_shared_ER_UPR"): -1,
}
# RIG-I(AC) has no bulk AC proxy in the frozen expected-pair set (reported descriptively).


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate(ranked * n / (np.arange(1, n + 1)))
    out = np.empty(n)
    out[order] = np.minimum(adj[::-1][::-1], 1.0)
    return out


def load_scores(tag):
    df = pd.read_csv(os.path.join(C3D, f"{tag}_ssgsea_wide.tsv"), sep="\t")
    df = df.set_index("Name")
    return df


def q1_cohort(df, cohort):
    rows = []
    for i, a in enumerate(A8):
        for b in A8[i + 1:]:
            rho, p = stats.spearmanr(df[a], df[b], nan_policy="omit")
            rows.append({"cohort": cohort, "module_a": a, "module_b": b, "spearman_rho": rho, "p": p})
    res = pd.DataFrame(rows)
    res["FDR_BH"] = bh(res["p"].values)
    exp = np.array([EXPECTED.get((r.module_a, r.module_b), EXPECTED.get((r.module_b, r.module_a), np.nan))
                    for r in res.itertuples()])
    res["expected_sign"] = exp
    res["direction_match"] = np.sign(res["spearman_rho"]) == exp
    res.loc[exp == 0, "direction_match"] = np.nan
    return res


def cohort_summary(res):
    paired = res[res["expected_sign"].notna()]
    hit = paired["direction_match"].sum()
    tot = paired["direction_match"].notna().sum()
    hit_fdr = paired[(paired["direction_match"]) & (paired["FDR_BH"] < 0.05)]
    # anti-expected significant = formal counter-evidence
    counter = paired[(~paired["direction_match"]) & (paired["FDR_BH"] < 0.05)]
    return {"paired_n": int(tot), "dir_match_n": int(hit),
            "dir_match_pct": round(hit / tot * 100, 1),
            "match_and_fdr05_n": int(len(hit_fdr)),
            "counter_sig_n": int(len(counter))}


def rank_group_test(df, score_cols, group):
    """group: pd.Series (index=sample) with 1/0 and NaN for excluded; test only labeled samples."""
    rows = []
    lab = group[group.notna() & group.index.isin(df.index)]
    for col in score_cols:
        x = df.loc[lab.index, col]
        g1 = lab.values == 1
        if g1.sum() < 3 or (~g1).sum() < 3:
            rows.append({"score": col, "n_group1": int(g1.sum()), "n_group0": int((~g1).sum()),
                         "median_diff_1minus0": np.nan, "mw_p": np.nan, "CI95_low": np.nan, "CI95_high": np.nan})
            continue
        v1, v0 = x[g1], x[~g1]
        delta = np.median(v1) - np.median(v0)
        rng = np.random.default_rng(20260909)
        b1 = rng.choice(v1.values, size=(5000, len(v1)), replace=True)
        b0 = rng.choice(v0.values, size=(5000, len(v0)), replace=True)
        boot_delta = np.median(b1, axis=1) - np.median(b0, axis=1)
        lo, hi = np.percentile(boot_delta, [2.5, 97.5])
        u, p = stats.mannwhitneyu(v1, v0, alternative="two-sided")
        rows.append({"score": col, "n_group1": int(g1.sum()), "n_group0": int((~g1).sum()),
                     "median_diff_1minus0": delta, "mw_p": p, "CI95_low": lo, "CI95_high": hi})
    return pd.DataFrame(rows)


def main():
    t = load_scores("tcga")
    g = load_scores("cgga")

    # ---- Q1 ----
    q1_t = q1_cohort(t, "TCGA")
    q1_g = q1_cohort(g, "CGGA")
    q1 = pd.concat([q1_t, q1_g], ignore_index=True)
    q1.to_csv(os.path.join(OUT, "C3_PRIMARY_RESULTS.tsv"), sep="\t", index=False)
    print("Q1 TCGA summary", cohort_summary(q1_t))
    print("Q1 CGGA summary", cohort_summary(q1_g))
    paired = q1[q1["expected_sign"].notna()]
    cross = paired.pivot_table(index=["module_a", "module_b"], columns="cohort",
                               values=["spearman_rho", "FDR_BH"])
    print("Q1 cross-cohort direction-consistent pairs:",
          ((np.sign(cross[("spearman_rho", "TCGA")]) == np.sign(cross[("spearman_rho", "CGGA")])) &
           (cross[("spearman_rho", "TCGA")].notna())).sum(), "/", len(cross))

    # ---- Q2 TCGA official subtype (MES vs non-MES), frozen Grade A anchors ----
    idh = pd.read_csv(os.path.join(C3D, "IDH1-Sample-Status.txt"), sep="\t")
    tmeta = pd.DataFrame({"Name": t.index})
    tmeta["CLID"] = tmeta["Name"]
    ann = tmeta.merge(idh, on="CLID", how="left")
    ann["is_MES"] = np.where(ann["Subtype"] == "MES", 1, np.where(ann["Subtype"].notna(), 0, np.nan))
    ann["is_CL"] = np.where(ann["Subtype"] == "CL", 1, np.where(ann["Subtype"].notna(), 0, np.nan))
    sub_ok = ann["Subtype"].notna()
    print("Q2 TCGA official subtype available:", int(sub_ok.sum()),
          ann.loc[sub_ok, "Subtype"].value_counts().to_dict())

    q2_cols = ["C4_UP", "C4_DOWN", "M_shared_antiproliferative", "M_shared_ER_UPR",
               "M_cholesterol_SREBP", "M_TNF_NFkB", "M_p53_apoptosis", "M_RIGI_typeI_IFN"]
    q2_mes = rank_group_test(t, q2_cols, ann.set_index("Name")["is_MES"])
    q2_mes["context"] = "TCGA_official_MES_vs_nonMES"
    q2_mes["FDR_BH"] = bh(q2_mes["mw_p"].fillna(1).values)
    q2_cl = rank_group_test(t, ["M_shared_antiproliferative", "M_shared_ER_UPR", "C4_DOWN"],
                            ann.set_index("Name")["is_CL"])
    q2_cl["context"] = "TCGA_official_CL_vs_nonCL(descriptive)"
    q2_cl["FDR_BH"] = bh(q2_cl["mw_p"].fillna(1).values)
    q2 = pd.concat([q2_mes, q2_cl], ignore_index=True)
    q2.to_csv(os.path.join(OUT, "C3_SECONDARY_RESULTS.tsv"), sep="\t", index=False)  # provisional name; will separate
    # Q2 primary table actually PRIMARY-level per pre-reg; file naming handled later.
    q2_mes.to_csv(os.path.join(C3D, "_q2_mes_intermediate.tsv"), sep="\t", index=False)
    print("Q2 MES table:")
    print(q2_mes.round(4).to_string())

    # ---- Q3 secondary (Grade B, module-score level; component-level NOT testable in bulk) ----
    q3 = rank_group_test(t, B2, ann.set_index("Name")["is_MES"])
    q3["context"] = "TCGA_official_MES_vs_nonMES"
    q3["FDR_BH"] = bh(q3["mw_p"].fillna(1).values)
    print("Q3 table:")
    print(q3.round(4).to_string())

    # ---- C5 exploratory (TCGA MES context; CGGA no subtype so described separately) ----
    c5 = rank_group_test(t, C5S, ann.set_index("Name")["is_MES"])
    c5["context"] = "TCGA_official_MES_vs_nonMES"
    c5["FDR_BH"] = bh(c5["mw_p"].fillna(1).values)
    print("C5 exploratory (TCGA MES):")
    print(c5.round(4).to_string())

    # C5 association with MES-axis proxy (C4_UP score) in BOTH cohorts (exploratory, non-inferential)
    c5ax_rows = []
    for tag, df in [("TCGA", t), ("CGGA", g)]:
        for c in C5S:
            rho, p = stats.spearmanr(df["C4_UP"], df[c], nan_policy="omit")
            c5ax_rows.append({"cohort": tag, "score": c, "C4UP_rho": rho, "p": p})
    c5ax = pd.DataFrame(c5ax_rows)
    c5ax["FDR_BH"] = bh(c5ax["p"].values)
    print("C5 vs C4_UP axis (both cohorts):")
    print(c5ax.round(4).to_string())

    # persist C5 exploratory + C5-axis
    c5.to_csv(os.path.join(OUT, "C3_C5_EXPLORATORY.tsv"), sep="\t", index=False)
    c5ax.to_csv(os.path.join(C3D, "_c5_axis_intermediate.tsv"), sep="\t", index=False)

    # Q3 secondary persist
    q3.to_csv(os.path.join(C3D, "_q3_intermediate.tsv"), sep="\t", index=False)
    print("done")


if __name__ == "__main__":
    main()

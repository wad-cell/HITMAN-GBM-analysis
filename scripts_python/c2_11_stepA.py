# C2.11 Step A: assemble Layer B pool + QC; compute cell-level AUC for LA tissue + LB pool on full Couturier universe
import os, time, numpy as np, pandas as pd, scipy.sparse as sp, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scRNA_lib import cell_auc_for_sigs

T = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(T, "scRNA", "3CA_Couturier", "Data_Couturier2020_Brain")
OUT = os.path.join(os.path.dirname(T), "output")
O = OUT
META = os.path.join(T, "couturier_meta_layers.tsv")
genes = [l.strip() for l in open(os.path.join(B, "Genes.txt"))]
meta = pd.read_csv(META, sep="\t")
idx_la = np.load(os.path.join(T, "couturier_layerA_idx.npy"))
la = meta.iloc[idx_la].copy()
LA_PATS = set(la["patient"].dropna().astype(str))
tissue = la[la["patient"].astype(str).isin(LA_PATS)].copy()   # 17,271 patient-attributed LA cells
print("LA patient-attributed cells:", len(tissue), "patients:", len(LA_PATS))

# ---------- 1) QC from per-sample CNV runs ----------
LB_SAMPLES = ["BT322","BT338_1of2","BT338_2of2","BT363_1of2","BT363_2of2",
              "BT364_2of2","BT390","BT397_1of2","BT397_2of2","BT407"]
LA_SAMPLES = set(tissue["sample"].unique())
rows = []
cnv_calls = []
for s in LB_SAMPLES:
    f = os.path.join(T, "layerb_scores", f"{s}_cnv.tsv")
    d = pd.read_csv(f, sep="\t")
    tum = d[d.cnv_group == "tumor"].copy()
    n_ref = int((d.cnv_group == "normal").sum())
    refs = d[d.cnv_group == "normal"]["cnv_score"]
    thr = float(refs.mean() + 3 * refs.std())
    n_mal = int((tum.lb_call == "malignant_cnv").sum())
    n_dip = int((tum.lb_call == "nonmalignant_diploid").sum())
    off = int((tum["mal_annot"] == True).sum()) if "mal_annot" in tum else int(tum["mal_annot"].sum()) if tum["mal_annot"].dtype.kind in "bi" else 0
    rows.append(dict(patient=str(meta.loc[meta["sample"] == s, "patient"].iloc[0]),
                     sample=s, ref_cells=n_ref, ref_mean=round(float(refs.mean()),5),
                     ref_sd=round(float(refs.std()),5), cnv_method="infercnvpy.tl.infercnv",
                     window_size=100, step=10, dynamic_threshold=1.5, lfc_clip=3,
                     threshold_mean_plus_3sd=round(thr,5),
                     cnv_malignant_n=n_mal, cnv_diploid_n=n_dip))
    c = tum[tum.lb_call == "malignant_cnv"].copy()
    c["source"] = "CNV_inferred"
    cnv_calls.append(c)
qc = pd.DataFrame(rows)
# official inherited small samples (not in Layer A sample list)
small = ["BT322","BT338_2of2","BT364_2of2","BT407"]
off_df = meta[(meta["sample"].isin(small)) & (meta["mal_annot"] == True)].copy()
off_df["lb_call"] = "inherited_official_malignant"; off_df["cnv_score"] = np.nan
off_df["cnv_group"] = "tumor"; off_df["source"] = "official_annotation"
# append official counts to QC
offc = off_df.groupby("sample").size()
qc["official_inherited_n"] = qc["sample"].map(offc).fillna(0).astype(int)
qc["pool_malignant_n"] = qc["cnv_malignant_n"] + qc["official_inherited_n"]
qc["note"] = np.where(qc["sample"].isin(small),
                      "CNV on annotation-missing cells; official malignant inherited",
                      "CNV on annotation-missing cells; official malignant cells already in Layer A (excluded from Layer B pool)")
qc.to_csv(os.path.join(O, "C2_LAYERB_CNV_QC.tsv"), sep="\t", index=False)
print("QC written; total cnv calls:", sum(qc.cnv_malignant_n), "official:", int(offc.sum()))

# ---------- 2) Layer B pool indices ----------
cnv_all = pd.concat(cnv_calls, ignore_index=True)
cnv_cell = cnv_all["cell_name"].tolist()
off_cell = off_df["cell_name"].tolist()
lb_pool_meta = pd.concat([cnv_all, off_df], ignore_index=True)
name2pos = {c: i for i, c in enumerate(meta["cell_name"])}
lb_pos = np.array([name2pos[x] for x in lb_pool_meta["cell_name"]], dtype=np.int64)
print("LB pool cells:", len(lb_pos), "patients:", sorted(meta.loc[lb_pos, "patient"].dropna().unique().tolist()))

# ---------- 3) signature variants (base + CCDEP) ----------
sig = {}
for l in open(os.path.join(O, "C2_2_FROZEN_SIGNATURE_GENES.tsv"), encoding="utf-8", errors="replace"):
    if l.startswith("signature"): continue
    p = l.rstrip("\n").split("\t")
    if len(p) >= 5 and p[4].strip():
        sig.setdefault(p[0].strip(), []).append(p[4].strip())
removal = set()
for l in open(os.path.join(O, "C2_3_CELLCYCLE_REMOVAL_LIST_FROZEN.tsv"), encoding="utf-8", errors="replace"):
    if l.startswith("symbol"): continue
    p = l.rstrip("\n").split("\t")
    removal.add(p[0].strip())
core_sigs = [k for k in sig if not k.endswith("_CCDEP")]
variants = {}
for k in core_sigs:
    variants[k] = sig[k]
    variants[k + "_CCDEP"] = [g for g in sig[k] if g not in removal]
print("sig variants:", len(variants))

# ---------- 4) scores on raw full-universe matrix ----------
t0 = time.time()
M = sp.load_npz(os.path.join(T, "couturier_matrix_csr.npz"))
print("loaded full M", M.shape, f"{time.time()-t0:.0f}s")
all_pos = np.concatenate([tissue.index.values.astype(np.int64), lb_pos])
uniq_pos, back = np.unique(all_pos, return_inverse=True)
info = meta.iloc[uniq_pos].copy()
info["pool"] = np.where(np.isin(uniq_pos, tissue.index.values), "LA_tissue",
                        np.where(np.isin(uniq_pos, lb_pos), "LB", "both"))
Ml = M.tocsc()[:, uniq_pos]
sp.save_npz(os.path.join(T, "c2_11_pool_csc.npz"), Ml)
auc = cell_auc_for_sigs(Ml, genes, variants, already_csc=True)
auc.columns = [c + ("" if c in core_sigs else "") for c in auc.columns]
np.save(os.path.join(T, "c2_11_pool_pos.npy"), uniq_pos)
np.save(os.path.join(T, "c2_11_pool_lb.npy"), lb_pos)
auc.to_csv(os.path.join(T, "c2_11_pool_auc.tsv"), sep="\t")
info.to_csv(os.path.join(T, "c2_11_pool_info.tsv"), sep="\t", index=False)
print("auc", auc.shape, f"{time.time()-t0:.0f}s")
# validation: recomputed vs existing core file (first rows approx)
old = pd.read_csv(os.path.join(T, "couturier_layerA_auc_core.tsv"), sep="\t", nrows=500)
print("old core head auc sample:", old.iloc[0, 0])
print("saved AUC head row C5_UP:", auc.iloc[0]["C5_UP"])
print("done step A")

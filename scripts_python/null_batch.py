import sys, os, time, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scRNA_lib import cell_auc_for_sigs
T = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(T, "..", "output")
sig_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(T, "null_sigs_A.txt")
out_tag = os.path.basename(sig_file).replace(".txt", "")
sigs_to_run = [l.strip() for l in open(sig_file) if l.strip()]
ss2 = pd.read_csv(os.path.join(T, "scRNA", "3CA_Neftel", "Data_Neftel2019_Brain", "SmartSeq2", "Cells.csv"))
genes = [l.strip() for l in open(os.path.join(T, "scRNA", "3CA_Neftel", "Data_Neftel2019_Brain", "SmartSeq2", "Genes.txt"))]
M_csc = sp.load_npz(os.path.join(T, "ss2_matrix_csr.npz")).tocsc()
am = ss2[(ss2.cell_type.eq("Malignant")) & (ss2.age_group.eq("Adult"))].copy()
am["MES"] = am[["MESlike1", "MESlike2"]].max(axis=1); am["NPC"] = am[["NPClike1", "NPClike2"]].max(axis=1)
def lbl(r):
    if r.cell_cycle_phase in ("G1/S", "G2/M"): return "Cycling"
    ax = {"MES": r.MES, "AC": r.AClike, "OPC": r.OPClike, "NPC": r.NPC}
    return max(ax, key=ax.get)
am["state6"] = [lbl(r) for r in am.itertuples()]
states = ["MES", "AC", "OPC", "NPC", "Cycling"]
cell_idx = am.index.values
# load signature genes
sig_map = {}
for l in open(os.path.join(O, "C2_2_FROZEN_SIGNATURE_GENES.tsv"), encoding="utf-8", errors="replace"):
    if l.startswith("signature"): continue
    p = l.rstrip("\n").split("\t")
    if len(p) >= 5:
        sig_map.setdefault(p[0].strip(), []).append(p[4].strip())
rem = set()
for l in open(os.path.join(O, "C2_3_CELLCYCLE_REMOVAL_LIST_FROZEN.tsv"), encoding="utf-8", errors="replace"):
    if l.startswith("symbol"): continue
    p = l.rstrip("\n").split("\t")
    if len(p) >= 1 and p[0].strip():
        rem.add(p[0].strip())
def make_variants_local():
    out = {}
    for k, syms in sig_map.items():
        out[k] = syms
        dep = [s for s in syms if s not in rem]
        out[k + "_CCDEP"] = dep
    return out
variants = make_variants_local()
gi = {g: i for i, g in enumerate(genes)}
expr_mean = np.asarray(M_csc.tocsr()[:, cell_idx].mean(axis=1)).ravel()
bq = pd.qcut(expr_mean, q=10, labels=False, duplicates="drop")
bins = np.array(bq, dtype=float); bins[np.isnan(bins)] = 0; bins = bins.astype(int)
cand = np.flatnonzero(expr_mean > 0); cand_b = bins[cand]
# per-state cells (cap 30 per patient for null), pooled per state
state_pool = {}
for st in states:
    idxs = []
    for p, g in am[am.state6.eq(st)].groupby("sample").groups.items():
        a = np.asarray(g)
        if len(a) > 30: a = np.random.default_rng(1).choice(a, 30, replace=False)
        idxs.append(a)
    state_pool[st] = np.concatenate(idxs) if idxs else np.array([], dtype=int)
rng = np.random.default_rng(123)
rows = []
t0 = time.time()
for k in sigs_to_run:
    base = k.replace("_CCDEP", "")
    syms = variants[k]
    idxs = np.array([gi[g] for g in syms if g in gi], dtype=np.int64)
    sigb = pd.Series(bins[idxs]).value_counts()
    obs = {st: float(np.mean(cell_auc_for_sigs(M_csc, genes, {k: syms}, cell_idx=state_pool[st], already_csc=True)[k].values)) if len(state_pool[st]) else np.nan for st in states}
    nset = {}
    for r in range(20):
        sel = []
        for b in range(10):
            pool = np.flatnonzero(cand_b == b); need = int(sigb.get(b, 0))
            if need > 0 and len(pool) > 0: sel.append(rng.choice(pool, size=min(need, len(pool)), replace=False))
        sel = np.concatenate(sel) if sel else np.array([], dtype=int)
        if len(sel) < len(idxs):
            sel = np.concatenate([sel, rng.choice(np.setdiff1d(cand, sel), size=len(idxs) - len(sel), replace=False)])
        gset = list(np.array(genes)[sel])
        for st in states:
            if len(state_pool[st]) == 0: continue
            a = cell_auc_for_sigs(M_csc, genes, {k: gset}, cell_idx=state_pool[st], already_csc=True)[k].values
            nset.setdefault(st, []).append(float(np.mean(a)))
    for st in states:
        null = np.array(nset.get(st, [])); o = obs.get(st, np.nan)
        emp = (np.sum(null >= o) + 1) / (len(null) + 1) if len(null) else np.nan
        z = (o - null.mean()) / (null.std() if len(null) and null.std() > 0 else np.nan)
        rows.append({"sig": k, "state": st, "n_sig_genes": len(idxs), "obs": o, "null_mean": null.mean() if len(null) else np.nan, "null_sd": null.std() if len(null) else np.nan, "z": z, "emp_p": emp})
    print(k, "done", round(time.time() - t0, 1), "s", flush=True)
ndf = pd.DataFrame(rows)
ndf.to_csv(os.path.join(T, f"ss2_cellnull_{out_tag}.tsv"), sep="\t", index=False)
print("saved", ndf.shape)

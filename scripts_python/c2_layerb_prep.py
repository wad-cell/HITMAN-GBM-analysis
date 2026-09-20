# C2.11 Layer B prep: build per-sample anndata scaffold + infercnv run (prototype on BT407)
import os, time, numpy as np, pandas as pd, scipy.sparse as sp, sys
T = r"${PROJECT_ROOT}\temp"
sys.path.insert(0, T)

def load_full_csc():
    # cache full csc filtered to mapped genes if not present
    fcache = os.path.join(T, "couturier_mapped_csc.npz")
    genes = [l.strip() for l in open(os.path.join(T, "scRNA", "3CA_Couturier", "Data_Couturier2020_Brain", "Genes.txt"))]
    coord = pd.read_csv(os.path.join(T, "C2_LB_gene_coords.tsv"), sep="\t")
    coord = coord[~coord.notfound.fillna(False)].copy()
    coord = coord[coord.chrom.astype(str).str.replace("chr", "", regex=False).isin([str(i) for i in range(1, 23)] + ["X"])].copy()
    coord["chrom"] = coord["chrom"].astype(str).str.replace("chr", "", regex=False)
    coord["start"] = coord["start"].astype(float); coord["end"] = coord["end"].astype(float)
    coord = coord.drop_duplicates("symbol")
    coord = coord.set_index("symbol")
    coord = coord.reindex(genes).dropna(subset=["chrom", "start", "end"])
    gene_keep = coord.index.tolist()
    keep_pos = [genes.index(x) for x in gene_keep]
    if os.path.exists(fcache):
        M = sp.load_npz(fcache)
        print("cached mapped csc", M.shape, "nnz", M.nnz, flush=True)
    else:
        A = sp.load_npz(os.path.join(T, "couturier_matrix_csr.npz"))
        print("full csr", A.shape, "nnz", A.nnz, flush=True)
        M = A.tocsr()[keep_pos, :].tocsc()
        sp.save_npz(fcache, M)
        print("built mapped csc", M.shape, "nnz", M.nnz, flush=True)
    return M, coord, gene_keep

def build_adata(M, coord, cell_idx, obs):
    X = M[:, cell_idx].T.tocsr().astype(np.float32)
    adata = anndata.AnnData(X=X)
    adata.obs = obs.reset_index(drop=True)
    vd = coord.copy()
    vd.index = range(len(vd))
    adata.var = vd[["chrom", "start", "end"]].copy()
    adata.var_names = coord.index.tolist()
    # normalise per cell then log1p (documented for infercnvpy on TPM)
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

if __name__ == "__main__":
    import anndata, scanpy as sc
    import infercnvpy as cnv
    t0 = time.time()
    M, coord, gene_keep = load_full_csc()
    print("load+mapped ok", time.time() - t0, flush=True)
    meta = pd.read_csv(os.path.join(T, "couturier_meta_layers.tsv"), sep="\t")
    meta["is_gsc"] = meta["sample"].str.contains("-GSC|_GSC", case=False, regex=True)
    ref_ok = meta[~meta.mal_annot & ~meta.is_gsc & meta.sample.str.startswith(("HFA", "NSC"))].copy()
    print("reference pool", len(ref_ok), ref_ok.sample.value_counts().to_dict(), flush=True)
    rng = np.random.RandomState(42)
    ref_sel = ref_ok.groupby("sample", group_keys=False).apply(lambda d: d.sample(min(len(d), 150), random_state=42))
    print("ref selected", len(ref_sel), flush=True)
    sample = "BT407"
    cells = meta[meta.sample == sample]
    print(sample, "cells", len(cells), "unannot", (~cells.mal_annot).sum(), flush=True)
    keep_obs = pd.concat([cells, ref_sel], ignore_index=True)
    keep_obs["cnv_group"] = np.where(keep_obs.sample == sample, "tumor", "normal")
    keep_idx = keep_obs.index.values
    adata = build_adata(M, coord, keep_idx, keep_obs)
    print("adata", adata.shape, "groups", adata.obs.cnv_group.value_counts().to_dict(), flush=True)
    t1 = time.time()
    cnv.tl.infercnv(adata, reference_key="cnv_group", reference_cat="normal",
                    window_size=100, step=10, dynamic_threshold=1.5, lfc_clip=3, n_jobs=2)
    print("infercnv done", time.time() - t1, flush=True)
    cnv.tl.cnv_score(adata)
    scores = adata.obs[["cell_name", "sample", "mal_annot", "cnv_group", "cnv_score"]].copy()
    print(scores.groupby("cnv_group")["cnv_score"].describe().to_string(), flush=True)
    scores.to_csv(os.path.join(T, "proto_bt407_cnv_scores.tsv"), sep="\t", index=False)
    print("prototype saved; total", time.time() - t0, flush=True)

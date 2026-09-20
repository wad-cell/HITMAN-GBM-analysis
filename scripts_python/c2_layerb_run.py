# C2.11 Layer B CNV run per sample (infercnvpy). Usage: python c2_layerb_run.py BT407 [prototype]
import os, sys, time, numpy as np, pandas as pd, scipy.sparse as sp
T = r"${PROJECT_ROOT}\temp"
OUT = os.path.join(T, "layerb_scores")
os.makedirs(OUT, exist_ok=True)

def load_mapped():
    fcache = os.path.join(T, "couturier_mapped_csc_hg38.npz")
    genes = [l.strip() for l in open(os.path.join(T, "scRNA", "3CA_Couturier", "Data_Couturier2020_Brain", "Genes.txt"))]
    coord = pd.read_csv(os.path.join(T, "C2_LB_gene_coords.tsv"), sep="\t")
    coord["chrom"] = coord["chrom"].astype(str).str.replace("chr", "", regex=False)
    coord = coord[coord["chrom"].isin([str(i) for i in range(1, 23)] + ["X"])].copy()
    coord = coord.drop_duplicates("symbol")
    coord["chr_rank"] = pd.to_numeric(coord["chrom"].replace({"X": 23}), errors="coerce").astype(int)
    coord = coord.sort_values(["chr_rank", "start"])
    keep = [g for g in coord["symbol"] if g in set(genes)]
    keep_pos = [genes.index(x) for x in keep]
    coord = coord.set_index("symbol").reindex(keep)
    if os.path.exists(fcache):
        M = sp.load_npz(fcache)
    else:
        A = sp.load_npz(os.path.join(T, "couturier_matrix_csr.npz"))
        M = A.tocsr()[keep_pos, :].tocsc()
        sp.save_npz(fcache, M)
    return M, coord

def main(sample, prototype=False):
    import anndata, scanpy as sc
    import infercnvpy as cnv
    t0 = time.time()
    M, coord = load_mapped()
    meta = pd.read_csv(os.path.join(T, "couturier_meta_layers.tsv"), sep="\t")
    meta["is_gsc"] = meta["sample"].str.contains("-GSC|_GSC", case=False, regex=True)
    ref_pool = meta[~meta.mal_annot & ~meta.is_gsc & meta["sample"].str.startswith(("HFA", "NSC"))].copy()
    n_ref = 100 if prototype else 150
    refs = ref_pool.groupby("sample", group_keys=False).apply(lambda d: d.sample(min(len(d), n_ref), random_state=42))
    # C2.11 rule: run CNV inference only on cells WITHOUT official malignant annotation.
    cells = meta[(meta["sample"] == sample) & (~meta["mal_annot"])].copy()
    keep_obs = pd.concat([cells, refs], ignore_index=True)
    keep_obs["cnv_group"] = np.where(keep_obs["sample"] == sample, "tumor", "normal")
    keep_idx = keep_obs.index.values
    X = M[:, keep_idx].T.tocsr().astype(np.float32)
    adata = anndata.AnnData(X=X)
    adata.obs = keep_obs.reset_index(drop=True)
    vd = coord.copy(); vd.index = range(len(vd))
    adata.var = vd[["chrom", "start", "end"]].rename(columns={"chrom": "chromosome"}).copy()
    adata.var["chromosome"] = "chr" + adata.var["chromosome"].astype(str)
    adata.var_names = coord.index.tolist()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("adata ready", adata.shape, time.time() - t0, flush=True)
    t1 = time.time()
    cnv.tl.infercnv(adata, reference_key="cnv_group", reference_cat="normal",
                    window_size=100, step=10, dynamic_threshold=1.5, lfc_clip=3, n_jobs=2)
    print("infercnv done", round(time.time() - t1, 1), flush=True)
    Xc = adata.obsm["X_cnv"]
    adata.obs["cnv_score"] = np.asarray(abs(Xc).mean(axis=1)).ravel()
    sc_ = adata.obs[["cell_name", "sample", "mal_annot", "cnv_group", "cnv_score"]].copy()
    refs_sc = sc_[sc_.cnv_group == "normal"]["cnv_score"]
    thr = float(refs_sc.mean() + 3 * refs_sc.std())
    sc_["lb_call"] = np.where(sc_.cnv_group == "tumor",
                              np.where(sc_["cnv_score"] > thr, "malignant_cnv", "nonmalignant_diploid"), "reference")
    print("ref mean sd thr", round(refs_sc.mean(), 4), round(refs_sc.std(), 4), round(thr, 4), flush=True)
    print(sc_[sc_.cnv_group == "tumor"].lb_call.value_counts().to_dict(), flush=True)
    if prototype:
        print(sc_[sc_.cnv_group == "tumor"].groupby(["sample","mal_annot","lb_call"])["cnv_score"].agg(["count","mean"]).to_string(), flush=True)
    sc_.to_csv(os.path.join(OUT, f"{sample}_cnv.tsv"), sep="\t", index=False)
    print("saved", sample, "total", round(time.time() - t0, 1), flush=True)

if __name__ == "__main__":
    sample = sys.argv[1]
    proto = len(sys.argv) > 2 and sys.argv[2] == "prototype"
    main(sample, prototype=proto)

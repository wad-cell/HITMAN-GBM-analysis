suppressMessages({ library(oligo); library(pd.clariom.s.human); library(clariomshumantranscriptcluster.db); library(Biobase); library(limma) })
options(warn=1)
base <- "${PROJECT_ROOT}"
outdir <- file.path(base, "output"); tmpdir <- file.path(base, "temp")
logf <- file.path(tmpdir, "C1_limma.log"); con <- file(logf, open="wt"); sink(con, type="output"); sink(con, type="message")
cat("start", format(Sys.time()), "\n"); flush(con)
eset <- readRDS(file.path(tmpdir, "C1_eset_rma.rds"))
cat("eset dim:", dim(exprs(eset)), "\n"); flush(con)
# ---- feature annotation: transcript cluster -> symbol ----
pids <- featureNames(eset)
anno <- AnnotationDbi::select(clariomshumantranscriptcluster.db, keys=pids, columns=c("SYMBOL","ENTREZID","GENENAME","GENETYPE"), keytype="PROBEID")
cat("anno rows:", nrow(anno), "unique probes:", length(unique(anno$PROBEID)), "\n"); flush(con)
# first match per probe (avoid duplicated annotation rows)
anno <- anno[!duplicated(anno$PROBEID), ]
rownames(anno) <- anno$PROBEID
fData(eset) <- anno[featureNames(eset), ]
saveRDS(eset, file.path(tmpdir, "C1_eset_rma_annotated.rds"))
cat("annotation saved. mapped symbols:", sum(!is.na(fData(eset)$SYMBOL)), "\n"); flush(con)
# ---- QC ----
expr <- exprs(eset)
pd <- pData(eset)
pd$group <- paste(pd$cell_line, pd$stimulus, sep="_")
cormat <- cor(expr)
png(file.path(tmpdir, "C1_QC_correlation.png"), width=2400, height=2200, res=160)
heatmap(cormat, symm=TRUE, labRow=pd$sample, labCol=pd$sample, col=colorRampPalette(c("navy","white","firebrick3"))(100), main="Sample correlation (RMA)")
dev.off()
pca <- prcomp(t(expr), center=TRUE, scale.=FALSE)
pve <- round(100*summary(pca)$importance[2,], 2)
png(file.path(tmpdir, "C1_QC_PCA.png"), width=2200, height=1800, res=160)
cols <- as.numeric(factor(pd$cell_line)); pchs <- as.numeric(factor(pd$stimulus))
plot(pca$x[,1], pca$x[,2], col=cols, pch=pchs, xlab=paste0("PC1 (",pve[1],"%)"), ylab=paste0("PC2 (",pve[2],"%)"), main="PCA RMA expr")
legend("topright", legend=levels(factor(pd$cell_line)), col=1:3, pch=15, title="cell line")
legend("bottomright", legend=levels(factor(pd$stimulus)), pch=1:3, title="stimulus")
dev.off()
hc <- hclust(dist(t(expr)))
png(file.path(tmpdir, "C1_QC_dendro.png"), width=1800, height=1200, res=150)
plot(hc, labels=pd$sample, main="hclust (Euclidean on RMA)")
dev.off()
# QC summary: corr within/between groups
qc_out <- data.frame(sample=pd$sample, cell_line=pd$cell_line, stimulus=pd$stimulus,
  mean_expr=round(colMeans(expr),3), median_corr_to_same=NA, median_corr_to_diff=NA)
grp <- pd$group
for(i in 1:nrow(qc_out)){
  same <- grp==grp[i]; same[i] <- FALSE
  qc_out$median_corr_to_same[i] <- round(median(cormat[i, same]),3)
  qc_out$median_corr_to_diff[i] <- round(median(cormat[i, !same & seq_len(nrow(qc_out))!=i]),3)
}
write.csv(qc_out, file.path(tmpdir, "C1_QC_summary.csv"), row.names=FALSE)
cat("QC saved\n"); flush(con)
# ---- limma per cell line ----
write_limma <- function(cl, stim, lv_control="ctrl"){
  keep <- pd$cell_line==cl & pd$stimulus %in% c(lv_control, stim)
  e2 <- eset[, keep]; pd2 <- pData(e2)
  pd2$stimulus <- factor(pd2$stimulus, levels=c(lv_control, stim))
  design <- model.matrix(~0 + stimulus, data=pd2)
  colnames(design) <- gsub("stimulus","",colnames(design))
  fit <- lmFit(exprs(e2), design)
  contrast <- makeContrasts(contrasts=paste0(stim, "-", lv_control), levels=design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, coef=1, number=Inf, sort.by="none")
  tt$probe <- rownames(tt)
  tt$symbol <- fData(e2)[rownames(tt), "SYMBOL"]
  tt$entrez <- fData(e2)[rownames(tt), "ENTREZID"]
  tt$rank_stat <- sign(tt$t) * (-log10(tt$P.Value))
  colnames(tt)[colnames(tt)=="logFC"] <- "logFC"
  tt <- tt[, c("probe","symbol","entrez","logFC","AveExpr","t","P.Value","adj.P.Val","B","rank_stat")]
  tt
}
res_list <- list()
for(cl in c("U87-MG","KNS-42","GIN-28")){
  for(stim in c("TTFields","DBS")){
    nm <- paste0(cl, "_", stim)
    tt <- write_limma(cl, stim)
    # collapse by symbol: keep top |t| probe per symbol
    tt <- tt[order(tt$symbol, -abs(tt$t)), ]
    tt <- tt[!duplicated(tt$symbol) & !is.na(tt$symbol), ]
    rownames(tt) <- NULL
    res_list[[nm]] <- tt
    fn <- file.path(outdir, paste0("C1_", gsub("-","",cl), "_", stim, "_vs_CTRL.tsv"))
    write.table(tt, fn, sep="\t", quote=FALSE, row.names=FALSE)
    cat("wrote", basename(fn), nrow(tt), "genes\n"); flush(con)
  }
}
cat("limma done", format(Sys.time()), "\n"); flush(con)
sink(type="message"); sink(type="output"); close(con)

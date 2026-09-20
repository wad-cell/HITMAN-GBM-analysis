suppressPackageStartupMessages({
  library(DESeq2)
})

cat("=== B3.5 R/DESeq2 engine verification ===\n", file=stdout())
cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")

# ---- Load rounded counts (TSV, gene x 20) ----
counts_df <- read.table("temp/B35_counts.tsv", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
gene_ids <- counts_df$gene_id
counts_mat <- as.matrix(counts_df[, -1, drop=FALSE])
rownames(counts_mat) <- gene_ids
stopifnot(all(!is.na(counts_mat)))
cat("Counts dim:", nrow(counts_mat), "x", ncol(counts_mat), "\n")

# ---- Load sample metadata ----
meta <- read.csv("temp/b3_sample_metadata.csv", stringsAsFactors=FALSE)
cat("Metadata samples:", nrow(meta), "\n")

sample_names <- meta$sample
MAN <- factor(meta$MAN, levels=c("ND","MAN"))
Field <- factor(meta$Field, levels=c("NF","F"))
colData <- data.frame(row.names=sample_names, MAN=MAN, Field=Field)

# Ensure count columns align with metadata order
colnames(counts_mat) <- sample_names
counts_mat <- counts_mat[, sample_names, drop=FALSE]

# ---- Build DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=colData, design=~MAN*Field)

cat("\n--- model.matrix() ---\n")
mm <- model.matrix(~MAN*Field, colData)
print(colnames(mm))
write.csv(mm, "temp/B35_model_matrix.csv")

cat("\n--- resultsNames after DESeq ---\n")

# ---- Run DESeq2 ----
cat("Running DESeq...\n")
dds <- DESeq(dds, quiet=FALSE)
cat("DESeq complete.\n")

rns <- resultsNames(dds)
print(rns)
writeLines(rns, "temp/B35_results_names.txt")

# ---- Extract C1-C5 using same coefficient vectors as PyDESeq2 ----
# Design columns (R order): Intercept, MAN_MAN_vs_ND, Field_F_vs_NF, MAN_MAN.Field_F
# PyDESeq2 used: C1=[0,0,1,0], C2=[0,1,0,0], C3=[0,0,1,1], C4=[0,1,1,1], C5=[0,0,0,1]
vecs <- list(
  C1 = c(0,0,1,0),
  C2 = c(0,1,0,0),
  C3 = c(0,0,1,1),
  C4 = c(0,1,1,1),
  C5 = c(0,0,0,1)
)

save_RDS <- function(obj, path) saveRDS(obj, path)

for (cn in names(vecs)) {
  cat("Extracting", cn, "...\n")
  res <- results(dds, contrast=vecs[[cn]], alpha=0.05)
  save_RDS(res, file.path("temp", paste0("B35_R_", cn, ".rds")))
  df <- as.data.frame(res)
  df <- data.frame(gene=rownames(df), baseMean=df$baseMean, log2FC=df$log2FoldChange,
                   lfcSE=df$lfcSE, stat=df$stat, pvalue=df$pvalue, padj=df$padj,
                   stringsAsFactors=FALSE)
  write.table(df, file.path("temp", paste0("B35_R_", cn, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
  n_sig <- sum(df$padj < 0.05, na.rm=TRUE)
  n_eval <- sum(!is.na(df$padj))
  cat(cn, "padj-evaluable:", n_eval, "FDR<0.05:", n_sig, "\n")
}

# ---- Independent filtering audit for C1-C5 ----
cat("\n--- Independent filtering audit ---\n")
audit_rows <- list()
for (cn in names(vecs)) {
  res <- readRDS(file.path("temp", paste0("B35_R_", cn, ".rds")))
  input_universe <- nrow(res)
  stat_available <- sum(!is.na(res$stat))
  # genes that have pvalue (tested)
  pvalue_available <- sum(!is.na(res$pvalue))
  # after independent filtering / p-adj computation
  padj_evaluable <- sum(!is.na(res$padj))
  padj_na <- sum(is.na(res$padj) & !is.na(res$pvalue))
  n_removed_by_filtering <- pvalue_available - padj_evaluable
  fdr_sig <- sum(res$padj < 0.05, na.rm=TRUE)
  audit_rows[[cn]] <- data.frame(
    contrast=cn,
    input_universe=input_universe,
    statistic_available=stat_available,
    pvalue_available=pvalue_available,
    padj_evaluable=padj_evaluable,
    padj_NA=padj_na,
    removed_by_independent_filtering=n_removed_by_filtering,
    FDR_005=fdr_sig,
    stringsAsFactors=FALSE
  )
  cat(cn, sprintf("input=%d stat=%d pval=%d padj_eval=%d padjNA=%d removed=%d FDR=%d\n",
      input_universe, stat_available, pvalue_available, padj_evaluable, padj_na,
      n_removed_by_filtering, fdr_sig))
}
audit <- do.call(rbind, audit_rows)
write.table(audit, "temp/B35_filtering_audit_r.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# ---- Save session info ----
sink("temp/B35_R_SESSIONINFO.txt")
cat("R version:\n"); print(R.version)
cat("\nsessionInfo():\n"); print(sessionInfo())
cat("\nDESeq2 version:\n"); print(packageVersion("DESeq2"))
cat("\nmodel.matrix():\n"); print(mm)
cat("\nresultsNames():\n"); print(rns)
sink()

cat("\n=== B3.5 R engine run complete ===\n")

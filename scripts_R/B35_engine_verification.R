# B35_ENGINE_VERIFICATION.R
# R/DESeq2 engine verification for GSE324656 B3
# design = ~ MAN * Field; rounded Salmon counts; 5 contrasts; BH FDR < 0.05
# Compares with pydeseq2 0.5.4 results saved in temp/b3_primary_artifacts.pkl

options(warn = 2)
options(stringsAsFactors = FALSE)
options(java.opts = "-Xmx6g")

R.version <- paste0(R.version$major, ".", R.version$minor)

library(DESeq2)
library(matrixStats)
library(gdata)

# --- Load rounded matrix ---
rmat <- readRDS("temp/b3_rounded_counts_rds.rds")
# rmat should be 62757 x 20 with gene symbols as rownames
if (is.null(rownames(rmat)) || any(nchar(rownames(rmat)) == 0)) {
  stop("Rounded matrix rownames are empty. Check b3_rounded_counts_rds.rds")
}

# --- Load sample metadata ---
meta <- readRDS("temp/b3_sample_metadata.rds")
# meta should have: Sample, group, Design (MAN=ND/MAN, Field=NF/F)

# --- Design matrix ---
dds <- DESeqDataSetFromMatrix(countData = rmat,
                               colData = meta,
                               design = ~ MAN * Field)
# Set references
dds$MAN <- relevel(dds$MAN, ref = "ND")
dds$Field <- relevel(dds$Field, ref = "NF")

mm <- model.matrix(dds)
cat("model.matrix() dimensions:", dim(mm), "\n")
cat("rank:", qr(mm)$rank, "\n")
cat("Condition number:", round(cod0eR:::condnum(mm), 2), "\n")
write.table(mm, "temp/B35_model_matrix.tsv", sep="\t", quote=FALSE)

# --- Save model matrix ---
cat("Model matrix columns:", paste(colnames(mm), collapse=", "), "\n")

# --- Save design and resultsNames ---
cat("Design:", deparse(dds@design), "\n")
cat("levels MAN:", levels(dds$MAN), "\n")
cat("levels Field:", levels(dds$Field), "\n")
cat("resultsNames(dds):", paste(resultsNames(dds), collapse=", "), "\n")
write.table(data.frame(sample=meta$Sample, group=meta$group, design=meta$Design, stringsAsFactors=FALSE),
            "temp/B35_sample_design.tsv", sep="\t", quote=FALSE)

# --- Run DESeq ---
dds <- DESeq(dds)

# --- Contrast definitions ---
contrast_map <- list(
  C1 = list(name="ND-F - ND-NF", vec=NULL, LHS=c(1), RHS=NULL, interaction=FALSE),
  C2 = list(name="MAN-NF - ND-NF", vec=NULL, LHS=c(2), RHS=NULL, interaction=FALSE),
  C3 = list(name="MAN-F - MAN-NF", vec=NULL, LHS=c(2,3), RHS=c(3), interaction=FALSE),
  C4 = list(name="MAN-F - ND-NF", vec=NULL, LHS=c(1,2,3), RHS=c(1,3), interaction=FALSE),
  C5 = list(name="(MAN-F - MAN-NF) - (ND-F - ND-NF) = MAN*Field",
            vec=c(0,1,-1,1), LHS=NULL, RHS=NULL, interaction=TRUE)
)

# --- Get results ---
all_results <- list()
all_pvals <- list()

for (cm in names(contrast_map)) {
  res <- results(dds, contrast=contrast_map[[cm]], alpha=0.05)
  all_results[[cm]] <- res
  pvals <- as.numeric(res$pvalue)
  pvals <- pvals[!is.na(pvals)]
  all_pvals[[cm]] <- pvals
}

# --- Save raw results as RDS for Python ---
for (cm in names(contrast_map)) {
  write.csv(as.data.frame(all_results[[cm]]),
            file.path("temp", paste0("B35_DESeq2_", cm, ".csv")),
            row.names=TRUE)
}

# --- Python comparison ---
cat("\n=== COMPARING WITH PyDESeq2 ===\n")

# Load pydeseq2 results
py_c4 <- readRDS("temp/b3_primary_artifacts.pkl")$results$C4
py_c5 <- readRDS("temp/b3_primary_artifacts.pkl")$results$C5
py_c1 <- readRDS("temp/b3_primary_artifacts.pkl")$results$C1
py_c2 <- readRDS("temp/b3_primary_artifacts.pkl")$results$C2
py_c3 <- readRDS("temp/b3_primary_artifacts.pkl")$results$C3

py_results <- list(C1=py_c1, C2=py_c2, C3=py_c3, C4=py_c4, C5=py_c5)

# Common genes (ENSG unversioned)
gene_names <- as.character(rownames(all_results$C4))
# PyDESeq2 gene names
py_gene_names <- as.character(py_c4$gene)
# Both should be ENSG unversioned

# For each contrast, compare
comp <- list()
for (cm in names(contrast_map)) {
  res <- all_results[[cm]]
  pyres <- py_results[[cm]]
  
  # Only use genes with padj in both
  r_avail <- !is.na(res$padj)
  py_avail <- !is.na(pyres$padj)
  
  # Get common genes with padj in both
  r_idx <- which(r_avail)
  py_idx <- which(py_avail)
  
  # Match by gene name
  r_genes <- gene_names[r_idx]
  py_genes <- as.character(pyres$gene)
  
  common <- intersect(r_genes, py_genes)
  if (length(common) == 0) {
    cat(cm, "no common genes with padj\n")
    next
  }
  
  # Get py results for common genes (in order)
  py_match_idx <- match(common, py_genes)
  py_match <- pyres[py_match_idx, , drop=FALSE]
  
  r_match <- res[common, , drop=FALSE]
  
  # Pearson/Spearman log2FC
  r_log2fc <- as.numeric(r_match$log2FoldChange)
  py_log2fc <- as.numeric(py_match$log2FoldChange)
  
  r_stat <- as.numeric(r_match$stat)
  py_stat <- as.numeric(py_match$stat)
  
  r_log2fc <- r_log2fc[is.finite(r_log2fc) & is.finite(py_log2fc)]
  py_log2fc <- py_log2fc[is.finite(r_log2fc) & is.finite(py_log2fc)]
  r_stat <- r_stat[is.finite(r_log2fc) & is.finite(py_log2fc)]
  py_stat <- py_stat[is.finite(r_log2fc) & is.finite(py_log2fc)]
  
  # Sign concordance
  r_sign <- sign(r_log2fc)
  py_sign <- sign(py_log2fc)
  sign_conc <- sum(r_sign == py_sign & r_sign != 0) / sum(r_sign != 0)
  
  # FDR<0.05 counts
  r_sig <- as.numeric(r_match$padj) < 0.05 & is.finite(as.numeric(r_match$padj))
  py_sig <- as.numeric(py_match$padj) < 0.05 & is.finite(as.numeric(py_match$padj))
  r_sig_count <- sum(r_sig)
  py_sig_count <- sum(py_sig)
  
  # Jaccard
  common_all <- intersect(which(r_sig), which(py_sig))
  union_all <- union(which(r_sig), which(py_sig))
  jaccard <- length(common_all) / length(union_all)
  
  # Top overlap (by |stat|)
  r_stat_abs <- abs(r_stat)
  py_stat_abs <- abs(py_stat)
  r_order <- order(-r_stat_abs)
  py_order <- order(-py_stat_abs)
  
  r_top <- common[which(r_sig)[r_order[1:min(250, sum(r_sig))]]]
  py_top <- common[which(py_sig)[py_order[1:min(250, sum(py_sig))]]]
  
  top250_inter <- length(intersect(r_top, py_top))
  top250_union <- length(union(r_top, py_top))
  top250_jac <- top250_inter / top250_union
  
  top50_inter <- length(intersect(common[which(r_sig)[r_order[1:50]]], common[which(py_sig)[py_order[1:50]]]))
  top50_union <- length(union(common[which(r_sig)[r_order[1:50]]], common[which(py_sig)[py_order[1:50]]]))
  top50_jac <- top50_inter / top50_union
  
  top100_inter <- length(intersect(common[which(r_sig)[r_order[1:100]]], common[which(py_sig)[py_order[1:100]]]))
  top100_union <- length(union(common[which(r_sig)[r_order[1:100]]], common[which(py_sig)[py_order[1:100]]]))
  top100_jac <- top100_inter / top100_union
  
  comp[[cm]] <- list(
    pearson_log2fc = round(cor(r_log2fc, py_log2fc, method="pearson"), 4),
    spearman_log2fc = round(cor(r_log2fc, py_log2fc, method="spearman"), 4),
    pearson_stat = round(cor(r_stat, py_stat, method="pearson"), 4),
    spearman_stat = round(cor(r_stat, py_stat, method="spearman"), 4),
    sign_concordance = round(sign_conc, 4),
    r_sig = r_sig_count,
    py_sig = py_sig_count,
    jaccard = round(jaccard, 4),
    top250_jac = round(top250_jac, 4),
    top100_jac = round(top100_jac, 4),
    top50_jac = round(top50_jac, 4),
    n_common = length(common),
    r_padj_na = sum(is.na(as.numeric(r_match$padj))),
    py_padj_na = sum(is.na(as.numeric(py_match$padj))),
    r_padj_evaluable = sum(is.finite(as.numeric(r_match$padj))),
    py_padj_evaluable = sum(is.finite(as.numeric(py_match$padj))),
    n_input = length(common)
  )
}

# --- C5-specific: check how many of primary C5 1891 replicate in DESeq2 ---
cat("\n=== C5 SPECIFIC: PRIMARY vs DESeq2 ===\n")
c5_common <- intersect(as.character(py_c5$gene), gene_names)
c5_py_idx <- which(py_c5$gene %in% c5_common)
c5_py_sig <- as.numeric(py_c5$padj[c5_py_idx]) < 0.05 & is.finite(as.numeric(py_c5$padj[c5_py_idx]))
cat("Primary C5 sig:", sum(c5_py_sig), "/ 1891\n")
c5_common_genes <- c5_common[c5_py_sig]

# How many of those 1891 replicate in DESeq2 with same direction?
c5_r_match <- all_results$C5[c5_common_genes, , drop=FALSE]
c5_r_padj <- as.numeric(c5_r_match$padj)
c5_r_sig <- c5_r_padj < 0.05 & is.finite(c5_r_padj)
replicated_sig <- sum(c5_r_sig)
cat("DESeq2 C5 sig among primary 1891:", replicated_sig, "/", length(c5_common_genes), "\n")

# Same direction?
c5_r_log2fc <- as.numeric(c5_r_match$log2FoldChange)
c5_py_log2fc <- as.numeric(py_c5$log2FoldChange[c5_py_idx])
c5_log2fc_both <- c5_r_log2fc * c5_py_log2fc
same_dir <- sum(c5_log2fc_both > 0 & c5_r_sig) / sum(c5_r_sig)
cat("Same direction among replicated:", round(same_dir, 4), "\n")

# --- Stable 1602 genes check ---
cat("\n=== STABLE 1602 GENE REPLICATION ===\n")
stable <- read.csv("temp/b3_stable_genes_c5.tsv", sep="\t", stringsAsFactors=FALSE)
stable_ensg <- stable$gene_id
# Convert ENSG versioned to unversioned
stable_unversioned <- gsub("\\.v[0-9]+$", "", stable_ensg)

# Match with DESeq2 gene names
stable_in_r <- stable_unversioned[stable_unversioned %in% gene_names]
cat("Stable 1602 genes present in DESeq2 universe:", length(stable_in_r), "\n")

# How many of those replicate in DESeq2 C5?
c5_stable_r <- all_results$C5[stable_in_r, , drop=FALSE]
c5_stable_r_padj <- as.numeric(c5_stable_r$padj)
c5_stable_r_sig <- c5_stable_r_padj < 0.05 & is.finite(c5_stable_r_padj)
c5_stable_same_dir <- c5_stable_r_log2fc <- as.numeric(c5_stable_r$log2FoldChange)
c5_stable_py_log2fc <- as.numeric(py_c5$log2FoldChange[match(stable_in_r, as.character(py_c5$gene))])
c5_stable_log2fc_both <- c5_stable_r_log2fc * c5_stable_py_log2fc
stable_replicated <- sum(c5_stable_r_sig)
stable_same_dir <- sum(c5_stable_log2fc_both > 0 & c5_stable_r_sig) / sum(c5_stable_r_sig)
cat("Stable 1602 replicated in DESeq2 C5:", stable_replicated, "/", length(stable_in_r), "(", round(100*stable_replicated/length(stable_in_r)), "%)\n")
cat("Stable same direction:", round(100*stable_same_dir), "%\n")

# --- Filtering audit ---
cat("\n=== FILTERING AUDIT ===\n")
filter_audit <- data.frame(
  contrast = character(),
  universe = integer(),
  stats_available = integer(),
  padj_evaluable = integer(),
  padj_na = integer(),
  fdr_sig = integer(),
  stringsAsFactors = FALSE
)

for (cm in names(contrast_map)) {
  res <- all_results[[cm]]
  universe <- nrow(res)
  stats_avail <- sum(!is.na(res$stat))
  padj_finite <- sum(is.finite(as.numeric(res$padj)))
  padj_na <- sum(is.na(as.numeric(res$padj)))
  sig <- sum(as.numeric(res$padj) < 0.05, na.rm=TRUE)
  filter_audit <- rbind(filter_audit, data.frame(
    contrast=cm, universe=universe, stats_available=stats_avail,
    padj_evaluable=padj_finite, padj_na=padj_na, fdr_sig=sig,
    stringsAsFactors=FALSE
  ))
}
write.csv(filter_audit, "temp/B35_filtering_audit.csv", row.names=FALSE)

# --- Save results for Python summary ---
comp_out <- list()
for (cm in names(comp)) comp_out[[cm]] <- comp[[cm]]
saveRDS(comp_out, "temp/B35_engine_comparison.rds")

# --- R sessionInfo ---
writeLines(capture.output(sessionInfo()), "temp/B35_sessioninfo.txt")

# --- R version ---
cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", packageVersion("DESeq2"), "\n")

cat("\n=== DONE ===\n")

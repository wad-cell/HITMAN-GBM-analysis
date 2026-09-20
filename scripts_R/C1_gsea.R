suppressMessages({ library(fgsea) })
options(warn=1)
base <- "${PROJECT_ROOT}"
outdir <- file.path(base, "output"); tmpdir <- file.path(base, "temp")
logf <- file.path(tmpdir, "C1_gsea.log"); con <- file(logf, open="wt"); sink(con, type="output"); sink(con, type="message")
read_gmt <- function(fn){
  lines <- readLines(fn); lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(ln){ sp <- strsplit(ln, "\t")[[1]]; list(name=sp[1], genes=unique(sp[-(1:2)][nzchar(sp[-(1:2)])])) })
  nm <- sapply(sets, `[[`, "name"); gl <- lapply(sets, `[[`, "genes"); names(gl) <- nm; gl
}
hall <- read_gmt(file.path(tmpdir,"b3_gslib_HALLMARK.gmt"))
react <- read_gmt(file.path(tmpdir,"b3_gslib_REACTOME.gmt"))
gobp <- read_gmt(file.path(tmpdir,"b3_gslib_GOBP.gmt"))
cat("gmt sizes", length(hall), length(react), length(gobp), "\n"); flush(con)
# ranked lists from GSE (R/DESeq2 frozen, collapsed representative per symbol)
mkrank_gse <- function(fn){
  d <- read.delim(fn, stringsAsFactors=FALSE)
  d <- d[!is.na(d$stat) & !is.na(d$log2FC), ]
  r <- d$stat; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
# ranked lists from limma comparisons
mkrank_lim <- function(cl, stim){
  fn <- file.path(outdir, paste0("C1_", cl, "_", stim, "_vs_CTRL.tsv"))
  d <- read.delim(fn, stringsAsFactors=FALSE)
  d <- d[!is.na(d$t), ]
  r <- d$t; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
lists <- list()
for(cl in c("U87MG","KNS42","GIN28")) for(st in c("TTFields","DBS")) lists[[paste0("EMTAB_",cl,"_",st)]] <- mkrank_lim(cl, st)
lists[["GSE_C4"]] <- mkrank_gse(file.path(tmpdir,"C1_gse_C4_symbol_rep.tsv"))
lists[["GSE_C5"]] <- mkrank_gse(file.path(tmpdir,"C1_gse_C5_symbol_rep.tsv"))
cat("lists:", names(lists), "\n"); flush(con)
for(nm in names(lists)) cat(nm, length(lists[[nm]]), "\n")
run_fgsea <- function(rl, sets, tag, cmp){
  set.seed(42)
  res <- fgsea(pathways=sets, stats=rl, minSize=5, maxSize=500, nproc=4, eps=0)
  res$tag <- tag
  res$comparison <- cmp
  res
}
out_tables <- list()
for(i in seq_along(lists)){
  nm <- names(lists)[i]; rl <- lists[[i]]
  cat("running", nm, format(Sys.time()), "\n"); flush(con)
  h <- run_fgsea(rl, hall, "HALLMARK", nm); h$pathway <- h$pathway
  r <- run_fgsea(rl, react, "REACTOME", nm)
  out_tables[[paste0(nm,"_HALLMARK")]] <- h[, .(pathway, pval, padj, ES, NES, size, leadingEdge)]
  out_tables[[paste0(nm,"_REACTOME")]] <- r[, .(pathway, pval, padj, ES, NES, size, leadingEdge)]
}
comb <- do.call(rbind, lapply(names(out_tables), function(nm){
  d <- out_tables[[nm]]; d$comparison <- sub("_(HALLMARK|REACTOME)$", "", nm); d$db <- ifelse(grepl("HALLMARK", nm), "HALLMARK", "REACTOME")
  d[, .(comparison, db, pathway, pval, padj, ES, NES, size, leadingEdge)]
}))
comb$leadingEdge <- sapply(comb$leadingEdge, function(x) paste(x, collapse="|"))
saveRDS(comb, file.path(tmpdir,"C1_gsea_all.rds"))
write.table(comb, file.path(tmpdir,"C1_gsea_all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("gsea done", format(Sys.time()), "\n"); flush(con)
sink(type="message"); sink(type="output"); close(con)

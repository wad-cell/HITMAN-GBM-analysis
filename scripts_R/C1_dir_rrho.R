suppressMessages({ library(fgsea); library(RRHO) })
options(warn=1)
base <- "${PROJECT_ROOT}"
outdir <- file.path(base, "output"); tmpdir <- file.path(base, "temp")
logf <- file.path(tmpdir, "C1_dir_rrho.log"); con <- file(logf, open="wt"); sink(con, type="output"); sink(con, type="message")
cat("start", format(Sys.time()), "\n"); flush(con)
read_gmt <- function(fn){
  lines <- readLines(fn); lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(ln){ sp <- strsplit(ln, "\t")[[1]]; list(name=sp[1], genes=unique(sp[-(1:2)][nzchar(sp[-(1:2)])])) })
  gl <- lapply(sets, `[[`, "genes"); names(gl) <- sapply(sets, `[[`, "name"); gl
}
dirsets <- read_gmt(file.path(tmpdir,"C1_directional_sets.gmt"))
mkrank_lim <- function(cl, stim){
  fn <- file.path(outdir, paste0("C1_", cl, "_", stim, "_vs_CTRL.tsv"))
  d <- read.delim(fn, stringsAsFactors=FALSE)
  d <- d[!is.na(d$t), ]
  r <- d$t; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
lists <- list()
for(cl in c("U87MG","KNS42","GIN28")) for(st in c("TTFields","DBS")) lists[[paste0("EMTAB_",cl,"_",st)]] <- mkrank_lim(cl, st)
set.seed(1)
res_dir <- list()
for(nm in names(lists)){
  rl <- lists[[nm]]
  fg <- fgsea(pathways=dirsets, stats=rl, minSize=5, maxSize=2000, nproc=4, eps=0)
  res_dir[[nm]] <- fg[, .(pathway, pval, padj, ES, NES, size)]
  cat(nm, "dirsets done\n"); flush(con)
}
dir_all <- do.call(rbind, lapply(names(res_dir), function(nm){ d <- res_dir[[nm]]; d$comparison <- nm; d[, .(comparison, pathway, pval, padj, ES, NES, size)] }))
write.table(dir_all, file.path(tmpdir,"C1_directional_gsea.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# RRHO: compare E-MTAB ranked lists vs GSE C4/C5 ranked by signed t/stat (abs-rank approach: use signed two-sided p-like ranking = order by stat desc gives top-up first; for RRHO we need list sorted by differential evidence; use -log10(p)*sign => p-like)
mkrank_lim_p <- function(cl, stim){
  fn <- file.path(outdir, paste0("C1_", cl, "_", stim, "_vs_CTRL.tsv"))
  d <- read.delim(fn, stringsAsFactors=FALSE); d <- d[!is.na(d$t), ]
  d$score <- sign(d$t) * -log10(d$P.Value)
  r <- d$score; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
mkrank_gse_p <- function(fn){
  d <- read.delim(fn, stringsAsFactors=FALSE); d <- d[!is.na(d$stat), ]
  d$score <- sign(d$stat) * -log10(d$pvalue)
  r <- d$score; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
# RRHO requires two vectors of ranked gene ids (sorted). We rank by absolute signed evidence ascending p? RRHO default: list1/list2 = gene IDs ranked by significance of up-regulation? Standard: order genes by p-value with direction ignored not matter; RRHO tests overlap of top-k across thresholds.
rrho_summary <- list()
for(gs in c("C4","C5")){
  gse_rank <- names(mkrank_gse_p(file.path(tmpdir, paste0("C1_gse_", gs, "_symbol_rep.tsv"))))
  for(cl in c("U87MG","KNS42","GIN28")) for(st in c("TTFields","DBS")){
    nm <- paste0("EMTAB_",cl,"_",st)
    lim_rank <- names(mkrank_lim_p(cl, st))
    common <- intersect(gse_rank, lim_rank)
    # take common-only order
    g1 <- gse_rank[gse_rank %in% common]; g2 <- lim_rank[lim_rank %in% common]
    # RRHO wants vectors with names? pass as is (order from most significant)
    rr <- tryCatch(RRHO(g1, g2, stepsize=100, labels=c(paste0("GSE_",gs), nm), alternative="enrichment"), error=function(e){cat("RRHO err",conditionMessage(e),"\n"); NULL})
    if(!is.null(rr)){
      # summary: Pearson correlation of the log p matrix, max overlap zone info not directly exported easily
      cat(gs, cl, st, "RRHO done\n"); flush(con)
    }
  }
}
cat("end", format(Sys.time()), "\n"); flush(con)
sink(type="message"); sink(type="output"); close(con)

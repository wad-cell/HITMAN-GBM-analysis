suppressMessages(library(fgsea))
options(warn=1)
base <- "${PROJECT_ROOT}"
outdir <- file.path(base, "output"); tmpdir <- file.path(base, "temp")
logf <- file.path(tmpdir, "C1_dir.log"); con <- file(logf, open="wt"); sink(con, type="output"); sink(con, type="message")
cat("start", format(Sys.time()), "\n"); flush(con)
read_gmt <- function(fn){
  lines <- readLines(fn); lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(ln){ sp <- strsplit(ln, "\t")[[1]]; list(name=sp[1], genes=unique(sp[-(1:2)][nzchar(sp[-(1:2)])])) })
  gl <- lapply(sets, `[[`, "genes"); names(gl) <- sapply(sets, `[[`, "name"); gl
}
dirsets <- read_gmt(file.path(tmpdir,"C1_directional_sets.gmt"))
mkrank_lim <- function(cl, stim){
  fn <- file.path(outdir, paste0("C1_", cl, "_", stim, "_vs_CTRL.tsv"))
  d <- read.delim(fn, stringsAsFactors=FALSE); d <- d[!is.na(d$t), ]
  r <- d$t; names(r) <- d$symbol; r <- r[!duplicated(names(r))]; sort(r, decreasing=TRUE)
}
lists <- list()
for(cl in c("U87MG","KNS42","GIN28")) for(st in c("TTFields","DBS")) lists[[paste0("EMTAB_",cl,"_",st)]] <- mkrank_lim(cl, st)
set.seed(1)
res_dir <- list()
for(nm in names(lists)){
  fg <- fgsea(pathways=dirsets, stats=lists[[nm]], minSize=5, maxSize=2000, nproc=4, eps=0)
  res_dir[[nm]] <- fg[, .(pathway, pval, padj, ES, NES, size, leadingEdge)]
  cat(nm, "dirsets done\n"); flush(con)
}
dir_all <- do.call(rbind, lapply(names(res_dir), function(nm){ d <- res_dir[[nm]]; d$comparison <- nm; d[, .(comparison, pathway, pval, padj, ES, NES, size, leadingEdge)] }))
dir_all$leadingEdge <- sapply(dir_all$leadingEdge, function(x) paste(x, collapse="|"))
write.table(dir_all, file.path(tmpdir,"C1_directional_gsea.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("dir done", format(Sys.time()), "\n"); flush(con)
sink(type="message"); sink(type="output"); close(con)

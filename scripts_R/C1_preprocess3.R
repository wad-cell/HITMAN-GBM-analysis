suppressMessages({ library(oligo); library(pd.clariom.s.human); library(Biobase) })
options(warn=1)
base <- "${PROJECT_ROOT}"
logf <- file.path(base, "temp", "C1_preprocess3.log"); con <- file(logf, open="wt"); sink(con, type="output"); sink(con, type="message")
cat("start", format(Sys.time()), "\n"); flush(con)
rawdir <- file.path(base, "temp", "EMTAB9043_raw")
smap <- read.csv(file.path(base, "temp", "C1_emt9043_sample_map.csv"), stringsAsFactors=FALSE)
cel_files <- file.path(rawdir, smap$cel_file)
pd <- data.frame(sample=basename(cel_files), cell_line=smap$cell_line, stimulus=smap$stimulus, stringsAsFactors=FALSE)
rownames(pd) <- pd$sample
half <- ceiling(nrow(pd)/2)
read_batch <- function(idx){
  pdx <- pd[idx,,drop=FALSE]
  e <- read.celfiles(cel_files[idx], phenoData=new("AnnotatedDataFrame", data=pdx))
  # align sample order between assayData and phenoData
  pData(e) <- pData(e)[colnames(exprs(e)), , drop=FALSE]
  e
}
es1 <- read_batch(1:half); cat("b1", dim(exprs(es1)), "aligned:", identical(colnames(exprs(es1)), rownames(pData(es1))), "\n"); flush(con)
es2 <- read_batch((half+1):nrow(pd)); cat("b2", dim(exprs(es2)), "aligned:", identical(colnames(exprs(es2)), rownames(pData(es2))), "\n"); flush(con)
es <- combine(es1, es2)
cat("combined", dim(exprs(es)), "\n"); flush(con)
saveRDS(es, file.path(base, "temp", "C1_raw_combined.rds"))
cat("saved_raw\n"); flush(con)
cat("RMA start", format(Sys.time()), "\n"); flush(con)
eset <- rma(es)
cat("RMA done", dim(exprs(eset)), format(Sys.time()), "\n"); flush(con)
saveRDS(eset, file.path(base, "temp", "C1_eset_rma.rds"))
write.csv(exprs(eset), file.path(base, "temp", "C1_expression_rma.csv"), quote=FALSE)
cat("ALL DONE\n"); flush(con)
sink(type="message"); sink(type="output"); close(con)

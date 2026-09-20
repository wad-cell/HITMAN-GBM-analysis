cat("R.version.string:", R.version.string, "\n")
pkgs <- c("oligo","limma","fgsea","msigdbr","clusterProfiler","AnnotationDbi","Biobase","oligoClasses","pd.clariomshuman")
for (p in pkgs) { ok <- suppressWarnings(requireNamespace(p, quietly=TRUE)); if (ok) cat(p, as.character(packageVersion(p)), "\n") else cat(p, "NA\n") }

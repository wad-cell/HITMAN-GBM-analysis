cat("R", R.version.string, "\n")
pk <- c("oligo","limma","clariomshumantranscriptcluster.db","Biobase","fgsea","data.table")
for(p in pk) cat(p, as.character(packageVersion(p)), "\n")

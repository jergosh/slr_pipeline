library(plyr)
setwd("~/Documents/projects/slr_pipeline/")

pdb_cluster_master <- read.table("data/pdb_clustering_clean.tab")

pval_products <- ddply(pdb_cluster_master, .variables=c("V1"), .fun=function(x) { min(1, prod(x$V4))})
hist(pval_products$V1, breaks=100)
hist(pdb_cluster_master$V4, breaks=10000, xlim=c(0, 1))

hist(subset(pdb_cluster_master, V2 > 1)$V4, breaks=1000, xlim=c(0, 1))

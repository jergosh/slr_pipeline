setwd("~/Documents/projects/slr_pipeline")

pvals_4 <- read.table("data/cluster/cluster_4_pvals.tab", header=F)
colnames(pvals_4) <- c("stable_id", "pdb_id", "len", "n_sign", "p.combined", "p.max", "p.number")
pvals_6 <- read.table("data/cluster/cluster_6_pvals.tab", header=F)
colnames(pvals_6) <- c("stable_id", "pdb_id", "len", "n_sign", "p.combined", "p.max", "p.number")
pvals_8 <- read.table("data/cluster/cluster_8_pvals.tab", header=F)
colnames(pvals_8) <- c("stable_id", "pdb_id", "len", "n_sign", "p.combined", "p.max", "p.number")

pvals_4_0.75 <- read.table("data/cluster/cluster_4_0.75_pvals.tab", header=F)
colnames(pvals_4_0.75) <- c("stable_id", "pdb_id", "len", "n_sign", "p.combined", "p.max", "p.number")
pvals_4_0.9 <- read.table("data/cluster/cluster_4_0.9_pvals.tab", header=F)
colnames(pvals_4_0.9) <- c("stable_id", "pdb_id", "len", "n_sign", "p.combined", "p.max", "p.number")


pdf("cluster_pvals.pdf", height=8, width=11)
par(mfrow=c(2, 3))
hist(pvals_4$p.max, breaks=100, main=expression(paste("Clustering p-values, largest cluster (4", ring(A), ")")), xlab="p-value")
hist(pvals_6$p.max, breaks=100, main=expression(paste("Clustering p-values, largest cluster (6", ring(A), ")")), xlab="p-value")
hist(pvals_8$p.max, breaks=100, main=expression(paste("Clustering p-values, largest cluster (8", ring(A), ")")), xlab="p-value")

hist(pvals_4$p.number, breaks=100, main=expression(paste("Clustering p-values, number of clusters (4", ring(A), ")")), xlab="p-value")
hist(pvals_6$p.number, breaks=100, main=expression(paste("Clustering p-values, number of clusters (6", ring(A), ")")), xlab="p-value")
hist(pvals_8$p.number, breaks=100, main=expression(paste("Clustering p-values, number of clusters (8", ring(A), ")")), xlab="p-value")
par(mfrow=c(1,1))
dev.off()

hist(subset(pvals_4, n_sign > 10 & len > 100)$p.number,
     breaks=100, main=expression(paste("Clustering p-values, number of clusters (", omega, "=1)")), xlab="p-value")
hist(subset(pvals_4_0.9, n_sign > 10 & len > 100)$p.number,
     breaks=100, main=expression(paste("Clustering p-values, number of clusters (", omega, "=0.9)")), xlab="p-value")
hist(subset(pvals_4_0.75, n_sign > 10 & len > 100)$p.number,
     breaks=100, main=expression(paste("Clustering p-values, number of clusters (", omega, "=0.75)")), xlab="p-value")

subset(pvals_4, p.max == 1 & n_sign >= 5)

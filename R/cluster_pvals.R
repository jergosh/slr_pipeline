setwd("~/Documents/projects/slr_pipeline/")
pvals <- read.table("pvals.txt", header=F)
colnames(pvals) <- c("stable_id", "pdb_id", "len", "n_sign", "pval_combined", "pval_max", "pval_number")

pvals <- subset(pvals, n_sign >= 5 & len > 100)

pdf("clustering_pvals.pdf", width=11, height=9)
hist(pvals$pval_number, breaks=100, main="Clustering of sites under positive selection",
     xlab="P-value", cex.main=1.5, cex.lab=1.5)
dev.off()

hist(pvals$pval_max, breaks=100)
# hist(pvals$pval_combined)
plot(pvals$pval_number, pvals$pval_max)
cor.test(pvals$pval_number, pvals$pval_max)

nrow(pvals)
sum(p.adjust(pvals$pval_max, method="BH") < 0.05)

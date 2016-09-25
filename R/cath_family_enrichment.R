133 1.20.1070.10 539
106  1.10.510.10 213
675  3.40.50.300 190
347   2.60.40.10 170
257  2.130.10.10 164
724   3.80.10.10 133
222   1.25.10.10 127
50   1.10.238.10 117
154 1.20.1250.20 115
559  3.30.70.330 112
667 3.40.50.1820  92
225   1.25.40.10  89
663  3.40.50.150  68
230   1.25.40.20  66
296   2.30.30.40  64
521   3.30.40.10  62
360  2.60.40.150  60
680  3.40.50.720  58
397  3.10.100.10  56
754  3.90.190.10  56
44   1.10.20.10 53
116 1.10.630.10 53
242  1.50.40.10 51
39  1.10.150.50 44
302  2.30.42.10 44
112 1.10.565.10 39
110 1.10.555.10 38
204 1.20.58.390 38
574 3.30.710.10 37
794 4.10.280.10 37

library(plyr)
library(gplots)

cath_gpcr <- subset(cath_all, cath_id == "1.20.1070.10")

gpcr_matrix <- daply(cath_gpcr, "stable_id", function(df) {
  r <- df$omega
  if (mean(r, na.rm=T) > 0.5) {
    print(paste(df$dataset[1], df$stable_id[1]))
  }
  # names(r) <- df$stable_id[1]
  r
})
sum(gpcr_matrix > 1, na.rm=T) / length(gpcr_matrix)

gpcr_matrix_clean <- gpcr_matrix[apply(gpcr_matrix, 1, function(x) { sum(is.na(x)) < 0.5*ncol(gpcr_matrix)  }), ]
my.fisher <- function(x) {
  m <- matrix(c(sum(x > 1.0, na.rm=T),
                length(x) - sum(x > 1.0, na.rm=T),
                sum(gpcr_matrix > 1.0, na.rm=T) - sum(x > 1.0, na.rm=T),
                length(gpcr_matrix)-length(x)-sum(gpcr_matrix > 1.0, na.rm=T)), nrow=2)
  fisher.test(m, alternative="greater")$p.value
}

sum_conserved <- apply(gpcr_matrix_clean, 2, function(x) { sum(x == 0, na.rm=T) })
plot(1:length(sum_conserved), sort(sum_conserved))
hist(sum_conserved, breaks=100)

sum(gpcr_matrix_clean > 1.0, na.rm=T)/length(gpcr_matrix)
sum_variable <- apply(gpcr_matrix_clean, 2, function(x) { sum(x > 1.0, na.rm=T) })
plot(1:length(sum_variable), sort(sum_variable))
hist(sum_variable, breaks=100)

##
## Barplot
##
site_types <- matrix(rbind(sum_conserved, nrow(gpcr_matrix_clean)-sum_conserved-sum_variable, sum_variable), nrow=3)

pdf("Cdh_barplot.pdf", height=7, width=9)
barplot(site_types, col=c("blue", "lightgray", "red"))
dev.off()

pdf("Cdh_barplot_sorted.pdf", height=7, width=9)
barplot(site_types[, order(site_types[3, ], decreasing=T)], col=c("blue", "lightgray", "red"))
dev.off()


pvals_variable <- apply(gpcr_matrix_clean, 2, my.fisher)
hist(pvals_variable, breaks=100)
pvals_variable_adj <- p.adjust(pvals_variable, method="BH")
sum(pvals_variable_adj < 0.001)

my.fisher.cons <- function(x) {
  m <- matrix(c(sum(x < 0.1, na.rm=T),
                length(x) - sum(x < 0.1, na.rm=T),
                sum(gpcr_matrix < 0.1, na.rm=T) - sum(x < 0.1, na.rm=T),
                length(gpcr_matrix)-length(x)-sum(gpcr_matrix < 0.1, na.rm=T)), nrow=2)
  fisher.test(m, alternative="greater")$p.value
}

pvals_conserved <- apply(gpcr_matrix_clean, 2, my.fisher.cons)
hist(pvals_conserved, breaks=100)
pvals_conserved_adj <- p.adjust(pvals_conserved, method="BH")
sum(pvals_conserved_adj < 0.05)

col_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

gpcr_matrix_nomissing <- gpcr_matrix
gpcr_matrix_nomissing[is.na(gpcr_matrix_nomissing)] <- 0
gpcr_matrix_nomissing[gpcr_matrix_nomissing > 5] <- 5
pca <- prcomp(t(gpcr_matrix_nomissing))
# summary(pca)

pdf("GPCR_PCA_sites_1.pdf", height=8, width=9)
plot(pca$x[,1],pca$x[,2], col=rainbow(length(pca$x[,1])))
text(pca$x[,1], pca$x[,2], colnames(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black")
dev.off()

pdf("GPCR_PCA_sites_2.pdf", height=8, width=9)
plot(pca$x[,2],pca$x[,3], col=rainbow(length(pca$x[,1])))
text(pca$x[,2], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 
dev.off()

pdf("GPCR_PCA_sites_3.pdf", height=8, width=9)
plot(pca$x[,1],pca$x[,3], col=rainbow(length(pca$x[,1])))
text(pca$x[,1], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 
dev.off()

pca <- prcomp((gpcr_matrix_nomissing))

pdf("GPCR_PCA_genes_1.pdf", height=8, width=9)
plot(pca$x[,1],pca$x[,2], col=rainbow(length(pca$x[,1])))
# text(pca$x[,1], pca$x[,2], colnames(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black")
dev.off()

pdf("GPCR_PCA_genes_2.pdf", height=8, width=9)
plot(pca$x[,2],pca$x[,3], col=rainbow(length(pca$x[,1])))
# text(pca$x[,2], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 
dev.off()

pdf("GPCR_PCA_genes_3.pdf", height=8, width=9)
plot(pca$x[,1],pca$x[,3], col=rainbow(length(pca$x[,1])))
# text(pca$x[,1], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 
dev.off()


pdf("Cdh_heatmap.pdf", height=7, width=9)
gpcr_matrix_clean[gpcr_matrix_clean > 1] <- 1
heatmap.2(gpcr_matrix_clean, scale="none", dendrogram="none", sepwidth=c(0, 0), trace="none",
          col=col_palette, na.color="black", Rowv=T, Colv=NULL, keysize=1, labRow=T, labCol=F)
dev.off()

gpcr_matrix[gpcr_matrix > 1] <- 1
heatmap.2(gpcr_matrix, scale="none", dendrogram="none", sepwidth=c(0, 0), trace="none",
          col=col_palette, na.color="black", Rowv=NULL, Colv=NULL, keysize=1, labRow=T, labCol=F)



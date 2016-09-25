library(gplots)
# Untangling the different GPCR correlations
# 

mystery_gpcr_cors <- subset(cath_cors_dist_gpcr, dist > 2 & dist < 3.8)
nrow(mystery_gpcr_cors)
mystery_gpcr <- unique(c(mystery_gpcr_cors$id_1, mystery_gpcr_cors$id_2))
sum(mystery_gpcr %in% annotated)
sum(mystery_gpcr %in% annotated_taste)

mystery_gpcr_smell <- subset(mystery_gpcr_cors, (id_1 %in% annotated) & (id_2 %in% annotated))
points(mystery_gpcr_smell$dist, mystery_gpcr_smell$cor, col="darkgreen")
nrow(mystery_gpcr_smell)
points()

paste(mystery_gpcr[!(mystery_gpcr %in% annotated)], collapse=" ")

x <- subset(cath_all, stable_id == "ENSP00000307259")$omega
y <- subset(cath_all, stable_id == "ENSP00000366506")$omega
plot(x, y)
lm.object <- lm(x ~ y)
abline(lm.object, col="red")
both <- x > 0 | y > 0
cor.test(x[boy])
plot(x[both], y[both])
cor.test(x[both], y[both], method="spearman")
plot(lm.object)

sum(x == 0, na.rm=T)
sum(y == 0, na.rm=T)

fisher.test(matrix(c(131, 204-131,
154-131, 313-131-(204-131)-(154-131)), nrow=2))

cath_gpcr <- subset(cath_all, cath_id == "3.40.50.300")
all_pairs <- combn(unique(cath_gpcr$stable_id), 2, simplify=F)
str(all_pairs)
all_pairwise <- ldply(all_pairs[1:20000], function(r) {
  data.frame(id_1=r[1], id_2=r[2], x=subset(cath_gpcr, stable_id == r[1])$omega, y=subset(cath_gpcr, stable_id == r[2])$omega)
})

all_pairwise <- subset(all_pairwise, x < 5 & y < 5)

x <- pmin(all_pairwise$x, all_pairwise$y)
y <- pmax(all_pairwise$x, all_pairwise$y)
plot(hexbin(x, y))
smoothScatter(all_pairwise$x, all_pairwise$y) 
smoothScatter(x, y) 


# 
gpcr_matrix <- daply(cath_gpcr, "stable_id", function(df) {
  r <- df$omega
  # names(r) <- df$stable_id[1]
  r
})

sum_conserved <- apply(gpcr_matrix, 2, function(x) { sum(x == 0, na.rm=T) })
plot(1:length(sum_conserved), sort(sum_conserved))

sum(gpcr_matrix > 1.0, na.rm=T)/length(gpcr_matrix)
sum_variable <- apply(gpcr_matrix, 2, function(x) { sum(x > 1.0, na.rm=T) })
plot(1:length(sum_variable), sort(sum_variable))

my.fisher <- function(x) {
  m <- matrix(c(sum(x > 1.0, na.rm=T),
                length(x) - sum(x > 1.0, na.rm=T),
                sum(gpcr_matrix > 1.0, na.rm=T) - sum(x > 1.0, na.rm=T),
                length(gpcr_matrix)-length(x)-sum(gpcr_matrix > 1.0, na.rm=T)), nrow=2)
  fisher.test(m, alternative="greater")$p.value
}

pvals_variable <- apply(gpcr_matrix, 2, my.fisher)
hist(pvals_variable, breaks=100)
pvals_variable_adj <- p.adjust(pvals_variable, method="BH")

my.fisher.cons <- function(x) {
  m <- matrix(c(sum(x < 0.1, na.rm=T),
                length(x) - sum(x < 0.1, na.rm=T),
                sum(gpcr_matrix < 0.1, na.rm=T) - sum(x < 0.1, na.rm=T),
                length(gpcr_matrix)-length(x)-sum(gpcr_matrix < 0.1, na.rm=T)), nrow=2)
  fisher.test(m, alternative="greater")$p.value
}

pvals_conserved <- apply(gpcr_matrix, 2, my.fisher.cons)
hist(pvals_conserved, breaks=100)
pvals_conserved_adj <- p.adjust(pvals_conserved, method="BH")
sum(pvals_conserved_adj < 0.05)


## Find the most variable positions
CV <- function(mean, sd){
  (sd/mean)
}
pos_sd <- apply(gpcr_matrix, 2, function(x) { sd(x, na.rm=T) })
pos_mean <- apply(gpcr_matrix, 2, function(x) { mean(x, na.rm=T) })
pos_cv <- CV(pos_mean, pos_sd)
hist(CV(pos_mean, pos_sd), breaks=100)
pos_cv[order(pos_cv, decreasing=T)][1:5]
pos_cv[order(pos_cv)][1:5]

hist(gpcr_matrix[gpcr_matrix[, 301] < 10, 301], breaks=100)
hist(gpcr_matrix[gpcr_matrix[, 162] < 10, 162], breaks=100)

hist(gpcr_matrix[, 311], breaks=100)
hist(gpcr_matrix[, 6], breaks=100)
hist(gpcr_matrix[, 137], breaks=100)

image(t(gpcr_matrix))
gpcr_matrix_clean <- gpcr_matrix[apply(gpcr_matrix, 1, function(x) { sum(is.na(x)) < 0.5*ncol(gpcr_matrix)  }), ]
gpcr_matrix_clean[gpcr_matrix_clean > 1] <- 1
col_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

gpcr_matrix_nomissing <- gpcr_matrix[, apply(gpcr_matrix, 2, function(x) { sum(is.na(x)) < 20 })]
gpcr_matrix_nomissing[is.na(gpcr_matrix_nomissing)] <- 0
gpcr_matrix_nomissing[gpcr_matrix_nomissing > 5] <- 5

pca <- prcomp((gpcr_matrix_nomissing))
summary(pca)

plot(pca$x[,1],pca$x[,2], col=rainbow(length(pca$x[,1])))
text(pca$x[,1], pca$x[,2], colnames(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black")

plot(pca$x[,2],pca$x[,3], col=rainbow(length(pca$x[,1])))
text(pca$x[,2], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 

plot(pca$x[,1],pca$x[,3], col=rainbow(length(pca$x[,1])))
text(pca$x[,1], pca$x[,3], 1:ncol(gpcr_matrix_nomissing), cex=0.7, pos=4, col="black") 

# Plot the same as a function of mean omega?

biplot(pca)


library(devtools)
install_github("vqv/ggbiplot")

library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)


pdf("GPCR_pseudoalignment.pdf", height=10, width=13)
heatmap.2(gpcr_matrix_clean, scale="none", dendrogram="none", sepwidth=c(0, 0), trace="none",
          col=col_palette, na.color="black", Rowv=T, Colv=NULL, keysize=1, labRow=T, labCol=F)
dev.off()

mean_omega <- apply(gpcr_matrix, 1, function(x) { mean(x, na.rm=T) })
mean_omega[order(mean_omega, decreasing=T)][1:20]

plot(NA, axes=FALSE, xlab="", ylab="")

boxplot(gpcr_matrix, ylim=c(0, 2))
gpcr_matrix[gpcr_matrix[gpcr_matrix > 20]]

hist(gpcr_matrix)
hist(apply(gpcr_matrix, 2, function(x) { sd(x)/mean(x) }))

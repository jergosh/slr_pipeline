library(GenomicRanges)
library(rtracklayer)
library(plyr)
library(adehabitat)
library(biomaRt)

setwd("~/Documents/projects/slr_pipeline")

Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- 1-pchisq(Xsq, df = 2*length(p))
  return(c(Xsq = Xsq, p.value = p.val))
}

Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c(Z = Z, p.value = p.val))
}

# slr_all <- read.table("data/slr_all.tab", header=T, sep="\t", comment="", stringsAsFactors=F)
slr_all <- read.table("data/slr_all_200716.tab", header=T, sep="\t", comment="", stringsAsFactors=F)
slr_all$upper <- as.numeric(slr_all$upper)
slr_raw <- slr_all
# For the 3.1.0 fail
# slr_all <- read.table("data/slr_all.tab", header=T, sep="\t", comment="", stringsAsFactors=F,
#                       colClasses=c("integer", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "character",
#                                    "character", "integer", "numeric", "character", "numeric"))
head(slr_all)

summary(slr_all$Omega)
plot(density(slr_all$Omega, na.rm=T), xlim=c(0,5))
hist(slr_all$Omega, breaks=100)

slr_filtered <- subset(slr_all, (Note != "Single char") & (Note != "! Single cha") & (Note != "All gaps") & (Note != "0 Single cha"))

# Significance threshold
thr = 0.05
 
slr_all <- subset(slr_all, is.finite(slr_all$upper))
slr_all$Pval[slr_all$Omega < 1] <- runif(sum(slr_all$Omega < 1))
# Divide by 2
slr_all$Pval[slr_all$Omega > 1] <- slr_all$Pval[slr_all$Omega > 1] / 2

p_adjusted = p.adjust(slr_all$Pval, "BH")
slr_all$Adj.Pval <- p_adjusted
sum(slr_all$Adj.Pval < thr)
pos_sel <- subset(slr_all, Omega > 1 & Adj.Pval < thr)


sum(as.numeric(as.character(slr_all$lower)) > 1, na.rm=T)

all_ids <- unique(slr_all$stable_id)
write.table(all_ids, file="data/all_ids.tab", quote=F, col.names=F, row.names=F)

# Ens annotation
# ens_annotation <- import.gff("~/Downloads/Homo_sapiens.GRCh37.73.gtf", format="gtf")
# unique(elementMetadata(ens_annotation)[, "type"])

# p_adjusted = p.adjust(pos_sel$Pval, "BH")
# pos_sel[, "Adj.Pval"] <- p_adjusted
# pos_sel_idx <- with(slr_all, Omega > 1 & !(lower == 0 & upper == 0))
# slr_all[pos_sel_idx, "Adj.Pval"] <- p_adjusted
# pos_sel = subset(slr_all, Omega > 1 & !(lower == 0 & upper == 0))
# 
# neg_sel = subset(slr_all, Omega < 1 & !(lower == 0 & upper == 0))
# neg_Pvals = 1 - pchisq(neg_sel$LRT_Stat, df=1)
# neg_p_adjusted <- p.adjust(neg_Pvals, "BH")
# neg_sel[, "Adj.Pval"] <- neg_p_adjusted
# neg_sel_idx <- with(slr_all, Omega < 1 & !(lower == 0 & upper == 0))
# slr_all[neg_sel_idx, "Adj.Pval"] <- neg_p_adjusted

# The number of sites under positive selection per gene.
pos_per_gene <- daply(pos_sel, .variables="stable_id", .fun=function(df) { sum(df[, "Adj.Pval"] < thr) })
hist(pos_per_gene, breaks=30)

top_genes <- tail(names(sort(pos_per_gene)), n=200)

pval.dist <- daply(subset(slr_all, stable_id %in% top_genes), "stable_id", function(x) {wawotest(x[, "Omega"] < 1 & x[, "Adj.Pval"] < thr)["p"] })
hist(pval.dist, breaks=100)

png("omega_dist.png", width=500, height=400)
par(mar=c(5.1,5.1,4.1,2.1))
hist(slr_filtered$Omega, xlim=c(0, 3), ylim=c(0, 20), breaks=10000, freq=F,
     main=expression(paste("Distribution of ", omega)), xlab=expression(omega), cex.lab=1.5, cex.main=1.5)
dev.off()

library(plotrix)
hdata <- hist(slr_filtered$Omega, xlim=c(0, 3), breaks=seq(0, ceiling(max(slr_filtered$Omega)), 0.05))

gap.plot(density(slr_filtered$Omega, from=0, to=3), gap=c(2, 27), xlim=c(0, 3))
pdf("omega_hist_broken.pdf", height=7, width=8)
gap.barplot(hdata$counts, gap=c(1e6,5.55e6), xtics=hdata$breaks, ytics=c(0, 5.6e5, 5.6e6),
            xlim=c(0, 3), col=c(rep(rgb(0, 0, 1, 0.4), 20), "white", rep(rgb(1, 0, 0, 0.4), 80)), xlab=expression(omega), cex.lab=1.8, cex.axis=1, cex.main=1.5, cex.sub=1.5,
            ylab="count",main=expression(paste("Distribution of ", omega)))
abline(v=1)
dev.off()
sum(slr_filtered$Omega < 1)/nrow(slr_filtered)
sum(slr_filtered$Omega > 1)/nrow(slr_filtered)
# Distribution of \omega point estimates + point estimates 
# What does Greg's plot look like? -- doesn't seem to be in the thesis?!

slr_nonzero <- subset(slr_filtered, Omega != 0)
slr_nonzero$lower <- as.numeric(slr_nonzero$lower)
slr_nonzero$upper <- as.numeric(slr_nonzero$upper)
slr_nonzero <- subset(slr_nonzero, !is.na(lower) & !is.na(upper) & !is.infinite(upper))

pdf("omega_summary.pdf", width=8, height=8)
par(mfrow=c(2, 2))
hist(slr_nonzero$Omega, breaks=2000, xlim=c(0, 6), col=rgb(0.0, 0.0, 0.0, 0.3),
     main=expression(paste("Histogram of ", omega)), xlab=expression(omega))
hist(slr_nonzero$lower, breaks=1000, xlim=c(0, 6), add=T, col=rgb(0.0, 0.0, 1.0, 0.3))
hist(slr_nonzero$upper, breaks=2000, xlim=c(0, 6), add=T, col=rgb(1.0, 0.0, 0.0, 0.3))
## Cumulative plots as well
plot(ecdf(slr_nonzero$Omega), xlim=c(0, 6), main=expression(paste("CDF of ", omega, " estimates")))
plot(ecdf(slr_nonzero$lower), col="blue", add=T)
plot(ecdf(slr_nonzero$upper), col="red", add=T)
legend("bottomright", legend=c(expression(omega), "CI lower", "CI upper"), fill=c("black", "blue", "red"))

##
## No. of significant sites as a function of FDR
## 
hist(subset(slr_nonzero, Omega > 1)$Pval, breaks=100)
thrs <- seq(0.0, 1.0, 0.05)
sums <- vapply(thrs, function(thr) sum(slr_nonzero$Adj.Pval < thr & slr_nonzero$Omega > 1),
               FUN.VALUE=0.0)

# This is not very informative as most p-values get pushed to 1 anyway.

plot(thrs, sums, main="Number of significant sites", xlab="FDR", ylab="Count")
par(mfrow=c(1, 1))
dev.off()

## Whole-gene analysis?
## Is it enough to just multiply p-values under < 0.05 
slr_nonzero_pos <- subset(slr_nonzero, Omega > 1.0)

fisher.pvals <- ddply(slr_nonzero_pos, "stable_id", function(df) {Fisher.test(df$Pval)})
stouffer.pvals <- ddply(slr_nonzero_pos, "stable_id", function(df) {Stouffer.test(df$Pval)})
min.pvals <- ddply(slr_nonzero_pos, "stable_id", function(df) {c(p.value=min(p.adjust(df$Pval, method="BH")))})

N <- length(unique(slr_all$stable_id))
fisher.sums <- vapply(thrs, function(thr) sum(p.adjust(fisher.pvals$p.value, method="BH", n=N) < thr),
               FUN.VALUE=0.0)

stouffer.sums <- vapply(thrs, function(thr) sum(p.adjust(stouffer.pvals$p.value, method="BH", n=N) < thr),
                      FUN.VALUE=0.0)

min.sums <- vapply(thrs, function(thr) sum(p.adjust(min.pvals$p.value, method="BH", n=N) < thr),
                      FUN.VALUE=0.0)

##
## Have the adjusted p-values by threshold as well?
## Maybe put the distribution + 
## What about the genes where there are no sites with \omega > 1?
##

pdf("PSG.pdf", height=11, width=8)
par(mfrow=c(3,2))
hist(fisher.pvals$p.value, breaks=100, main="Gene p-values (Fisher's method)", xlab="p-value")
plot(thrs, fisher.sums)
hist(stouffer.pvals$p.value, breaks=100, main="Gene p-values (Stouffer's method)", xlab="p-value")
plot(thrs, stouffer.sums)
hist(min.pvals$p.value, breaks=100, main="Gene p-values (BH-corrected minimum)", xlab="p-value")
plot(thrs, min.sums)
par(mfrow=c(1,1))
dev.off()

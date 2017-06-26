library(GenomicRanges)
library(rtracklayer)

library(adehabitat)
library(biomaRt)
library(plyr)

setwd("~/Documents/projects/slr_pipeline")

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

slr_all <- subset(slr_all, (Note != "Single char") & (Note != "! Single cha") & (Note != "All gaps") & (Note != "0 Single cha"))
slr_all$ens_pos <- slr_all$human_idx

all_ids <- unique(slr_all$stable_id)
write.table(all_ids, file="data/all_ids.tab", quote=F, col.names=F, row.names=F)

# Significance threshold
thr = 0.05
 
slr_all <- subset(slr_all, is.finite(slr_all$upper))
# slr_all$Pval[slr_all$Omega < 1] <- runif(sum(slr_all$Omega < 1))
# Divide by 2
# slr_all$Pval[slr_all$Omega > 1] <- slr_all$Pval[slr_all$Omega > 1] / 2

# slr_all <- subset(slr_all, Pval <= 1)

p_adjusted <- p.adjust(slr_all$Pval, "BH")
slr_all$Adj.Pval <- p_adjusted
sum(slr_all$Adj.Pval < thr & slr_all$Omega > 1)

pos_sel <- subset(slr_all, Omega > 1 & Adj.Pval < thr)
sum(as.numeric(as.character(slr_all$lower)) > 1, na.rm=T)

## Disorder predictions
disorder_pred <- read.table("data/iupred.tab", sep="\t", header=F, stringsAsFactors=F)
head(disorder_pred)
colnames(disorder_pred) <- c("stable_id", "ens_pos", "iupred")

slr_all <- join(slr_all, disorder_pred, match="first")
sum(is.na(slr_all$iupred))
slr_all$disorder <- slr_all$iupred > 0.5
hist(slr_all$iupred, breaks=100)
table(slr_all$disorder)

sum(is.na(slr_all$disorder))

sum(slr_all$Adj.Pval < 0.05 & slr_all$Omega > 1.0 & slr_all$disorder == TRUE, na.rm=T)/sum(slr_all$disorder == TRUE, na.rm=T)
sum(slr_all$Adj.Pval < 0.05 & slr_all$Omega > 1.0 & slr_all$disorder == FALSE, na.rm=T)/sum(slr_all$disorder == FALSE, na.rm=T)

mean(subset(slr_all, disorder == TRUE)$Omega)
mean(subset(slr_all, disorder == FALSE)$Omega)

ggplot(slr_all, aes(y=fr_aligned, x=disorder)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("") +
  ylab("Fraction aligned") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()


slr_all <- ddply(slr_all, "stable_id", function(df) {
  df[1:(nrow(df)-1), ]
})

aln_stats <- read.table("data/aln_stats.tab", header=F, stringsAsFactors=F, sep="\t")
colnames(aln_stats) <- c("dataset", "Site", "n_aligned", "fr_aligned", "entropy")
slr_stats <- join(slr_all, aln_stats)
slr_all <- slr_stats

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

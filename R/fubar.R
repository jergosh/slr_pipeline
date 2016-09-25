library(plyr)
library(dplyr)
library(modeest)
library(numDeriv)
library(ggplot2)
setwd("~/Documents/projects/slr_pipeline/")

dNdSplot <- function(df){
  ord = order(df$alpha, decreasing=T)
  plot(df$alpha[ord], main=df$ID[1], ylim=c(0, max(df$alpha, df$beta)+0.1))
  points(df$beta[ord], col="indianred")
  abline(h=mean(df$alpha))
  abline(h=mean(df$beta), col="indianred")
}

fubar_table <- read.table("data/fubar.tab", header=T)
# colnames(fubar_table) <- c("ID", "pos", "alpha", "beta", "diff", "Palphasmaller", "Palphagreater", "BayesFactor", "PSRF", "Neff", "Nseqs")
# fubar_table_filtered <- subset(fubar_table, Nseqs > 1)

fubar_table_filtered <- ddply(fubar_table, .variables=c("ID"), .fun=function(df) {
  mode <- mlv(df$alpha, method = "mfv")$M
  subset(df, abs(alpha - mode) > 1e-6)
})

dim(fubar_table_filtered)
  
fubar_summary <- ddply(fubar_table, .variables=c("ID"), .fun=function(df)
{ data.frame(mean_alpha=mean(df$alpha), mean_beta=mean(df$beta), sd_alpha=sd(df$alpha), sd_beta=sd(df$beta)) })

plot(density(fubar_summary$mean_alpha))
lines(density(fubar_summary$mean_beta), col="indianred")

# I suppose the overall substitution rate is 1
mean(fubar_summary$mean_alpha)

plot(density(fubar_summary$sd_alpha))
lines(density(fubar_summary$sd_beta), col="indianred")

pdf("dS_plots.pdf")
daply(fubar_table_filtered, .variables=c("ID"), .fun=dNdSplot)
dev.off()


#
# Is it better to estimate dS genewise of sitewise?
#

# beta - alpha as a function of alpha - mean(alpha)
alpha_diff <- ddply(fubar_table, c("ID"), function(df) {
    df$alpha_diff <- df$alpha - mean(df$alpha)
    df
})
alpha_diff$dNdS <- alpha_diff$beta - alpha_diff$alpha

ggplot(subset(alpha_diff, alpha_diff < 1), aes(x=alpha_diff, y=dNdS)) +
  geom_hex()

pdf("density2d.pdf")
ggplot(subset(alpha_diff, alpha_diff < 1), aes(x=alpha_diff, y=dNdS)) +
  stat_density2d(aes(fill = ..level..), geom="polygon")
dev.off()

ggplot(alpha_diff, aes(x=alpha_diff, y=dNdS))
  

dS_means <- dlply(fubar_table_filtered, .variables=c("ID"), .fun=function(df) { mean(df$alpha) })
fubar_table_pos <- subset(fubar_table_filtered, Palphasmaller > 0.95)
dS_dist <- daply(fubar_table_pos, .variables=c("ID", "pos"), .fun=function(r) { r$alpha[1] > dS_means[[r$ID[1]]] })
dS_sums <- table(dS_dist)
# above_fraction_pos <- dS_sums$

# overall_dist <- daply(fubar_table_filtered[1:1000, ], .variables=c("ID", "pos"), .fun=function(r) { r$alpha[1] > dS_means[[r$ID[1]]] })
# overall_dist <- apply(fubar_table_filtered, 1, function(r) { r["alpha"] > dS_means[[ r["ID"] ]] })

# Another go with dplyr
is_above <- function(alpha, ID) { sum(alpha > dS_means[[ID[1]]]) }
is_below <- function(alpha, ID) { sum(alpha < dS_means[[ID[1]]]) }

above_sums <- fubar_table_filtered %.% group_by(ID) %.%
  summarise(above = is_above(alpha, ID))
below_sums <- fubar_table_filtered %.% group_by(ID) %.%
  summarise(below = is_below(alpha, ID))
above_fraction <- sum(above_sums$above) / nrow(fubar_table_filtered)
below_fraction <- sum(below_sums$below) / nrow(fubar_table_filtered)

# What about confidence intervals?
above_sums_pos <- fubar_table_pos %.% group_by(ID) %.%
  summarise(above = is_above(alpha, ID))
below_sums_pos <- fubar_table_pos %.% group_by(ID) %.%
  summarise(below = is_below(alpha, ID))

above_fraction_pos <- sum(above_sums_pos$above) / sum(above_sums$above)
below_fraction_pos <- sum(below_sums_pos$below) / sum(below_sums$below)

pdf("above_below.pdf")
barplot(matrix(c(above_fraction, below_fraction), nrow=2), beside=T, width=0.5, xlim=c(0, 2), ylim=c(0, 1),
        main="Fraction above/below mean dS", names.arg=c("Above", "Below"))
dev.off()
pdf("above_below_pos.pdf")
barplot(matrix(c(above_fraction_pos, below_fraction_pos), nrow=2), beside=T, width=0.5, xlim=c(0, 2), ylim=c(0, 0.002),
        main="Fraction under pos. selection", names.arg=c("Above mean", "Below mean"))
dev.off()

# Things to try:
# - Overlap with SLR
# -- overlap in each case (under/above the mean)
# - counts normalised by number of sites under/below the mean



#
# Scatterplot of dN vs. dS?
# can do all in one plot!

library(hexbin)
pdf("FUBAR_hexbin.pdf")
plot(hexbin(log(fubar_table$alpha), log(fubar_table$beta)))
abline(a=1, b=0)
dev.off()

library(lattice)
panel.hexbinplot(log(fubar_table$alpha), log(fubar_table$beta), aspect = 1, bins=50, 
           xlab = expression(alpha), ylab = expression(beta), 
           # style = "nested.centroids",
           panel = function(...) {
             panel.hexbinplot(...)
             panel.abline(h=0)
           })

head(fubar_table)

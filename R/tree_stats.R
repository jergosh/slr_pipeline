library(gridExtra)
library(ggplot2)

setwd("~/Documents/projects/slr_pipeline/")

trees_raw_eutheria <- read.csv("trees_full_mammals.csv", header=T)
trees_raw <- read.csv("tree_stats_raw.csv", header=T)
trees_raw_63 <- read.csv("tree_stats_raw_old.csv", header=T)
trees_old <- read.csv("tree_stats.csv", header=T)
trees_greg <- read.csv("tree_stats_greg.csv", header=T)
trees_greg_63 <- read.csv("tree_stats_greg_63.csv", header=T, stringsAsFactors=F)
baseline_greg_63 <- read.csv("baseline_stats_greg_63.csv", header=T, stringsAsFactors=F)
human_baseline_greg_63 <- subset(baseline_greg_63, X..of.Human.seqs == 1)
good_baseline_greg_63 <- subset(human_baseline_greg_63, X..of.Species >= 21)
trees_new <- read.csv("tree_stats_new.csv", header=T)
trees_new_04 <- read.csv("tree_stats_new_0.4.csv", header=T)

tree_summary <- function(df) {
  data.frame(`Number of trees`=nrow(df), `Mean leaves`=mean(df$X..of.Leaves), `Mean species`=mean(df$X..of.Species), `Human sequences`=sum(df$X..of.Human.seqs), `0 human`=sum(df$X..of.Human.seqs == 0),
             `1 human`=sum(df$X..of.Human.seqs == 1), `2+ human`=sum(df$X..of.Human.seqs > 1), `Mean paralog fraction`=mean(df$X..of.Paralogs/df$X..of.Leaves))
}

tree_stats <- function(stats, prefix="") {
  hist(stats$X..of.Leaves, breaks=max(stats$X..of.Leaves)/7, xlim=c(0, max(400)), main=paste("Histogram of number of taxa in", prefix),
       xlab="Number of taxa")  
  hist(subset(stats, X..of.Leaves <= 100)$X..of.Leaves, breaks=100, xlim=c(0, 100), main=paste("Histogram of number of taxa in", prefix),
       xlab="Number of taxa")
  hist(stats$X..of.Paralogs/stats$X..of.Leaves, breaks=100, main=paste("Histogram of paralog fraction in", prefix),
       xlab="Paralog fraction")
  
  print(sum(stats$X..of.Paralogs/stats$X..of.Leaves > 0.2)/nrow(stats))
  print(nrow(stats))
  print(sum(stats$X..of.Human.seqs))
  print(summary(stats$X..of.Leaves))
  
  hist_top <- ggplot() + geom_histogram(data=stats, aes(X..of.Leaves))
  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  scatter <- ggplot()+geom_point(data=stats, aes(x=X..of.Leaves, y=X..of.Paralogs/X..of.Leaves))
  hist_right <- ggplot() + geom_histogram(data=stats, aes(X..of.Paralogs/X..of.Leaves)) +
    coord_flip()
  
  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4),
               top=prefix)
  
  ggplot(stats,aes(x=X..of.Leaves, y=X..of.Paralogs/X..of.Leaves))+
    stat_density2d(aes(alpha=..level..), geom="polygon") +
    scale_alpha_continuous(breaks=seq(0,0.2,by=0.025))+ # limits=c(0,0.2),
    #geom_point(colour="red",alpha=0.02)+
    ggtitle(prefix) +
    theme_bw()
  # k <- with(stats, MASS:::kde2d(X..of.Leaves, X..of.Paralogs/X..of.Leaves, n=500))
  # filled.contour(k)
}

tree_hist <- function(df, n_species, limit, yend, main) {
  df$X..of.Leaves[df$X..of.Leaves > limit] <- limit
  df <- subset(df, X..of.Leaves > 1)

  vlines <- data.frame(y=0, yend=yend, x=seq(n_species, limit, n_species), xend=seq(n_species, limit, n_species))
  ggplot(df, aes(x=X..of.Leaves)) +
    geom_histogram(bins=limit/2) +
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend), linetype="longdash", data=vlines, colour="grey") +
    ggtitle(main) +
  xlab("Number of leaf nodes") +
    ylab("Count") +
    theme_bw()
}

n_human_hist <- function(df, limit, yend, main) {
  df$X..of.Human.seqs[df$X..of.Human.seqs > limit] <- limit
  df <- subset(df, X..of.Human.seqs > 1)
  
  # vlines <- data.frame(y=0, yend=yend, x=seq(n_species, limit, n_species), xend=seq(n_species, limit, n_species))
  p <- ggplot(df, aes(x=X..of.Human.seqs)) +
    geom_histogram(aes(y=..density..), bins=max(df$X..of.Human.seqs)) +
    # geom_segment(aes(x=x, xend=xend, y=y, yend=yend), linetype="longdash", data=vlines, colour="grey") +
    ggtitle(main) +
    xlab("Number of human genes") +
    ylab("Count") +
    theme_bw()
  
  pg <- ggplot_build(p)
  sc <- as.numeric(max(pg$data[[1]]$x))

  p + stat_bin(aes_string(y=cumsum(..count..)/sum(..count..))/(sc), geom="line", color="grey")
}

tree_summary(trees_raw)
tree_summary(trees_raw_63)
tree_summary(trees_greg_63)
tree_summary(trees_greg)
tree_summary(trees_new)

trees_greg$fraction.paralogs <- trees_greg$X..of.Paralogs/trees_greg$X..of.Leaves
mean(trees_greg$fraction.paralogs[order(trees_greg$fraction.paralogs, decreasing=T)[1:100]])

trees_new$fraction.paralogs <- trees_new$X..of.Paralogs/trees_new$X..of.Leaves
mean(trees_new$fraction.paralogs[order(trees_new$fraction.paralogs, decreasing=T)[1:100]])

sum(trees_new$fraction.paralogs < 0.1)
sum(trees_greg$fraction.paralogs < 0.1)

## Number of duplications
hist(trees_new$X..of.duplications)

trees_new_pf <- trees_new[order(trees_new$fraction.paralogs, decreasing=T)[1:200], ]
trees_new_pf <- subset(trees_new, fraction.paralogs > 0.5)
hist(trees_new_pf$X..of.duplications, breaks=500)

# Histogram / eCDF for the number of human genes


pdf("results/chapter3/trees_raw_78.pdf", width=7, height=3)
tree_hist(trees_raw, 39, 200, 1400, "Distribution of raw tree sizes (Ensembl 78)")
dev.off()

pdf("results/chapter3/trees_raw_63.pdf", width=3.5, height=3)
tree_hist(trees_raw_63, 35, 200, 1300, "Distribution of raw tree sizes (Ensembl 63)")
dev.off()

pdf("results/chapter3/trees_new.pdf", width=4.5, height=4)
tree_hist(trees_new, 39, 80, 3400, "Distribution of split tree sizes (new scheme)")
dev.off()

pdf("results/chapter3/trees_greg.pdf", width=4.5, height=4)
tree_hist(trees_greg, 39, 80, 3400, "Distribution of split tree sizes (old scheme)")
dev.off()



pdf("results/chapter3/n_human_78.pdf", width=3.5, height=3)
n_human_hist(trees_raw, 200, 1400, "Human genes per tree (Ens 78)")
dev.off()
sum(trees_raw$X..of.Human.seqs)
summary(trees_raw$X..of.Human.seqs)

pdf("results/chapter3/n_human_63.pdf", width=3.5, height=3)
n_human_hist(trees_raw_63, 200, 1400, "Human genes per tree (Ens 63)")
dev.off()
summary(trees_raw_63$X..of.Human.seqs)

df <- trees_greg
limit <- 120
yend <- 70

df$X..of.Leaves[df$X..of.Leaves > limit] <- limit
  df <- subset(df, X..of.Leaves > 1)
  
  vlines <- data.frame(y=0, yend=yend, x=seq(n_species, limit, n_species), xend=seq(n_species, limit, n_species))
  ggplot(df, aes(x=X..of.Leaves)) +
    geom_histogram(bins=limit/2) +
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend), linetype="longdash", data=vlines, colour="grey") +
    theme_bw()

tree_hist(trees_greg, 39, 200, 2300)

tree_hist(trees_new, 39, 200, 1300)



pdf("results/tree_stats/unsplit_eutheria.pdf")
tree_stats(trees_raw_eutheria, prefix="unsplit trees (Eutheria)")
dev.off()
pdf("results/tree_stats/unsplit.pdf")
tree_stats(trees_raw, prefix="unsplit trees")
dev.off()
pdf("results/tree_stats/gregj.pdf")
tree_stats(trees_greg, prefix="split trees (GregJ)")
dev.off()
pdf("results/tree_stats/gregj_63.pdf")
tree_stats(trees_greg_63, prefix="split trees (GregJ Ens 63)")
dev.off()
pdf("results/tree_stats/baseline_gregj_63.pdf")
tree_stats(baseline_greg_63, prefix="split trees (GregJ Ens 63)")
dev.off()
pdf("results/tree_stats/baseline_gregj_63.pdf")
tree_stats(good_baseline_greg_63, prefix="split trees (GregJ Ens 63)")
dev.off()
pdf("results/tree_stats/old.pdf")
tree_stats(trees_old)
dev.off()
pdf("results/tree_stats/new.pdf")
tree_stats(trees_new)
dev.off()
pdf("results/tree_stats/new_0.4.pdf")
tree_stats(trees_new_04)
dev.off()


mean(trees_new$X..of.Paralogs/trees_new$X..of.Leaves)
mean(trees_greg$X..of.Paralogs/trees_greg$X..of.Leaves)
mean(trees_new_04$X..of.Paralogs/trees_new_04$X..of.Leaves)


sum(trees_new$X..of.Paralogs == 0)
sum(trees_greg$X..of.Paralogs == 0)

sum(trees_new$X..of.Paralogs/trees_new$X..of.Leaves > 0.1)
sum(trees_greg$X..of.Paralogs/trees_greg$X..of.Leaves > 0.1)

mean(trees_new$X..of.Leaves)
mean(trees_greg$X..of.Leaves)

sum(trees_new$X..of.Human.seqs)
sum(trees_new_04$X..of.Human.seqs)
sum(trees_greg$X..of.Human.seqs)
sum(subset(trees_raw, X..of.Leaves > 39)$X..of.Human.seqs)

# What's happening in the raw trees 
nrow(trees_raw)

## Where do the human IDs go?
table(trees_greg$X..of.Human.seqs)
table(trees_new$X..of.Human.seqs)
table(trees_new_04$X..of.Human.seqs)

sum(trees_greg$X..of.Human.seqs > 1)
sum(trees_new$X..of.Human.seqs > 1)

# Subset the raw trees to work out which human IDs get lost

# What would be the most indicative of cases where Greg's scheme is better are trees where there is only one human sequence which is not present 
# under my scheme.


# TODO Put the marginal distribution 
# TODO Replace this with a 2d density plot
# Really need a density plot here...
plot(trees_new$X..of.Leaves, trees_new$X..of.Paralogs/trees_new$X..of.Leaves)


positive_counts <- ddply(pdb_master, c("stable_id", "pdb_id"), function(df) {
  data.frame(count=sum(df$omega > 1.0))
})
positive_counts[order(positive_counts$count, decreasing=T)[1:50], ]


# Omega / fraction under positive selection as a function of paralogs
# (ideally the entire density)

# Cluster top 10 residues or whatever?!

# Percentage of amino-acid identity as a local measure of quality of structural mapping

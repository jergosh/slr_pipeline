library(ggplot2)
library(stringr)
library(plyr)
# Need to join or break down the data depending on the various features of the trees

trees_old$dataset <- str_extract(trees_old$Directory, "[0-9]+_[0-9]+")
trees_old$paralog_frac <- trees_old$X..of.Paralogs / trees_old$X..of.Leaves
trees_old$n_taxa <- trees_old$X..of.Leaves
trees_old$treelen <- trees_old$Length

tree_info <- trees_old[, c("dataset", "paralog_frac", "n_taxa", "treelen")]

slr_stats <- join(slr_all, tree_info)
slr_stats_nonzero <- subset(slr_stats, Omega != 0.0)

taxa_deciles <- quantile(slr_stats$n_taxa, probs=seq(0.1, 1.0, 0.1))
slr_stats$taxa_bin <- sapply(slr_stats$n_taxa, function(n) { sum(n > taxa_deciles) })
# slr_stats$paralog_bin <- 

table(subset(slr_stats, taxa_bin == 9)$n_taxa)
for (i in 0:9) {
  obs <- subset(slr_stats, taxa_bin == i)
  print(sum(obs$Omega == 0.0) / nrow(obs))
}

for (i in 0:9) {
  obs <- subset(slr_stats_nonzero, taxa_bin == i)
  print(mean(obs$Omega))
}

for (i in 0:9) {
  obs <- subset(slr_stats_nonzero, taxa_bin == i)
  print(sum(obs$Omega > 1.0) / nrow(obs))
}

ggplot(slr_stats_nonzero, aes(x=Omega, colour=as.factor(taxa_bin))) +
  coord_cartesian(xlim=c(0, 2), ylim=c(0, 6)) +
  geom_density()

# Does the ortholog fraction paint a clearer picture?
plot(ecdf(tree_info$paralog_frac))
# paralog_deciles <- quantile(slr_stats$paralog_frac, probs=seq(0.1, 1.0, 0.1))
paralog_breaks <- seq(0.0, 0.2, 0.02)
slr_stats$paralog_bin <- sapply(slr_stats$paralog_frac, function(n) { sum(n > paralog_breaks) })

ggplot(slr_stats_nonzero, aes(x=Omega)) +
  # coord_cartesian(xlim=c(0, 2), ylim=c(0, 6)) +
  geom_density()

ggplot(slr_stats_nonzero, aes(x=Omega, colour=as.factor(paralog_bin))) +
  # coord_cartesian(xlim=c(0, 2), ylim=c(0, 6)) +
  geom_density()

for (i in 0:9) {
  obs <- subset(slr_stats_nonzero, paralog_bin == i)
  print(sum(obs$Omega > 1.0) / nrow(obs))
}

# Tree length
treelen_deciles <- quantile(slr_stats$treelen, probs=seq(0.1, 1.0, 0.1))
slr_stats$treelen_bin <- sapply(slr_stats$treelen, function(n) { sum(n > treelen_deciles) })

for (i in 0:9) {
  obs <- subset(slr_stats, treelen_bin == i)
  print(sum(obs$Omega == 0.0) / nrow(obs))
}

# Is this a statement about fluctuating selective pressure (or lack there of)?
# Contra what people believe 
for (i in 0:9) {
  obs <- subset(slr_stats, treelen_bin == i)
  print(sum(obs$Omega > 1.0 & obs$Adj.Pval < 0.05) / nrow(obs))
}


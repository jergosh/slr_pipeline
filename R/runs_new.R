library(randtests)
library(adehabitatLT)
library(hexbin)
library(Hmisc)
require(cowplot)
library(ggplot2)
library(ggExtra)

make.plots <- function(df, title="") {
  for (thr in c(3, 5, 10)) {
    df.subset <- subset(df, nsites >= thr)
    hist(df.subset$p, breaks=1000, main=paste0("WaW p-value distribution (", nrow(df.subset), " ", title, " >=", thr," sites)"),
         xlab="", ylab="")
    hist(df.subset$p.adj, breaks=1000, main=paste0("WaW p-value distribution (", nrow(df.subset), " ", title, " >=", thr," sites)"),
         xlab="", ylab="")    
    
    plot(df.subset$len, df.subset$p, main=capitalize(title),
         xlab="Protein length", ylab="P-value")
    
    # df.subset.trimmed <- subset(df.subset, len > 100)
    # plot(hexbin(df.subset.trimmed$len, df.subset.trimmed$p), main=capitalize(title),
    #      xlab="Protein length", ylab="P-value")
    plot(df.subset.trimmed$len, df.subset.trimmed$p, main=paste(capitalize(title), "(length > 100)"),
        xlab="Protein length", ylab="P-value")
  }
}

# process.runs <- function(df, stat, op, stat_thr) {
#   opfun <- match.fun(FUN = op) 
#   v <- as.numeric(opfun(df[, stat], stat_thr) & df[, "Omega"] > 1) 
#   nsites <- sum(v)
#   if (nsites > 1) {
#     res <- runs.test(v, "left.sided", threshold=0.5)
#     data.frame(p=res$p.value, nsites=nsites, len=length(v))
#   } else {
#     data.frame()
#   }
# }

process.runs <- function(df, stat, op, stat_thr) {
  opfun <- match.fun(FUN = op)
  
  v <- as.numeric(opfun(df[, stat], stat_thr) & df[, "Omega"] > 1)
  nsites <- sum(v)

  if (nsites > 1) {
    # res <- wawotest(v, alter="greater")
    # data.frame(p=res["p"], nsites=nsites, len=length(v))
        
    res <- runs.test(v, alternative="left.sided", pvalue="exact", plot=F, threshold=0.5)
    data.frame(p=res$p.value, nsites=nsites, len=length(v))
  } else {
    data.frame()
  }
}

## NEW
runs.all_0.05 <- ddply(slr_structure_all, c("stable_id", "pdb_id", "pdb_chain"), process.runs, "Adj.Pval", "<", 0.05)
runs.all_0.1 <- ddply(slr_structure_all, c("stable_id", "pdb_id", "pdb_chain"), process.runs, "Adj.Pval", "<", 0.1)
runs.all_0.2 <- ddply(slr_structure_all, c("stable_id", "pdb_id", "pdb_chain"), process.runs, "Adj.Pval", "<", 0.2)
runs.all_0.5 <- ddply(slr_structure_all, c("stable_id", "pdb_id", "pdb_chain"), process.runs, "Adj.Pval", "<", 0.5)

runs.all_0.05 <- subset(runs.all_0.05, p == p)
runs.all_0.1 <- subset(runs.all_0.1, p == p)
runs.all_0.2 <- subset(runs.all_0.2, p == p)
runs.all_0.5 <- subset(runs.all_0.5, p == p)



pval_hist <- function(data) {
  h <- ggplot(data, aes(x=p)) +
    theme_minimal(base_size=8) +
    theme(plot.margin=margin(t=0.1, r=0.1, l=0.1, b=0.1, unit="cm"),
          panel.spacing=margin(0, 0, 0, 0),
          #a xis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x=element_line(size=0.25),
          axis.line.y=element_line(size=0.25)) +
    # ggtitle("") +
    xlab("P-value") +
    ylab("count") +
    geom_histogram(binwidth=0.01)
  
  h
}

p1 <- qqplot_mar(runs.all_0.05$p) #  & nsites > 10
save_plot("results/WaW/waw_0_05_qq.pdf", p1, base_aspect_ratio=1.1, base_height=6/2.54, base_width=6.5/2.54)

p2 <- qqplot_mar(runs.all_0.1$p) #  & nsites > 10
save_plot("results/WaW/waw_0_1_qq.pdf", p2, base_aspect_ratio=1.1, base_height=6/2.54, base_width=6.5/2.54)

p3 <- qqplot_mar(runs.all_0.2$p) #  & nsites > 10
save_plot("results/WaW/waw_0_2_qq.pdf", p3, base_aspect_ratio=1.1, base_height=6/2.54, base_width=6.5/2.54)

p4 <- qqplot_mar(subset(runs.all_0.5, p != 0)$p) #  & nsites > 10
save_plot("results/WaW/waw_0_5_qq.pdf", p4, base_aspect_ratio=1.1, base_height=6/2.54, base_width=6.5/2.54)


h1 <- pval_hist(runs.all_0.05)
save_plot("results/WaW/waw_0_05_hist.pdf", h1, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h2 <- pval_hist(runs.all_0.1)
save_plot("results/WaW/waw_0_1_hist.pdf", h2, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h3 <- pval_hist(runs.all_0.2)
save_plot("results/WaW/waw_0_2_hist.pdf", h3, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h4 <- pval_hist(runs.all_0.5)
save_plot("results/WaW/waw_0_5_hist.pdf", h4, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)


h1 <- pval_hist(subset(runs.all_0.05, nsites >= 5))
save_plot("results/WaW/waw_hist_0_05_filtered.pdf", h1, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h2 <- pval_hist(subset(runs.all_0.1, nsites >= 5))
save_plot("results/WaW/waw_hist_0_1_filtered.pdf", h2, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h3 <- pval_hist(subset(runs.all_0.2, nsites >= 5))
save_plot("results/WaW/waw_hist_0_2_filtered.pdf", h3, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)

h4 <- pval_hist(subset(runs.all_0.5, nsites >= 5))
save_plot("results/WaW/waw_hist_0_5_filtered.pdf", h4, base_aspect_ratio=1.1, base_height=5/2.54, base_width=6.5/2.54)


## Case studies
runs.all_0.2[order(runs.all_0.2$p)[1:4], ]



## Old way -- likely wouldn't work with marginal distributions
# old_par <- par(mfrow=c(2, 2), oma=c(2.0, 2.0, 3.0, 0.0) + 1.0,
#                mar = c(2.0, 2.0, 0.0, 0.0) + 1.1)
# ggd.qqplot(subset(runs.all_0.05, p == p & nsites > 10)$p, main="Significance threshold 0.05")
# ggd.qqplot(subset(runs.all_0.1, p == p & nsites > 10)$p, main="Significance threshold 0.1")
# ggd.qqplot(subset(runs.all_0.2, p == p & nsites > 10)$p, main="Significance threshold 0.2")
# ggd.qqplot(subset(runs.all_0.5, p == p & nsites > 5)$p, main="Significance threshold 0.5")
# par(old_par)
# qq(subset(runs.all_0.05, p == p & nsites > 10)$p)


runs.all <- ddply(slr_all, c("stable_id"), process.runs, "Omega", ">", 1)

runs.all$p.adj <- p.adjust(runs.all$p, method="BH")
sum(runs.all$p.adj < 0.05, na.rm=T)

pdf("WaW_all.pdf", height=7, width=9)
make.plots(runs.all, title="all")
dev.off()

runs.all.mtc <- ddply(slr_all, c("stable_id"), process.runs, "lower", ">", 1)

runs.all.mtc$p.adj <- p.adjust(runs.all.mtc$p, method="BH")
sum(runs.all.mtc$p.adj < 0.05, na.rm=T)

pdf("WaW_all_genewise_MTC.pdf", height=7, width=9)
make.plots(runs.all.mtc, title="all (genewise MTC)")
dev.off()

runs.all.gmtc <- ddply(slr_all, c("stable_id"), process.runs, "Adj.Pval", "<", 0.05)

runs.all.gmtc$p.adj <- p.adjust(runs.all.gmtc$p, method="BH")
sum(runs.all.gmtc$p.adj < 0.05, na.rm=T)

pdf("WaW_all_genomewise_MTC.pdf", height=7, width=9)
make.plots(runs.all.gmtc, title="all (genewise MTC)")
dev.off()


# Example plots
runs.all.subset <- subset(runs.all.gmtc, nsites >= 5)
order(runs.all.subset$p)[1:10]
runs.all.subset[order(runs.all.subset$p)[1:10], ]

pdf("example_plots.pdf", height=7, width=9)
runs.test(as.numeric(subset(slr_all, stable_id == "ENSP00000016913")$Adj.Pval < 0.05), "left.sided", threshold=0.5, plot=T)
runs.test(as.numeric(subset(slr_all, stable_id == "ENSP00000350418")$Adj.Pval < 0.05), "left.sided", threshold=0.5, plot=T)
runs.test(as.numeric(subset(slr_all, stable_id == "ENSP00000308279")$Adj.Pval < 0.05), "left.sided", threshold=0.5, plot=T)
runs.test(as.numeric(subset(slr_all, stable_id == "ENSP00000437563")$Adj.Pval < 0.05), "left.sided", threshold=0.5, plot=T)
runs.test(as.numeric(subset(slr_all, stable_id == "ENSP00000438038")$Adj.Pval < 0.05), "left.sided", threshold=0.5, plot=T)
dev.off()

process.runs(subset(slr_all, stable_id == "ENSP00000480132"))

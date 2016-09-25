library(plyr)

slr_yeast <- read.table("data/yeast_all.tab", sep="\t", header=T, stringsAsFactors=F)
nrow(slr_yeast)
slr_yeast <- subset(slr_yeast, !grepl("Single", Note))
nrow(slr_yeast)

hist(subset(slr_yeast, Omega > 1)$Pval, breaks=100)

slr_sane <- subset(slr_yeast, Omega < 5)
hist(slr_sane, breaks=100)

pdf("yeast_omega_dist.pdf", height=7, width=9)
ddply(slr_sane, "stable_id", function(df) {
  hist(df$Omega, breaks=100, main=df$stable_id[1], xlab=expression(omega))
})
dev.off()

sum(slr_yeast$Omega > 1)

yeast_pos <- subset(slr_yeast, Omega > 1 & Adj.Pval < 0.05)
subset(yeast_pos, Pval < 0.05)

hist(yeast_pos$Pval, breaks=100, )
hist(yeast_pos$Omega, breaks=100)

p.adj <- p.adjust(yeast_pos$Pval, method="BH")
sum(p.adj < 0.1)

cn <- c("Neutral", "Optimal", "omega", "lower", "upper", "LRT_Stat", "Pval", "Adj.Pval", "Q-value", "Result", "Note")
slr_results <- read.fwf("~/Downloads/paml_comparison/YDR233C_matched.res", widths=c(9, 8, 9, 9, 9, 7, 9, 10, 11, 11, 7, 12),
                        header=F, row.names=1, stringsAsFactors=F)
colnames(slr_results) <- cn
slr_results$omega <- as.numeric(substr(slr_results$omega, 1, 6))
slr_results$lower <- as.numeric(substr(slr_results$lower, 1, 6))
slr_results$upper <- as.numeric(substr(slr_results$upper, 1, 6))


# slr_results <- subset(slr_results, Note != " Single char")

hist(slr_results$omega, breaks=100)

tracks.ylim <- c(0, 2)
scores <- rep(0.5, nrow(slr_results))
scores[slr_results$lower > 1] <- 1
scores[slr_results$upper < 1] <- 0
slr_track <- data.frame(pos=1:nrow(slr_results),
                        score=scores,
                        y_hi=pmin(slr_results$upper, tracks.ylim[2])/2,
                        y_lo=pmin(slr_results$lower, tracks.ylim[2])/2,
                        id="SLR",
                        layout="below",
                        color.gradient="blue,gray,red",
                        height=6)

sim <- PhyloSim()
readAlignment(sim, "~/Downloads/paml_comparison/YDR233C_prank.best.fas")

readTree(sim, "~/Downloads/paml_comparison/RAxML_bestTree.YDR233C")

aln <- sim$.alignment
aln <- aln.tx(aln)
# aln <- aln[, slr_results$SiteNumber]

# if (!is.null(trim)) {
#   aln <- aln[, trim[1]:trim[2]]
#   slr_track <- slr_track[trim[1]:trim[2], ]
# }

axis.text.size <- 4

sim$.alignment <- aln 
aln.len <- getAlignmentLength(sim)

plot.PhyloSim(sim,
              plot.chars=F,
              plot.labels=T,
              # aln.char.text.size=1.8,
              # num.pages=num.pages,
              layout.ancestors=F,
              aln.to.tree.size=5,
              axis.text.size=axis.text.size,
              tracks=list(SLR=slr_track),
              tracks.ylim=tracks.ylim,
              # tree.xlim=xlims[[symbol]],
              plot.legend=T,
              color.scheme='protein'
)

subset(slr_results, omega > 1)

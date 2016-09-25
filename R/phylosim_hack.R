# Need to hack phylosim to add arrows around SNPs?
# 
## 
# Add an annotation type to tracks?
# - Arrows can be specified as label + position
# - Domains need ranges, else they come fragmented
# Arrows can be plotted as 'points' of the type 24/25
# library(phylosim)
library(devtools)

aln.tx <- function(aln) {
  assign("PSIM_FAST", TRUE, envir=.GlobalEnv)
  # Translates a codon alignment into proteins.
  ca <- CodonAlphabet()
  pep.aln <- matrix(nrow=nrow(aln), ncol=ncol(aln)/3)
  rownames(pep.aln) <- rownames(aln)
  for (i in 1:nrow(pep.aln)) {
    codons <- paste(aln[i,],collapse='')
    for (j in 1:ncol(pep.aln)) {
      cdna.lo <- (j-1)*3+1
      cdna.hi <- (j-1)*3+3
      codon <- substr(codons, cdna.lo, cdna.hi)
      if (codon == '---') {
        pep.aln[i,j] <- '-'
      } else if (length(grep('N', codon)) > 0) {
        pep.aln[i,j] <- 'X'
      } else {
        aa <- translateCodon(ca, codon)
        if (length(aa) > 0) {
          pep.aln[i,j] <- aa
        } else {
          pep.aln[i,j] <- 'X'
        }      
      }
    }
  }
  return(pep.aln)
}

aln_file <- "~/Documents/projects/sarah/data/541_1_prank.best.fas"
tree_file <- "~/Documents/projects/sarah/data/541_1_slr.nwk"
slr_file <- "~/Documents/projects/sarah/data/ENSP00000085219_541_1_results.csv"

tracks.ylim <- c(0, 2)

slr_results <- read.table(slr_file, sep=",", row.names=1, header=T)
scores <- rep(0.5, nrow(slr_results))
scores[slr_results$SLR_Lower > 1] <- 1
scores[slr_results$SLR_Upper < 1] <- 0

slr_track <- data.frame(pos=1:nrow(slr_results),
                        score=scores,
                        y_hi=pmin(slr_results$SLR_Upper, tracks.ylim[2]),
                        y_lo=pmin(slr_results$SLR_Lower, tracks.ylim[2]),
                        id="SLR",
                        layout="above",
                        color.gradient="blue,gray,red",
                        height=6)

tracks <- list(SLR=slr_track)


sim <- PhyloSim()
readAlignment(sim, aln_file)

readTree(sim, tree_file)

aln <- sim$.alignment
aln <- aln.tx(aln)
aln <- aln[, slr_results$SiteNumber]

## if (!is.null(trim)) {
##   aln <- aln[, trim[1]:trim[2]]
##   slr_track <- slr_track[trim[1]:trim[2], ]
## }

axis.text.size <- 4
num.pages <- 3

sim$.alignment <- aln 
aln.len <- getAlignmentLength(sim)

markers <- data.frame(pos=1:aln.len, color="red", layout="both")
ranges <- data.frame(id=c("JAM1", "Jnk1"), from=c(1, 400), to=c(300, 650),
                     color=c("red", "green"))
## ranges <- data.frame(id=c("JAM1"), from=c(1), to=c(200),
##                      color=c("red"))


load_all("~/Documents/projects/phylosim")
plot.PhyloSim(sim,
              plot.chars=F,
              plot.labels=T,
              # aln.char.text.size=1.8,
              # num.pages=num.pages,
              layout.ancestors=F,
              aln.to.tree.size=5,
              axis.text.size=axis.text.size,
              # tracks=tracks,
              # tracks.ylim=tracks.ylim,
              # range.tracks=ranges,
              # markers=markers,
              # tree.xlim=xlims[[symbol]],
              plot.legend=T,
              color.scheme='protein'
)
 


#
## Misc tests
# 
gg_test <- data.frame(x=c(rep(c(1,2,3), 2)), y=c(rep(1, 3), rep(2,3)))
gg_test$xmin <- gg_test$x - 0.5
gg_test$xmax <- gg_test$x + 0.5
gg_test$ymin <- gg_test$y - 0.5
gg_test$ymax <- gg_test$y + 0.5

gg_test2 <- data.frame(x=c(rep(seq(1,3), 2)), y=c(rep(3, 6)))

plot_test <- ggplot() +
  geom_rect(data=gg_test, aes(x=x, y=y, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) + 
  geom_point(aes(x=x, y=y), data=gg_test2, shape=25, fill=1, size=18) +
  theme_bw()
  
print(plot_test)


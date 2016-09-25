# PDB master table stuff

aln.plots <- function(aln_file, tree_file, slr_file, trim=NA, num.pages=3, tracks.ylim=c(0, 1)) {
  slr_results <- read.table(slr_file, sep="\t", row.names=1, header=T, comment.char="")
  scores <- slr_results$Omega
  # scores <- rep(0.5, nrow(slr_results))
  # scores[slr_results$lower > 1] <- 1
  # scores[slr_results$upper < 1] <- 0

  slr_track <- data.frame(pos=1:nrow(slr_results),
                          score=scores/max(scores),
                          y_hi=scores/max(scores),
                          y_lo=0,
                          id="SLR",
                          layout="below",
                          color.gradient="white,red",
                          height=6)
  
  tracks <- list(SLR=slr_track)
  sim <- PhyloSim()
  readAlignment(sim, aln_file)
  
  if(!is.na(tree_file)) {
    readTree(sim, tree_file)
  }
  
  aln <- sim$.alignment
  aln <- aln.tx(aln)
  aln <- aln[, as.numeric(rownames(slr_results))]

#   if (!any(is.na(trim))) {
#     aln <- aln[, trim[1]:trim[2]]
#     slr_track <- slr_track[trim[1]:trim[2], ]
#   }
  
  axis.text.size <- 4
  
  sim$.alignment <- aln 
  aln.len <- getAlignmentLength(sim)
  
  plot.PhyloSim(sim,
                plot.chars=F,
                plot.labels=T,
                # aln.char.text.size=1.8,
                num.pages=num.pages,
                layout.ancestors=F,
                aln.to.tree.size=5,
                axis.text.size=axis.text.size,
                tracks=tracks,
                tracks.ylim=tracks.ylim,
                # tree.xlim=xlims[[symbol]],
                plot.legend=T,
                color.scheme='protein'
  )
}

mean_pdb_omega <- ddply(pdb_master, "stable_id", function(df) {
  data.frame(omega=mean(df$omega))
})

tail(mean_pdb_omega[order(mean_pdb_omega$omega), "stable_id"], n=30)

pdf("example_aln.pdf")
aln.plots(aln_file="4085_1_prank.best.fas",
          tree_file="4085_1_slr.nwk",
          slr_file="ENSP00000294008_4085_1_matched.res",
          # trim=c(1, 100), # Subset of the alignment to plot
          # This controls the range for the tracks. I just hacked it together, probably doesn't work in every case
          # In particular, I wouldn't set the lower bound to anything but 0 ;)
          tracks.ylim=c(0, 1.05),
          num.pages=1)
dev.off()
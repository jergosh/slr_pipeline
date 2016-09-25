#
# 
#
library(phylosim)
library(stringr)
setwd("~/Documents/projects/slr_pipeline")

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

aln.plots <- function(df, trim=NA, num.pages=3, tracks.ylim=c(0, 1)) {
  ens <- as.character(df$ID[1])
  print(ens)
  aln_file <- paste(c("data/fubar/Eutheria/", substr(ens, 1, 2), "/", ens, ".fa"), collapse="")
  tree_file <- paste(c("data/fubar/Eutheria/", substr(ens, 1, 2), "/", ens, ".nwk"), collapse="")
  print(aln_file)
  print(tree_file)

  # Scale 
  alpha <- df$alpha
  beta <- df$beta
  max_rate <- max(alpha, beta)
  alpha <- alpha / max_rate
  beta <- beta / max_rate
  
  alpha_track <- data.frame(pos=1:nrow(df),
                            y_hi=alpha,
                            y_lo=0,
                            id="dS",
                            layout="below",
                            color.gradient="blue,gray,red",
                            height=6)
  
  beta_track <- data.frame(pos=1:nrow(df),
                            y_hi=beta,
                            y_lo=0,
                            id="dS",
                            layout="below",
                            color.gradient="blue,gray,red",
                            height=6)
  
  
  tracks <- list(alpha=alpha_track, beta=beta_track)
  sim <- PhyloSim()
  readAlignment(sim, aln_file)
  
  if(!is.na(tree_file)) {
    readTree(sim, tree_file)
  }
  
  aln <- sim$.alignment
  aln <- aln.tx(aln)
  aln <- aln[, df$pos]
  
  if (!any(is.na(trim))) {
    aln <- aln[, trim[1]:trim[2]]
    alpha_track <- alpha_track[trim[1]:trim[2], ]
    beta_track <- beta_track[trim[1]:trim[2], ]
  }
  
  axis.text.size <- 4
  
  sim$.alignment <- aln 
  aln.len <- getAlignmentLength(sim)
  
  plot.PhyloSim(sim,
                plot.chars=F,
                plot.labels=T,
                # aln.char.text.size=1.8,
                num.pages=num.pages,
                layout.ancestors=F,
                aln.to.tree.size=7,
                axis.text.size=axis.text.size,
                tracks=tracks,
                tracks.ylim=tracks.ylim,
                # tree.xlim=xlims[[symbol]],
                plot.legend=T,
                color.scheme='protein'
  )
}

aln.plots.apply <- function(df) {
  # pdf()
  aln.plots(df, num.pages="auto")
}

daply(fubar_table[1:1000, ], .variables=c("ID"), .fun=aln.plots.apply)

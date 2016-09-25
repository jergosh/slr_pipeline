library(phylosim)
library(stringr)
setwd("~/Documents/projects/slr_pipeline/")

aln.plots <- function(symbol, aln_file, tree_file, gene_slr, lcr_pos=NA, tx=TRUE, idx_subset=NA) {
  
  plotAln <- function(inF, tree=NA, tx=FALSE, tracks=NA, num.pages=1, idx_subset=NA) {
    sim <- PhyloSim()
    readAlignment(sim, inF)
    
    if(!is.na(tree)) {
      readTree(sim, tree)
    }
    
    aln <- sim$.alignment
    if (tx) {
      aln <- aln.tx(aln)
    }
    
    if (!is.na(idx_subset)) {
      aln <- aln[, idx_subset]
    }
    
    sim$.alignment <- aln
    aln.len <- getAlignmentLength(sim)
    print(aln.len)
    extra <- scale_x_continuous(breaks=seq(from=1, to=100, by=30))
    
    if (num.pages > 1) {
      axis.text.size = 10
    } else {
      axis.text.size = 7
    }
    
    plot.PhyloSim(sim,
                  plot.chars=F,
                  plot.labels=T,
                  # aln.char.text.size=1.8,
                  num.pages=num.pages,
                  layout.ancestors=F,
                  aln.to.tree.size=5,
                  aln.extras=extra,
                  axis.text.size=axis.text.size,
                  tracks=tracks,
                  tree.xlim=c(0, 0.25),
                  plot.legend=F,
                  plot.tree=T,
                  color.scheme='protein'
    )
  }
  
  low <- "#0000FF"
  med <- "#888888"
  hi <- "#FF0000"
  scores = rep(0.5, nrow(gene_slr))
  scores[gene_slr$Omega > 1 & gene_slr$Adj.Pval < thr] = 1
  scores[gene_slr$Omega < 1 & gene_slr$Adj.Pval < thr] = 0
  
  slr_track = data.frame(pos=gene_slr$Site, score=scores, y_hi=pmin(gene_slr$Omega, 0.95), id="SLR", layout="below", color.gradient="blue,gray,red")
  tracks=list(SLR=slr_track)
  
  if (!is.na(lcr_pos)) {
    lcr_track = data.frame(pos=lcr_pos, scores=0.6, y_hi=0.6, id="Low complexity", layout="above", color="orange")
    tracks$LCR = lcr_track
  }
  
  
  # pdf(str_c("figures/", symbol, ".pdf"), height=27*0.8, width=29.7*0.8)
  # pdf(str_c(symbol, ".pdf"), height=27*0.8, width=29.7*0.8)
  pdf(str_c(symbol, ".pdf"), height=8, width=29.7*0.8)
  plotAln(aln_file, tx=tx, tracks=tracks, tree=tree_file, idx_subset=idx_subset) # , 
  dev.off()
}

aln_file = "data/20320_3_matched.fa"
tree_file = "data/20320_3_matched.nh"
gene_slr = subset(slr_all, stable_id == "ENSP00000379823")
gene_lcr = subset(feature_data, type=="seg" & stable_id=="ENSP00000379823")
lcr_pos = ddply(gene_lcr, c("start", "end"), function(x) {data.frame(pos=seq(x[1, "start"], x[1, "end"]))})$pos
lcr_pos_aln = subset(slr_all, stable_id == "ENSP00000379823")$Site[lcr_pos]

aln.plots("ENSP00000379823", aln_file, tree_file, gene_slr, lcr_pos_aln)

aln_file = "data/20322_2_matched.fa"
tree_file = "data/20322_2_matched.nh"
gene_slr = subset(slr_all, stable_id == "ENSP00000379866")
gene_lcr = subset(feature_data, type=="seg" & stable_id=="ENSP00000379866")
lcr_pos = ddply(gene_lcr, c("start", "end"), function(x) {data.frame(pos=seq(x[1, "start"], x[1, "end"]))})$pos
lcr_pos_aln = subset(slr_all, stable_id == "ENSP00000379866")$Site[lcr_pos]

aln.plots("ENSP00000379866", aln_file, tree_file, gene_slr, lcr_pos_aln)

library(ggplot2)
library(seqinr)

setwd("~/Documents/projects/slr_pipeline")

# The dataset is 101_7
symbol = "TRIM21"
human_id = "ENSP00000254436"
res_file = "101_7_matched.res"
aln_file = "101_7_prank.best.fas"
tree_file = "101_7.nh"

plot.SLR(aln.file=aln_file,
         tree.file=tree_file,
         res.file=res_file,
         ref.id=human_id,
         trim=T,
         slr.ylim=c(0, 5),
         main="TRIM21",
         id.translate=F,
         num.pages="auto")

sim <- PhyloSim()
readAlignment(sim, aln_file)
readTree(sim, tree_file)
tree <- sim$.phylo
tree$edge.length[is.na(tree$edge.length)] <- 0.01
tree$edge.length <- rep(1, length(tree$edge.length))
sim$.phylo <- tree

aln <- sim$.alignment
ref.idx <- aln[human_id, ] != "-"
aln <- aln[, ref.idx]
write.dna(aln, file=)

plot.PhyloSim(sim,
              plot.chars=F,
              plot.labels=T,
              aln.char.text.size=1.8,
              num.pages="auto",
              layout.ancestors=F,
              aln.to.tree.size=5,
              #aln.extras=extra,
              axis.text.size=10,
              #tracks=list(SLR=slr_track),
              #tracks.ylim=slr.ylim,
              #track.height=7,
              #plot.legend=T,
              color.scheme='protein'
) 


# pfam.ann  <- data.frame(type="Pfam",
#                         id=c("B-box-type zinc finger", "SPRY domain"),
#                         from=c(91, 356),
#                         to=c(128, 489),
#                         colour=domain.pal[1:2],
#                         label=TRUE
# )


pfam.ann <- data.frame(type="Pfam",
                       id=c("B-box-type zinc finger", "SPRY-associated", "SPRY domain"),
                       from=c(88, 289, 339),
                       to=c(128, 337, 455),
                       colour=domain.pal[1:3],
                       label=TRUE
)


slr.tmp <- subset(slr_all, stable_id == "ENSP00000254436")[1:475, ]
# slr.tmp <- subset(slr_all, stable_id == "ENSP00000369373")[1:493, ]

# slr.pos <- slr.tmp$Pval[slr.tmp$Omega > 1]
# slr.tmp$Adj.Pval[slr.tmp$Omega > 1] <- p.adjust(slr.pos, method="BH")
slr.tmp$Pval[slr.tmp$Omega < 1] <- runif(sum(slr.tmp$Omega < 1))
hist(slr.tmp$Pval, breaks=100)
slr.tmp$Adj.Pval <- p.adjust(slr.tmp$Pval, method="BH")
subset(slr.tmp, Adj.Pval < 0.2 & Omega > 1)
write.csv(slr.tmp, file="TRIM21.csv")


y.lim <- c(0, 5)
seq <- "MASAARLTMMWEEVTCPICLDPFVEPVSIECGHSFCQECISQVGKGGGSVCPVCRQRFLLKNLRPNRQLANMVNNLKEISQEAREGTQGERCAVHGERLHLFCEKDGKALCWVCAQSRKHRDHAMVPLEEAAQEYQEKLQVALGELRRKQELAEKLEVEIAIKRADWKKTVETQKSRIHAEFVQQKNFLVEEEQRQLQELEKDEREQLRILGEKEAKLAQQSQALQELISELDRRCHSSALELLQEVIIVLERSESWNLKDLDITSPELRSVCHVPGLKKMLRTCAVHITLDPDTANPWLILSEDRRQVRLGDTQQSIPGNEERFDSYPMVLGAQHFHSGKHYWEVDVTGKEAWDLGVCRDSVRRKGHFLLSSKSGFWTIWLWNKQKYEAGTYPQTPLHLQVPPCQVGIFLDYEAGMVSFYNITDHGSLIYSFSECAFTGPLRPFFSPGFNDGGKNTAPLTLCPLNIGSQGSTDY"

slr.tmp$human_idx <- 1:nrow(slr.tmp)
  
slr.tmp$colour <- "grey"
# slr.tmp$colour[slr.tmp$lower > 1 & slr.tmp$Adj.Pval < 0.05] <- "red"
slr.tmp$colour[slr.tmp$lower > 1] <- "red"
slr.tmp$colour[slr.tmp$upper < 1] <- "blue"
  
# TODO This should be moved to plot.gene(), once y.lim is passed as an argument
# TODO These should be manipulated so that 
slr.tmp$Omega <- pmin(y.lim[2], slr.tmp$Omega)
slr.tmp$lower <- pmin(y.lim[2], slr.tmp$lower)
slr.tmp$upper <- pmin(y.lim[2], slr.tmp$upper)
  
  
slr.ann <- data.frame(type="Selective constraint",
                      pos=slr.tmp$human_idx,
                      value=slr.tmp$Omega,
                      y_hi=slr.tmp$upper,
                      y_lo=slr.tmp$lower,
                      y_min=y.lim[1],
                      y_max=y.lim[2],
                      colour=slr.tmp$colour)

plot.gene(seq,
          page.len=160,
          # range.ann=rbind(pfam.ann, disorder.ann, tm.ann),
          # marker.ann=marker.ann,
          range.ann=pfam.ann,
          track.ann=slr.ann, # rsa.ann),
          title="TRIM21")

ggsave(filename="TRIM_21.pdf",
       height=9, width=8, units="in", limitsize=F)


seq.trim7 <- "MAAVGPRTGPGTGAEALALAAELQGEATCSICLELFREPVSVECGHSFCRACIGRCWERPGAGSVGAATRAPPFPLPCPQCREPARPSQLRPNRQLAAVATLLRRFSLPAAAPGEHGSQAAAARAAAARCGQHGEPFKLYCQDDGRAICVVCDRAREHREHAVLPLDEAVQEAKELLESRLRVLKKELEDCEVFRSTEKKESKELLKQMAAEQEKVGAEFQALRAFLVEQEGRLLGRLEELSREVAQKQNENLAQLGVEITQLSKLSSQIQETAQKPDLDFLQEFKSTLSRCSNVPGPKPTTVSSEMKNKVWNVSLKTFVLKGMLKKFKEDLRGELEKEEKVELTLDPDTANPRLILSLDLKGVRLGERAQDLPNHPCRFDTNTRVLASCGFSSGRHHWEVEVGSKDGWAFGVARESVRRKGLTPFTPEEGVWALQLNGGQYWAVTSPERSPLSCGHLSRVRVALDLEVGAVSFYAVEDMRHLYTFRVNFQERVFPLFSVCSTGTYLRIWP"
trim7.ann <- data.frame(type="Pfam",
                        id=c("B-box-type zinc finger", "SPRY-associated", "SPRY domain"),
                        from=c(88, 289, 339),
                        to=c(128, 337, 455),
                        colour=domain.pal[4:6],
                        label=TRUE)

slr.trim7 <- subset(slr_all, stable_id == "ENSP00000274773")

slr.trim7$colour <- "grey"
slr.trim7$colour[slr.trim7$lower > 1] <- "red"
slr.trim7$colour[slr.trim7$upper < 1] <- "blue"

trim7.track <- data.frame(type="Selective constraint",
                      pos=slr.trim7$human_idx,
                      value=slr.trim7$Omega,
                      y_hi=slr.trim7$upper,
                      y_lo=slr.trim7$lower,
                      y_min=y.lim[1],
                      y_max=y.lim[2],
                      colour=slr.trim7$colour)

plot.gene(seq.trim7,
          page.len=171,
          range.ann=trim7.ann,
          # marker.ann=marker.ann,
          track.ann=trim7.track, # rsa.ann),
          title="TRIM7")

ggsave(filename="TRIM_7.pdf",
       height=9, width=8, units="in", limitsize=F)



# disorder.ann <- data.frame(type="Disorder",
#                            id="",
#                            from=start(ranges(disorder.tmp)),
#                            to=end(ranges(disorder.tmp)),
#                            colour="#BBBBBB",
#                            label=FALSE
# )


slr.tmp.ighg1 <- subset(slr_all, stable_id == "ENSP00000481725")
slr.tmp.ighg1$colour <- "grey"
# slr.tmp$colour[slr.tmp$lower > 1 & slr.tmp$Adj.Pval < 0.05] <- "red"
slr.tmp.ighg1$colour[slr.tmp.ighg1$lower > 1] <- "red"
slr.tmp.ighg1$colour[slr.tmp.ighg1$upper < 1] <- "blue"

slr.tmp.ighg1$Pval[slr.tmp.ighg1$Omega < 1] <- runif(sum(slr.tmp.ighg1$Omega < 1))
slr.tmp.ighg1$Adj.Pval <- p.adjust(slr.tmp.ighg1$Pval, method="BH")
write.csv(slr.tmp.ighg1, file="IGHG1.csv")



ighg1.track <- data.frame(type="Selective constraint",
                          pos=slr.tmp.ighg1$human_idx,
                          value=pmin(slr.tmp.ighg1$Omega, y.lim[2]),
                          y_hi=pmin(slr.tmp.ighg1$upper, y.lim[2]),
                          y_lo=pmax(slr.tmp.ighg1$lower, y.lim[1]),
                          y_min=y.lim[1],
                          y_max=y.lim[2],
                          colour=slr.tmp.ighg1$colour)

seq.ighg1 <- "MELGLCWVFLVAILEGVQCEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYEMNWVRQAPGKGLEWVSYISSSGSTIYYADSVKGFTVTRDNAKNSLYLQMNSLRSEDTAVYYCARQCENPHPESVRNPRYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
ighg1.ann <- data.frame(type="Pfam",
                        id=c(rep("Immunoglobulin C1-set", 3), "Immunoglobulin V-set domain"),
                        from=c(160, 275, 381, 26), # c(26, 160, 275, 283, 381, 384), 
                        to=c(237, 361, 463, 143), # c(143, 237, 361, 357, 463, 461), 
                        colour=c(rep(domain.pal[5], 3), domain.pal[6]),
                        label=TRUE)

plot.gene(seq.ighg1,
          page.len=159,
          range.ann=ighg1.ann,
          track.ann=ighg1.track,
          title="IGHG1")

ggsave(filename="IGHG1.pdf",
       height=9, width=8, units="in", limitsize=F)


## 
## Str. mapping
## 
pdb.colour <- data.frame(pdb_id="2iwg", pdb_chain="B", pdb_pos=3:181, color=slr.tmp$colour[287:465], group="TRIM21")
write.table(pdb.colour, file="TRIM21_colour.txt", sep="\t", quote=F, row.names=F)
  2:182
3:181 287:465

# IGHG1
265:471
237:443
# slr.tmp.ighg1 <- subset(slr_all, stable_id == "ENSP00000481725" & human_idx >= 265 & human_idx <= 471)

pdb.colour.2 <- data.frame(pdb_id="2iwg", pdb_chain="A", pdb_pos=237:443, color=subset(slr.tmp.ighg1, human_idx >= 265 & human_idx <= 471)$colour, group="IGHG1")
write.table(rbind(pdb.colour, pdb.colour.2), file="TRIM21_colour.txt", sep="\t", quote=F, row.names=F)
slr.tmp.ighg1[46, ]


## TRIM21 is 101_7
## IGHG1 is 1406_2
## 

trim21_aa <- "~/Documents/projects/selection/trim21/data/TRIM21_aln.fa"
aln <- read.fasta(file=trim21_aa, seqtype="AA", forceDNAtolower=F)
ref.idx <- aln[[human_id]] != "-"
aln <- lapply(aln, function(x) {x[ref.idx]})
write.fasta(aln, file.out="TRIM21_aln_trimmed.fa", names=names(aln))

ighg1_aa <- "~/Documents/projects/selection/trim21/data/IGHG1_aln.fa"
aln <- read.fasta(file=ighg1_aa, seqtype="AA", forceDNAtolower=F)
ref.idx <- aln[["ENSP00000481725"]] != "-"
aln <- lapply(aln, function(x) {x[ref.idx]})
write.fasta(aln, file.out="IGHG1_aln_trimmed.fa", names=names(aln))


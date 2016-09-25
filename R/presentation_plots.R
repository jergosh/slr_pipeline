res_file <- "251_5_matched.res"
tree_file <- "251_5.nh"
aln_file <- "251_5_prank.best.fas"

plot.SLR(aln_file,
         tree_file,
         res_file,
         ref.id="ENSP00000369814",
         trim=FALSE, 
         main=NA,
         id.translate=FALSE,
         slr.ylim=c(0, 4),
         mult.test=TRUE,
         sign.cutoff = 0.1,
         num.pages=1
)

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

library(RColorBrewer)
head(pdb_master_globular)
setwd("~/Documents/projects/slr_pipeline")
results <- "results/2015-04-08"
setwd(results)

pdf("Histograms_sec_buried.pdf", height=10, width=20)
par(mfrow=c(2, 4))
ddply(pdb_master_globular, c("sec_simple", "buried"), function(df) {
  
  hist(df$omega, main=paste(df$sec_simple[1], df$buried[1]),
       xlim=c(0, 5), ylim=c(0, 20), breaks=seq(0, 100, 0.05), freq=FALSE,
       col=c(rep(rgb(0, 0, 1, 0.4), 20), "white", rep(rgb(1, 0, 0, 0.4), 80)),
       xlab="Omega", cex.lab=1.5, cex.main=2)

  TRUE
})
dev.off()

pdf("Histograms_sec_buried_domain.pdf", height=20, width=24)
par(mfrow=c(4, 4))
ddply(pdb_master_globular, c("sec_simple", "buried", "is_domain"), function(df) {
  hist(df$omega, main=paste(df$sec_simple[1], df$buried[1], df$is_domain[1]),
       xlim=c(0, 2), ylim=c(0, 35), breaks=seq(0, 100, 0.025), freq=FALSE)
  nrow(df)
})
dev.off()

par(mfrow=c(1, 1))
hists.data <- dlply(pdb_master_globular, c("sec_simple", "buried"), function(df) {
  hist(df$omega, main=paste(df$sec_simple[1], df$buried[1], df$is_domain[1]),
       xlim=c(0, 2), ylim=c(0, 20), breaks=seq(0, 100, 0.025), freq=FALSE)
})

# FIXME: these could also be melted for use with ggplot2
hist.labels <- attr(hists.data, "split_labels")
hist.pal <- brewer.pal(4, "Set2")
i <- 1
plot(NA, xlim=c(0, 2), ylim=c(0, 2))
for (it in hists.data) {
  # Dashed and solid for buriedness
  # Markers for domain and outside
  # Colours for different kinds of structure
  if (hist.labels$buried[i]) {
    ltype <- 1
  } else {
    ltype <- 2
  }
#   if (hist.labels$is_domain[i]) {
#     symbol <- 0    
#   } else {
#     symbol <- 2    
#   }
  lcolor <- hist.pal[hist.labels$sec_simple[i]]
  lines(it$mids, it$density, type="l", lty=ltype, col=lcolor)
  i <- i + 1
}

# TODO Add the means
means.data <- ddply(pdb_master_globular, c("sec_simple", "buried"), function(df) {
  data.frame(mean=mean(df$omega), median=median(df$omega))
})
abline(v=means.data$mean, col=hist.pal[means.data$sec_simple], lty=3)

legend("topright", legend=levels(hist.labels$sec_simple), fill=hist.pal)

## Another idea would be to use a violin plot with the mean and SD marked 
## (like Karin and Anders did in their paper)
library(vioplot)
library(lattice)

ggplot(pdb_master_globular, aes(y=omega, x=sec_simple)) +
  coord_cartesian(ylim=c(0, 1)) + 
  geom_violin()

vioplot(subset(pdb_master_globular, sec_simple == "Helix" & buried == TRUE)$omega,
        subset(pdb_master_globular, sec_simple == "Beta sheet" & buried == TRUE)$omega,
        subset(pdb_master_globular, sec_simple == "Coil" & buried == TRUE)$omega,
        subset(pdb_master_globular, sec_simple == "Turn" & buried == TRUE)$omega,
        subset(pdb_master_globular, sec_simple == "Helix" & buried == FALSE)$omega,
        subset(pdb_master_globular, sec_simple == "Beta sheet" & buried == FALSE)$omega,
        subset(pdb_master_globular, sec_simple == "Coil" & buried == FALSE)$omega,
        subset(pdb_master_globular, sec_simple == "Turn" & buried == FALSE)$omega,
        names=c("Helix", "Beta sheet", "Coil", "Turn", "Helix", "Beta sheet", "Coil", "Turn"),
        h=0.5, ylim=c(0, 5))

ddply(pdb_master_globular, "sec_simple", function(df) {mean(df$omega, na.rm=T)})
ddply(pdb_master_globular, "sec_simple", function(df) {nrow(df)/nrow(pdb_master_globular) } )



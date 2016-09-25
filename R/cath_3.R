library(hexbin)
library(lattice)

dists <- read.table("data/cath_dists_ml.tab", sep="\t",stringsAsFactors=F, header=T)

hist(dists$dist, breaks=10000, xlim=c(0, 50))

dists_2 <- dists
tmp <- dists_2$id_1
dists_2$id_1 <- dists_2$id_2
dists_2$id_2 <- tmp
dists <- rbind(dists, dists_2)
str(dists)


cath_cors_dists <- join(cath.cors, dists)
# cath_cors_dists <- subset(cath_cors_dists, cath_id == "2.60.40.10" & complete.cases(cath_cors_dists))
cath_cors_dists <- subset(cath_cors_dists, n_obs > 50 & dist < 10)
cath_cors_dists <- cath_cors_dists[order(cath_cors_dists$dist), ]

nrow(cath_cors_dists)
plot(cath_cors_dists$dist, cath_cors_dists$cor) #, xlim=c(-2, 1))

hexbin_dists <- hexbin(cath_cors_dists$dist, cath_cors_dists$cor, ybnds=c(-0.505, 0.8))
hexbin_gpcr <- hexbin(cath_cors_dist_gpcr$dist, cath_cors_dist_gpcr$cor)

hexbinplot(cath_cors_dists$cor~cath_cors_dists$dist)

pdf("cors_all_hexbin.pdf", height=7, width=8.5)
hvp <- plot(hexbin_dists, xlab="Distance", ylab="Correlation")
hexVP.loess(hbin=hexbin_dists, hvp=hvp$plot.vp, span=0.4, n=20000, col="black")
# hexVP.loess(hbin=hexbin_gpcr, hvp=hvp$plot.vp, span=0.4, n=20000)
dev.off()

pdf("cors_gpcr_hexbin.pdf", height=7, width=8.5)
hvp <- plot(hexbin_gpcr, xlab="Distance", ylab="Correlation")
hexVP.loess(hbin=hexbin_dists, hvp=hvp$plot.vp, span=0.4, n=20000, col="black")
hexVP.loess(hbin=hexbin_gpcr, hvp=hvp$plot.vp, span=0.4, n=20000)
dev.off()

plot(cath_cors_dists$dist, cath_cors_dists$cor)

cath_cors_dist_subset <- subset(cath_cors_dists, cath_id != "1.20.1070.10" & dist < 10)
cath_cors_dist_gpcr <- subset(cath_cors_dists, cath_id == "1.20.1070.10" & dist < 10)
plot(cath_cors_dist_subset$dist, cath_cors_dist_subset$cor) #, xlim=c(0, 50))
points(cath_cors_dist_gpcr$dist, cath_cors_dist_gpcr$cor, col="red")

cath_cors_dist_50 <- subset(cath_cors_dist_subset, dist < 50)
plot(hexbin(cath_cors_dist_50$dist, cath_cors_dist_50$cor))
layer(panel.xyplot(x=cath_cors_dists$dist, y=pr.loess, col="blue"))

hexbinplot(cath_cors_dists$cor~cath_cors_dists$dist, main = "",
           xlab="wind",ylab="other wind", style="colorscale", type=c("r"), col.line = "red", lwd="3")

cor.loess <- loess(cath_cors_dists$cor~cath_cors_dists$dist, span=0.1)
pr.loess <- predict(cor.loess)
lines(pr.loess~cath_cors_dists$dist, col="blue")

# Do the per-family plot
plot(NA, xlim=c(0, 50), ylim=c(-1, 1))
ddply(cath_cors_dists, "cath_id", function(df) {
  if (nrow(df) > 100) {
    df <- df[order(df$dist), ]
    print(c(df$cath_id[1], nrow(df)))
    cor.loess <- loess(df$cor~df$dist, span=0.1)
    pr.loess <- predict(cor.loess)
    plot(pr.loess~df$dist, xlim=c(0, 5), ylim=c(-1, 1), main=df$cath_id[1])
    # lines(pr.loess~df$dist, col="blue")
  }  
})


plot(cath_cors_dist_gpcr$dist, cath_cors_dist_gpcr$cor, col="red", xlim=c(0, 55))
cor.loess <- loess(cath_cors_dist_gpcr$cor~cath_cors_dist_gpcr$dist, span=0.02)
pr.loess <- predict(cor.loess)
lines(pr.loess~cath_cors_dist_gpcr$dist, col="blue")

## Laundry list of things to correlate on
## SSG
## Mean omega

cath_cors_dist_smell <- subset(cath_cors_dist_gpcr, (id_1 %in% annotated) & (id_2 %in% annotated))
plot(cath_cors_dist_gpcr$dist, cath_cors_dist_gpcr$cor, xlim=c(0, 55))
points(cath_cors_dist_smell$dist, cath_cors_dist_smell$cor, col="green")

cor.loess <- loess(cath_cors_dist_smell$cor~cath_cors_dist_smell$dist, span=0.1)
pr.loess <- predict(cor.loess)
lines(pr.loess~cath_cors_dist_smell$dist, col="blue")


# The distances in the first cloud are already substantial!

# Large distance
large_dists <- subset(cath_cors_dists, dist > 40 & cath_id == "1.20.1070.10")
nrow(large_dists)
head(large_dists)

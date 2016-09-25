library(plyr)

rng <- c(0, 5)
step <- 0.1
s <- seq(rng[1], rng[2], step)

pdf("Gamma_fits.pdf")
plot(NA, xlim=rng, ylim=c(0, 1.5))
ddply(slr_all, "stable_id", function(df) {
  g <- fitdistr(df$Omega+0.001, "gamma")
  lines(s, dgamma(s, shape=g$estimate["shape"], rate=g$estimate["rate"]))
})
dev.off()

hist(pdb_master_globular$omega, breaks=100)

pdf("Gamma_fits_globular.pdf")
plot(NA, xlim=rng, ylim=c(0, 1.5))
ddply(pdb_master_globular, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
  g <- fitdistr(df$omega+0.001, "gamma")
  lines(s, dgamma(s, shape=g$estimate["shape"], rate=g$estimate["rate"]))
})
dev.off()


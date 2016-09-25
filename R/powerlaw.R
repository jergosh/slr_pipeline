library("poweRlaw")

# pdf("Pos_sel_hist.pdf")
png("pos_sel_hist.png", width=500, height=400)
par(mar=c(5.1,5.1,4.1,2.1))
hist(pos_per_gene, breaks=70, xlim=c(0, 25), freq=F, xlab="Number of sites",
     main=expression(paste("Sites with ", omega, " > 1 per gene")), cex.main=1.5, cex.lab=1.5)
# h <- hist(pos_per_gene, breaks=50)
# plot(h$counts, log="xy")
dev.off()

sel_dist <- pos_per_gene
sel_dist[sel_dist == 0] = 0.01

m_m = displ$new(pos_per_gene+1)
estimate_pars(m_m)
plot(m_m)
bs = bootstrap(m_m, no_of_sims = 1000, threads = 4)

source("http://tuvalu.santafe.edu/~aaronc/powerlaws/plfit.r")
plfit(sel_dist)

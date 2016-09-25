events <- list(100:120)
bounds <- list(c(1,300))
gs <- graphscan_1d(data=events, normalisation_factor=bounds, n_simulation=10000)
res <- cluster(gs, memory_size = 10000)
summary(res)
res@cluster

events <- list(c(1:10))
bounds <- list(c(1,300))
gs <- graphscan_1d(data=events, normalisation_factor=bounds, n_simulation=10000)
res <- cluster(gs, memory_size = 10000)
summary(res)
res@cluster

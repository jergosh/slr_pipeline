multihist <- function(dists, breaks=NULL, cols=NULL, ...) {
  hists <- list()
  hist_all <- hist(unlist(dists), breaks=breaks, plot=FALSE)
  n_breaks <- length(hist_all$mids)
  print(length(hist_all$breaks))
  
  max_count <- 0
  for (dist in dists) {
    hists[[length(hists)+1]] <- hist(dist, plot=FALSE, breaks=hist_all$breaks, ...)
    if (max(hists[[length(hists)]]$density) > max_count)
      max_count <- max(hists[[length(hists)]]$density)
  }
  
  ytop <- matrix(NA, nrow=length(hists), ncol=n_breaks)
  ybottom <- matrix(NA, nrow=length(hists), ncol=n_breaks)
  for (i in 1:length(hists)) {
    print(hists[[i]]$density)
    ytop[i, ] <- hists[[i]]$density
  }
  ybottom[1, ] <- 0.0
  for (i in 1:n_breaks) {
    dist_order <- order(ytop[, i])
    ybottom[2:length(hists), i] <- ytop[dist_order[1:(length(hists)-1)], i]
  }
  print(ybottom)
  print(ytop)

  xleft <- hist_all$breaks[1:n_breaks]
  xright <- hist_all$breaks[2:(n_breaks+1)]
  
  # possibly better to 
  plot(NA, xlim=c(min(hist_all$breaks), max(hist_all$breaks)), ylim=c(0, max_count))
  for (j in 1:n_breaks) {
    bin_order <- order(ytop[, j])
    # print(bin_order)
    rect(xleft[j], ybottom[, j], xright[j], ytop[bin_order, j], col=cols[bin_order])      
  }
}

colours <- c(rgb(1, 0, 0, alpha=0.5), rgb(0, 1, 0, alpha=0.5), rgb(0, 0, 1, alpha=0.5), rgb(1, 1, 0, alpha=0.5))
multihist(list(rnorm(1000), rnorm(10000, 1, 5), rnorm(1000, -5, 3), rnorm(10000, 10, 9)), breaks=100, cols=colours)

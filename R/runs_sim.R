library(tseries)

seq <- c(rep(0, 100), runif(20) < 0.25, rep(0, 100))

res <- acf(seq)
Box.test(seq, type="Box-Pierce")
Box.test(seq, type="Ljung-Box")

seq <- c(rep(0, 100), sample(c(rep(TRUE, 5), rep(FALSE, 35)), 40, replace=F), rep(0, 100))
Box.test(seq, type="Box-Pierce")
Box.test(seq, type="Ljung-Box")

pdf("Pval_hists_WaW.pdf")
for (w_1 in c(10, 20, 30, 40, 50)) {
  for (w_0 in c(50, 200)) {
    for (n_0 in c(0, 10)) {
      for (n in 6:10) {
        ps <- unlist(lapply(1:1000, function(x) {
          seq_0 <- sample(c(rep(1, n_0), rep(0, w_0)))
          seq <- sample(c(rep(1, n), rep(0, w_1)))
          seq_0_2 <- sample(c(rep(1, n_0), rep(0, w_0)))
          
          # Box.test(c(seq_0, seq, seq_0_2), type="Box-Pierce")$p.value
          runs.test(as.factor(c(seq_0, seq, seq_0_2)))$p.value
        }))
        
        hist(ps, breaks=100, main=paste("w_1", w_1, "w_0", w_0, "n", n, "n_0", n_0))
      }
    }
  }
}
dev.off()


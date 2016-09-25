setwd("~/Documents/phd_seminar_2014/")

t <- 0:100  # time
sig2 <- 0.01
nsim = 50
X <- matrix(rnorm(mean = -0.0005, n = nsim * (length(t) - 1), sd = sqrt(sig2/4)), 
            nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-1, 3), type = "l")
# apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

png("allele_freqs_neg.png", width=600, height=300)
plot(NA, xlab = "Time (generations)", ylab = "Allele frequency", xlim=c(0, 100), ylim=c(0, 1), type="l",
     col="darkgreen", lwd=2)

for (i in seq(1, nsim)) {
  if (length(which(X[i, ] > 1))) {
    thr = which(X[i, ] > 1)[1]
    X[i, thr:ncol(X)] = 1
  }
  
  if (length(which(X[i, ] < 0))) {
    thr = which(X[i, ] < 0)[1]
    X[i, thr:ncol(X)] = 0
  }
  
  lines(t, X[i, ], col="red", lwd=2)
}
dev.off()

plot(t, X[1, ], xlab = "Time (generations)", ylab = "Allele frequency", ylim=c(0, 1), type="l",
     col="darkgreen", lwd=2)


Y <- matrix(rnorm(mean = -0.001, n = nsim * (length(t) - 1), sd = sqrt(sig2/4)), 
            nsim, length(t) - 1)
Y <- cbind(rep(0, nsim), t(apply(Y, 1, cumsum)))
thr = which(Y[1, ] < 0)[1]
Y[1, thr:ncol(Y)] = 0

Z <- matrix(rnorm(mean = 0, n = nsim * (length(t) - 1), sd = sqrt(sig2/4)), 
            nsim, length(t) - 1)
Z <- cbind(rep(0, nsim), t(apply(Z, 1, cumsum)))
thr = which(Y[1, ] < 0)[1]

png("allele_freqs.png", width=600, height=300)
plot(t, X[1, ], xlab = "Time (generations)", ylab = "Allele frequency", ylim=c(0, 1), type="l",
     col="darkgreen", lwd=2)
lines(t+20, Y[1, ], col="indianred", lwd=2)
legend("topright", legend=c("beneficial mutation", "deleterious mutation"), fill=c("darkgreen", "indianred"))
dev.off()


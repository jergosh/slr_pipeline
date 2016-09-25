library(vioplot)
gamma1 <- rgamma(10000, shape=0.05, rate=1)
gamma2 <- rgamma(10000, shape=0.055, rate=0.98)
gamma3 <- rgamma(10000, shape=0.06, rate=0.95)

c(mean(gamma1), mean(gamma2), mean(gamma3))

ks.test(gamma1, gamma2)
ks.test(gamma1, gamma3)
ks.test(gamma2, gamma3)

png("3_hists.png", height=200, width=600)
par(mfrow=c(1, 3))
hist(log(gamma1), breaks=100)
hist(log(gamma2), breaks=100)
hist(gamma3, breaks=100)
par(mfrow=c(1, 1))
dev.off()

png("violin_example.png", height=400, width=480)
vioplot(gamma1, gamma2, gamma3)
dev.off()


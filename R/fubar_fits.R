##
## Curve fitting
##

## Another idea:
# Numerical derivative < -1?
# Could be done by fitting a spline so that we can evaluate at arbitrary points

test.df <- subset(fubar_table_filtered, ID == "10_1")
ord = order(test.df$alpha, decreasing=T)
test.df <- test.df[ord, ]
test.df <- test.df[1:200, ]
test.df$x <- 1:nrow(test.df)

spline <- smooth.spline(x=test.df$x, y=test.df$alpha)
spline_f <- splinefun(x=test.df$x, y=test.df$alpha, method="hyman")
which(spline_f(test.df$x, deriv=T) > -1)
plot(test.df$x, spline_f(test.df$x, deriv=T))

plot(test.df$x, test.df$alpha)

weights <- exp(seq(nrow(test.df), 1, -1) / 20)
fit <- nls(log(alpha) ~ log(b) -l*x, data=test.df[1:200, ], start=list(b=1, l=1), weights=)
lines(test.df$x[1:200], exp(predict(fit)))
f <- function(x) {  }

abline(v=max(which(cumsum(exp(predict(fit) ) < test.df$alpha[1:200]) == 1:200)))

plotFit <- function(df) {
  ord = order(df$alpha, decreasing=T)
  df <- df[ord, ]
  df <- df[1:min(nrow(df), 200), ]
  df$x <- 1:nrow(df)
  
  plot(df$x, df$alpha)
  
  weights <- c(rep(100,  50), rep(1, nrow(df)-50))
  fit <- nls(log(alpha) ~ log(b) -l*x, data=df, start=list(b=1, l=1), weights=weights)
  lines(df$x, exp(predict(fit)))
  
  last_val <- max(which(cumsum(exp(predict(fit)) < df$alpha) == 1:nrow(df)))
  abline(v=last_val)
  df[1:last_val, ]
}

pdf("fit_plots.pdf")
fubar_table_top <- ddply(fubar_table_filtered[1:200000, ], .variables=c("ID"), .fun=plotFit)
dev.off()

x <- 1:100
d <- 0.0
# x <- test.df$alpha
deriv <- numericDeriv(quote(x+d^2), "d")
lines(attr(deriv, "gradient"))


#
# Fitting gamma distributions
#
dat <- subset(fubar_table_filtered, ID == "1018_3")
dgamma_wrap <- function(x, shape, scale) {
  print(c(shape, scale))
  print((x))
  print(dgamma(x, shape, scale))
  dgamma(x, shape=shape, scale=scale)
}

fit_gamma <- function(dat) {
  print(dat$ID[1])
  x <- dat$alpha / mean(dat$alpha)
  
  hist(x, freq=F, breaks=100, main=dat$ID[1])
  m2 <- MASS::fitdistr(x, "gamma",
                       start=list(shape=2,scale=2), lower=c(0.01,0.01))
  
  lines(seq(0, 6, 0.1),
        dgamma(seq(0, 6, 0.1), shape=m2$estimate[1], scale=m2$estimate[2]),
        col="red")
}

fit_gamma(dat)
pdf("gamma_fits.pdf")
ddply(fubar_table_filtered, .variables=c("ID"), .fun=fit_gamma)
dev.off()
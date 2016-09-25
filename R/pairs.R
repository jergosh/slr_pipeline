library(scales)
library(beanplot)

qqplot(subset(pdb_master_globular, sec_simple == "Helix" & buried == TRUE & omega < 5 & is_domain == TRUE)$omega,
       subset(pdb_master_globular, sec_simple == "Coil" & buried == TRUE & omega < 5 & is_domain == TRUE)$omega)
abline(a=0,b=1)


pairs.sec <- function(df) {
  sec_simple.split <- split(df$omega, df$sec_simple)
  
  max.length <- max(unlist(lapply(sec_simple.split, length)))
  sec_simple.split <- llply(sec_simple.split, function(l) {
    c(l, rep(NA, max.length-length(l)))
  })
  
  sec_simple.df <- as.data.frame(sec_simple.split)
  
  pairs(sec_simple.df, panel=function(x, y, ...) {
    sx <- sort(x[!is.na(x)])
    sy <- sort(y[!is.na(y)])
    
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
      sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
      sy <- approx(1L:leny, sy, n = lenx)$y
    points(sx, sy, ...)
    abline(a=0, b=1)
  })
}

# Buried and non-buried, domains and non-domains
pdb_master_globular_buried <- subset(pdb_master_globular, buried == TRUE)
pairs.sec(subset(pdb_master_globular_buried, omega < 20))
pairs.sec(subset(pdb_master_globular_buried, omega < 10))

pairs.sec(subset(pdb_master_globular_buried, omega < 10 & is_domain == TRUE))
pairs.sec(subset(pdb_master_globular_buried, omega < 10 & is_domain == FALSE))


# Since it doesn't look like any variant of pairs will work, we can settle for the pairwise plots
pairs.sec.helper <- function(x, y, ...) {
  lenx <- length(x)
  leny <- length(y)
  
  if (leny < lenx) 
    x <- approx(1L:lenx, x, n = leny)$y
  if (leny > lenx) 
    y <- approx(1L:leny, y, n = lenx)$y
  
  points(x, y, ...)
}

plot.pairs <- function(df, title="") {
  opar <- par(mfrow=c(length(levels(df$sec_simple)), length(levels(df$sec_simple))), mar=rep(2.6, 4),
              oma=c(0,0,2,0))
  
  x.range <- range(df$omega)
  y.range <- range(df$omega)  
  
  for (level_x in levels(df$sec_simple)) {
    sx <- subset(df, sec_simple == level_x)
    sx.domains <- sort(subset(sx, is_domain == TRUE)$omega)
    sx.rest <- sort(subset(sx, is_domain == FALSE)$omega)
    
    for (level_y in levels(df$sec_simple)) {
      x_n <- which(levels(df$sec_simple) %in% level_x)
      y_n <- which(levels(df$sec_simple)  %in% level_y)
      
      if (x_n == y_n) {
        plot.new()
         box()
        text(0.5, 0.5, level_x)
         next
      }
      if (x_n > y_n) {
        plot.new()
        box()
        domains.p <- ks.test(sx.domains, sy.domains)$p
        rest.p <- ks.test(sx.rest, sy.rest)$p
        
        next
      }
  
      sy <- subset(df, sec_simple == level_y)
      sy.domains <- sort(subset(sy, is_domain == TRUE)$omega)
      sy.rest <- sort(subset(sy, is_domain == FALSE)$omega)      
             
      plot(NA, xlim=x.range, ylim=y.range, xlab=level_x, ylab=level_y)
      pairs.sec.helper(sx.domains, sy.domains, col=alpha("red", 0.5), cex=0.25)
      pairs.sec.helper(sx.rest, sy.rest, col=alpha("grey", 0.5), cex=0.25)
      # legend("topleft", c(paste0("Domains (", domains.p, ")"), paste0("Rest (", rest.p, ")")), fill=c("red", "grey"))
      abline(a=0, b=1)
    }
  }
  
  title(title, outer=TRUE)
  par(opar)
}

pdf("qqplot_pairs_buried_10.pdf")
plot.pairs(subset(pdb_master_globular, buried == TRUE & omega < 10), title="Buried")
dev.off()
pdf("qqplot_pairs_exposed_10.pdf")
plot.pairs(subset(pdb_master_globular, buried == FALSE & omega < 10), title="Exposed")
dev.off()

pdf("qqplot_pairs_buried_5.pdf")
plot.pairs(subset(pdb_master_globular, buried == TRUE & omega < 5), title="Buried")
dev.off()
pdf("qqplot_pairs_exposed_5.pdf")
plot.pairs(subset(pdb_master_globular, buried == FALSE & omega < 5), title="Exposed")
dev.off()

# TODO Make the paired plots separately to make sure the result is the same
for (buriedness in c("buried", "exposed")) {
  x <- 1
  subset_buried <- subset(pdb_master_globular, buried == buriedness)
    
  for (level_x in levels(subset_buried$sec_simple)) {
    subset_x <- subset(subset_buried, omega < 10 & sec_simple == level_x)

    y <- 1
    for (level_y in levels(subset_buried$sec_simple)) {
      if (x < y) {
        next
      }
      subset_y <- subset(subset_buried, omega < 10 & sec_simple == level_y)
      
      for (domain_status in c("domain", "non-domain")) {
        subset_x_domain <- subset(subset_x, is_domain == domain_status)
        subset_y_domain <- subset(subset_y, is_domain == domain_status)        
        
        if (domain_status == "domain") {
          color <- "red"
        } else {
          color <- "black"
        } 
        
        pdf(paste("qqplots_separate/qqplot", level_x, level_y,
                  buriedness, domain_status, ".pdf", sep="_"), height=11, width=10)
        qqplot(subset_x_domain$omega, subset_y_domain$omega, 
               main=paste(buriedness, domain_status), xlab=level_x, ylab=level_y, col=color)
        abline(a=0, b=1)
        dev.off() 
      }
      y <- y + 1
    }
    x <- x + 1
  }
}

# Domain and non-domain together
omega_limit <- 5
pdf(paste0("qqplots_paired_", omega_limit, ".pdf"), height=11, width=10)
for (buriedness in c("buried", "exposed")) {
  x <- 1
  subset_buried <- subset(pdb_master_globular, buried == buriedness)
  
  for (level_x in levels(subset_buried$sec_simple)) {
    subset_x <- subset(subset_buried, omega < omega_limit & sec_simple == level_x)
    
    y <- 1
    for (level_y in levels(subset_buried$sec_simple)) {
      if (x <= y) {
        next
      }
      subset_y <- subset(subset_buried, omega < omega_limit & sec_simple == level_y)
      
      
      subset_x_domain <- sort(subset(subset_x, is_domain == "domain")$omega)
      subset_y_domain <- sort(subset(subset_y, is_domain == "domain")$omega)
      l_x_domain <- length(subset_x_domain)
      l_y_domain <- length(subset_y_domain)
            
      if (l_y_domain < l_x_domain)
        subset_x_domain <- approx(1L:l_x_domain, subset_x_domain, n=l_y_domain)$y      
      if (l_x_domain < l_y_domain)
        subset_y_domain <- approx(1L:l_y_domain, subset_y_domain, n=l_x_domain)$y
            
      subset_x_nondomain <- sort(subset(subset_x, is_domain == "non-domain")$omega)
      subset_y_nondomain <- sort(subset(subset_y, is_domain == "non-domain")$omega) 
      l_x_nondomain <- length(subset_x_nondomain)
      l_y_nondomain <- length(subset_y_nondomain)
      
      if (l_y_nondomain < l_x_nondomain)
        subset_x_nondomain <- approx(1L:l_x_nondomain, subset_x_nondomain, n=l_y_nondomain)$y      
      if (l_x_nondomain < l_y_nondomain)
        subset_y_nondomain <- approx(1L:l_y_nondomain, subset_y_nondomain, n=l_x_nondomain)$y
      
      plot(subset_x_nondomain, subset_y_nondomain, xlim=c(0, omega_limit), ylim=c(0, omega_limit),
             main=paste(buriedness), xlab=level_x, ylab=level_y, col="black")
      points(subset_x_domain, subset_y_domain, 
             main=paste(buriedness, domain_status), col="blue")
      rect(1, 1, omega_limit, omega_limit, col=rgb(1, 0, 0, 0.25), border=NA)
      abline(a=0, b=1)
      p_nondomains <- ks.test(subset_x_nondomain, subset_y_nondomain)$p
      p_domains <- ks.test(subset_x_domain, subset_y_domain)$p
      legend("topleft", legend=paste(c("domain", "non-domain"), c(p_domains, p_nondomains)), fill=c("blue", "black"))      
      
      y <- y + 1
    }
    x <- x + 1
  }
}
dev.off() 

# Actually qqplots and correlations are only half of the picture -- what about beanplots?
# Need to have identical range and set them up with mfrow
# And also possibly log them (which might in turn require adding a small epsilon to the zeroes)
# Beanplots

y.range <- range(pdb_master_globular$omega)

pdb_master_globular$omega_asinh <- asinh(pdb_master_globular$omega)
pdb_master_globular$omega_eps <- pdb_master_globular$omega + 1e-16

beanplot(subset(pdb_master_globular, sec_simple == "Helix" & buried == TRUE & omega < 5)$omega,
#          subset(pdb_master_globular, sec_simple == "Beta sheet" & buried == TRUE)$omega_eps,
#          subset(pdb_master_globular, sec_simple == "Coil" & buried == TRUE)$omega_asinh,
#          subset(pdb_master_globular, sec_simple == "Turn" & buried == TRUE)$omega_asinh,
#          subset(pdb_master_globular, sec_simple == "Helix" & buried == FALSE)$omega_asinh,
#          subset(pdb_master_globular, sec_simple == "Beta sheet" & buried == FALSE)$omega_asinh,
#          subset(pdb_master_globular, sec_simple == "Coil" & buried == FALSE)$omega_asinh,
#          subset(pdb_master_globular, sec_simple == "Turn" & buried == FALSE)$omega_asinh,
#          names=c("Helix", "Beta sheet", "Coil", "Turn", "Helix", "Beta sheet", "Coil", "Turn")
# log="y",
bw="nrd0"
)
# Sec str. x buried x domain?
for (buriedness in c(TRUE, FALSE)) {
  
}

# log

# asinh
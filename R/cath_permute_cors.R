# Permute values 
sample_ids <- sample(unique(cath_gpcr$stable_id), size=10, replace=F)
sample_pairs <- combn(sample_ids, 2, simplify=F)

laply(sample_pairs, function(pair) {
  x <- subset(cath_gpcr, stable_id == pair[1])$omega
  y <- subset(cath_gpcr, stable_id == pair[2])$omega
  
  pvals <- aaply(1:10000, 1, function(n) {
    cor.test(x, sample(y, length(y), replace=F), method="spearman")$p.value
  })
  
  hist(pvals, breaks=100)  
})

hist(laply(sample_pairs, function(pair) {
  x <- subset(cath_gpcr, stable_id == pair[1])$omega
  y <- subset(cath_gpcr, stable_id == pair[2])$omega
  
  ks.test(x, y, exact=T)$p.value
}), breaks=100)

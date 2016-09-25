omega_1 <- subset(cath_all, stable_id == "ENSP00000345163")$omega
omega_2 <- subset(cath_all, stable_id == "ENSP00000318956")$omega
plot(omega_1, omega_2)
cor.test(omega_1, omega_2, method="pearson")
cor.test(omega_1, omega_2, method="spearman")
 

cath.cors.same <- ddply(cath.cors, c("id_1", "id_2"), function(r) {
  if (strsplit(id2dataset[[ r$id_1[1] ]], "_", fixed=T)[[1]][1] != strsplit(id2dataset[[ r$id_2[1] ]], "_", fixed=T)[[1]][1]) {
    return(data.frame())
  } else {
    return(r)
  }
})

pdf("domain_correlations_same.pdf", width=11, height=9)
hist(as.numeric(cath.cors.same$pval), breaks=100, main="Recent gene duplication", xlab="P-value")
dev.off()

pdf("domain_correlations_different.pdf", width=11, height=9)
hist(as.numeric(cath.cors.different$pval), breaks=100, main="Other", xlab="P-value")
dev.off()

ploop <- subset(cath.cors, cath_id == "3.40.50.300")
nrow(ploop)
pdf("ploop_hist.pdf", height=9, width=11)
hist(as.numeric(ploop$pval), breaks=100, main="CATH superfamily 3.40.50.300", xlab="P-value")
dev.off()


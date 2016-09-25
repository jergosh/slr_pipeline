domain.cors <- function(df) {
  print(paste(df$cath_id[1], nrow(df)))
  if (length(unique(df$stable_id)) < 2) {
    return(data.frame())
  }
  cath_id <- df$cath_id[1]
  res <- data.frame(cath_id=character(), id_1=character(), id_2=character(), corr=numeric(), pval=numeric(), stringsAsFactors=F)
  for (comb in combn(as.character(unique(df$stable_id)), 2, simplify=F)) {
    if (comb[1] == comb[2]) {
      print("IDs identical!")
    }
    omega_1 <- subset(df, stable_id == comb[1])$omega
    omega_2 <- subset(df, stable_id == comb[2])$omega
    
    # Make sure there is actually stuff to correlate
    if (sum((!is.na(omega_1)) & (!is.na(omega_2))) < 10) {
      res[nrow(res)+1, ] <- c(cath_id, comb[1], comb[2], NA, NA)
    } else {
      corr <- cor.test(omega_1, omega_2, use="complete.obs", method="spearman")
      res[nrow(res)+1, ] <- c(cath_id, comb[1], comb[2], corr$estimate, corr$p.value)
    }
  }
  
  res
}

cath_tmp <- read.table("data/fam_tmp.tab", sep="\t", stringsAsFactors=F, header=T)
tmp_all <- ddply(cath_tmp, .variables="cath_id", .fun=domain.cors)


hist(as.numeric(tmp_all$pval), breaks=100)
hist(as.numeric())
tmp_all$p.adjusted <- p.adjust(tmp_all$pval, method="BH")
subset_cors_sign <- subset(tmp_all, p.adjusted < 1e-05)

cors_sif <- data.frame(P1=subset_cors_sign$id_1, int="cor", P2=subset_cors_sign$id_2)
write.table(cors_sif, file="data/cath/cytoscape/1.20.58.390.sif", sep="\t", quote=F, row.names=F)
